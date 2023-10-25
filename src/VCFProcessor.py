

import os
import allel
import pandas as pd
from tqdm import tqdm
import time

class VCFProcessor:
    def __init__(self, log):
        self.log = log
        self.identifier_set = set()

    @staticmethod
    def should_process_file(file_name, exclude_patterns):
        result = file_name.endswith(".vcf") and not any(x in file_name.upper() for x in exclude_patterns)
        return result

    @staticmethod
    def extract_identifier(vcf_file):
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    identifier = line.strip().split('\t')[-1]
                    if "full_variant" in vcf_file.lower() or "full-variant" in vcf_file.lower():
                        identifier = identifier.split('_')[1]
                    try:
                        identifier = identifier.split('_')[0]
                        identifier = identifier.replace("-", "/")
                    except Exception as e:
                        pass

                    return identifier
        return None

    def process_vcf_file(self, vcf_file):
        splitted_path = vcf_file.split(os.sep)
        tissue = splitted_path[-2].upper()
        identifier = self.extract_identifier(vcf_file)
        
        if identifier is None:
            self.log.write_log(f"Unable to extract identifier from {vcf_file}", level="ERROR")
            return None

        if identifier in self.identifier_set:
            self.log.write_log(f"Duplicate identifier found: {identifier}, File path: {vcf_file}", level="WARNING")
            return None

        self.identifier_set.add(identifier)
        if not identifier.upper().startswith("FC") and not splitted_path[-1].upper().startswith("FC") and not splitted_path[-2].upper().startswith("FC"):
            temp_df = allel.vcf_to_dataframe(vcf_file, fields='*', alt_number=1)
            
            if temp_df is not None:
                temp_df['name'] = identifier
                
                if identifier.upper().startswith("BRCA"):
                    temp_df['CTYPE'] = "BRCA"
                    temp_df['TISSUE'] = tissue
                
                elif identifier.upper().startswith("HC"):
                    temp_df['CTYPE'] = "HC"
                    temp_df['TISSUE'] = "GERMLINE"
                
                else:
                    if splitted_path[-1].upper().startswith("BRCA"):
                        temp_df['CTYPE'] = "BRCA"
                        temp_df['TISSUE'] = tissue
                        self.log.write_log(f"CTYPE and identifier mismatch, identifier: {identifier}, CTYPE: BRCA, File path: {vcf_file}", level="WARNING")
                    
                    elif splitted_path[-1].upper().startswith("HC"):
                        temp_df['CTYPE'] = "HC"
                        temp_df['TISSUE'] = "GERMLINE"
                        self.log.write_log(f"CTYPE and identifier mismatch, identifier: {identifier}, CTYPE: HC, File path: {vcf_file}", level="WARNING")
                    
                    else:
                        temp_df['CTYPE'] = "UNKNOWN"
                        temp_df['TISSUE'] = tissue
                        self.log.write_log(f"Unknown identifier: {identifier}, File path: {vcf_file}", level="ERROR")
                
                return temp_df
        else:
            return None

    def get_dataframe(self, path="Data/VCF", exclude_patterns=None):
        if exclude_patterns is None:
            exclude_patterns = ["POS", "NEG", "PROVA", "VEQ", "EM", "CTRL"]
        
        data_folder = os.path.join(os.getcwd(), path)
        df_list = []

        for root, dirs, files in os.walk(data_folder):
            print(f"Processing: {root}")
            self.log.write_log(f"Starting processing directory {root}", level="INFO")
            start_time = time.time()

            for file in tqdm(files, leave=False):
                if self.should_process_file(file, exclude_patterns):
                    vcf_file = os.path.join(root, file)
                    df_processed = self.process_vcf_file(vcf_file)
                    if df_processed is not None:
                        df_list.append(df_processed)

            self.log.write_log(f"Finished processing directory {root}", level="SUCCESS")
            self.log.write_log(f"Processing time: {time.time() - start_time} seconds", level="DEBUG")
            print("Done!\n")

        df = pd.concat(df_list, ignore_index=True)
        return df

    @staticmethod
    def load_data(path="Data/processed_data.csv",log=None ):
        try:
            df = pd.read_csv(path)
        except Exception as e:
            if log is not None:
                log.write_log(f"Program finished with error: {e}", level="CRITICAL")
            exit()
        return df
    
class DataCleaner:
    def __init__(self, log):
        self.log = log
    
    def drop_nan_columns(self, dataframe):
        df = dataframe.copy()
        for col in df.columns:
            if df[col].isnull().all():
                df.drop(col, axis=1, inplace=True)
        return df

    def clean_dataframe(self, dataframe):
        df = dataframe.copy()
        pre_num_cols = len(df.columns)
        df = self.drop_nan_columns(df)
        self.log.write_log(f"Dropped {pre_num_cols - len(df.columns)} NaN columns", level="DEBUG")
        self.log.write_log("Dropping FILTER columns", level="DEBUG")
        pre_num_cols = len(df.columns)
        for column in df.columns:
            if "FILTER" in column:
                df.drop(column, axis=1, inplace=True)
        self.log.write_log(f"Dropped {pre_num_cols - len(df.columns)} FILTER columns", level="DEBUG")
        return df