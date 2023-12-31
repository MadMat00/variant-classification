

import os
import allel
import pandas as pd
from tqdm import tqdm
import time
import numpy as np

class VCFProcessor:
    def __init__(self, log):
        self.log = log
        self.identifier_set = set()
        #self.cleaner = DataCleaner(log)
        self.duplicated_files = []

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
                        if (identifier.startswith("BRCA") and identifier[4] == "/") or (identifier.startswith("HC") and identifier[2] == "/"):
                            identifier = identifier.replace("/", "", 1)

                        sub_identifier = identifier.split("/")
                        sub_identifier[1] = sub_identifier[1][:2]
                        identifier = sub_identifier[0] + "/" + sub_identifier[1]

                    except Exception as e:
                        pass

                    return identifier
        return None

    def process_vcf_file(self, vcf_file):
        splitted_path = vcf_file.split(os.sep)
        origin_folder = splitted_path[-2].upper()
        identifier = self.extract_identifier(vcf_file)
        
        if identifier is None:
            if self.log is not None: self.log.write_log(f"Unable to extract identifier from {vcf_file}", level="ERROR")
            return None

        if identifier in self.identifier_set:
            if self.log is not None: self.log.write_log(f"Duplicate identifier found: {identifier}, File path: {vcf_file}", level="WARNING")
            self.duplicated_files.append(vcf_file)
            return None

        self.identifier_set.add(identifier)
        if not identifier.upper().startswith("FC") and not splitted_path[-1].upper().startswith("FC") and not splitted_path[-2].upper().startswith("FC"):
            temp_df = allel.vcf_to_dataframe(vcf_file, fields='*', alt_number=1)
            if temp_df is not None:
                temp_df['NAME'] = identifier
                if origin_folder == "HC" or origin_folder == "GERMLINE":
                    temp_df["TISSUE"] = "GERMLINE"
                elif origin_folder == "SOMATIC":
                    temp_df["TISSUE"] = "SOMATIC"
                else:
                    temp_df["TISSUE"] = "UNKNOWN"

                if "BRCA" in splitted_path[-1].upper() and "HC" not in splitted_path[-1].upper():
                    temp_df["CTYPE"] = "BRCA"

                elif "HC" in splitted_path[-1].upper() and "BRCA" not in splitted_path[-1].upper():
                    temp_df["CTYPE"] = "HC"
                elif "HC" in splitted_path[-1].upper() and "BRCA" in splitted_path[-1].upper():
                    if splitted_path[-1].upper().index("BRCA") < splitted_path[-1].upper().index("HC"):
                        temp_df["CTYPE"] = "BRCA"
                        if self.log is not None: self.log.write_log(f"found both so choosing BRCA {vcf_file}", level="DEBUG")
                    else:
                        temp_df["CTYPE"] = "HC"
                        if self.log is not None: self.log.write_log(f"found both so choosing HC {vcf_file}", level="DEBUG")
                else:
                    temp_df["CTYPE"] = np.nan
                    if self.log is not None: self.log.write_log(f"Unable to extract CTYPE from {vcf_file}", level="WARNING")
                
                with open(vcf_file, 'r') as f:
                    index = 0
                    for line in f.readlines():
                        if not line.startswith("#"):
                            try:
                                gt = line.split('GT:')[1].split('/')
                                gt1 = gt[0][-1]
                                gt2 = gt[1][0]
                                temp_df.loc[index, "GT"] = f"{gt1}/{gt2}"
                                index += 1
                            except Exception as e:
                                if self.log is not None: self.log.write_log(f"Unable to extract GT from {vcf_file}", level="WARNING")
                                temp_df.loc[index, "GT"] = np.nan
                return temp_df
        else:
            return None

    def get_dataframe(self, path="Data/VCF", exclude_patterns=None):
        if exclude_patterns is None:
            exclude_patterns = ["POS", "NEG", "PROVA", "VEQ", "EM", "CTRL","DB","DM"]
        
        data_folder = os.path.join(os.getcwd(), path)
        df_list = []

        for root, dirs, files in os.walk(data_folder):
            print(f"Processing: {root}")
            if self.log is not None: self.log.write_log(f"Starting processing directory {root}", level="INFO")
            start_time = time.time()

            for file in tqdm(files, leave=False):
                if self.should_process_file(file, exclude_patterns):
                    vcf_file = os.path.join(root, file)
                    df_processed = self.process_vcf_file(vcf_file)
                    if df_processed is not None:
                        df_list.append(df_processed)

            if self.log is not None: self.log.write_log(f"Finished processing directory {root}", level="SUCCESS")
            if self.log is not None: self.log.write_log(f"Processing time: {time.time() - start_time} seconds", level="DEBUG")
            print("Done!\n")

        df = pd.concat(df_list, ignore_index=True)
        df["GENEINFO"] = df["GENEINFO"].apply(lambda x: x.split(":")[0] if pd.notnull(x) else x)
        df.rename(columns={"CLINVARPAT": "RIS"}, inplace=True)
        if self.log is not None: self.log.write_log(f"Found {len(self.duplicated_files)} duplicates", level="WARNING")
        self.log.write_log(f"Dataset grezzo: {df.shape}", level="SUCCESS")
        cleaner = DataCleaner(self.log)
        cleaned_df = cleaner.clean_dataframe(df, columns_to_keep=["CHROM", "POS", "REF", "ALT", "AF", "GENEINFO", "NAME", "TISSUE", "CTYPE","GT"], duplicated_files=self.duplicated_files)
        self.log.write_log(f"Dataset pulito: {cleaned_df.shape}", level="SUCCESS")
        return cleaned_df


    def load_data(self,path="Data/processed_data.csv",log=None ):
        try:
            df = pd.read_csv(path)
            self.log.write_log(f"Data loaded form {path}: {df.shape}", level="INFO")
        except Exception as e:
            if log is not None:
                log.write_log(f"Program finished with error: {e}", level="CRITICAL")
            exit()
        return df
    
class DataCleaner:
    def __init__(self, log):
        self.log = log
        self.processor = VCFProcessor(None)
    
    def drop_nan_columns(self, dataframe):
        df = dataframe.copy()
        for col in df.columns:
            if df[col].isnull().all():
                df.drop(col, axis=1, inplace=True)
        return df

    def clean_dataframe(self, dataframe, columns_to_keep=["CHROM", "POS", "REF", "ALT", "AF", "GENEINFO", "NAME", "TISSUE", "CTYPE","GT"], duplicated_files=[]):

        df = dataframe.copy()
        pre_num_cols = len(df.columns)
        self.log.write_log(f"Going trough {len(duplicated_files)} files", level="DEBUG")
        for file in duplicated_files:
            old_df = df[df["NAME"] == self.processor.extract_identifier(file)]
            new_df = self.processor.process_vcf_file(file)

            if new_df is not None and new_df.shape[0] > 0 and new_df.shape[0] > old_df.shape[0]: 
                df = df[df["NAME"] != self.processor.extract_identifier(file)]
                df = pd.concat([df, new_df], ignore_index=True)
                if self.log is not None: self.log.write_log(f"Replaced {self.processor.extract_identifier(file)} with {file}", level="DEBUG")

        df = self.drop_nan_columns(df)
        if self.log is not None: self.log.write_log(f"Dropped {pre_num_cols - len(df.columns)} NaN columns", level="DEBUG")
        pre_num_cols = len(df.columns)
        if self.log is not None: self.log.write_log("Dropping FILTER columns", level="DEBUG")
        for column in df.columns:
            if "FILTER" in column or column not in columns_to_keep:
                df.drop(column, axis=1, inplace=True)

        df["CHROM"] = df["CHROM"].apply(lambda x: x.replace("chr", "") if pd.notnull(x) else x)
        
        if self.log is not None: self.log.write_log(f"Dropped {pre_num_cols - len(df.columns)} columns", level="DEBUG")
        return df.astype(str)