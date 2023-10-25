from loguru import logger
import pendulum
import os
import allel
import pandas as pd
import os
from tqdm import tqdm
import warnings
import time

warnings.filterwarnings('ignore')

def check_integrity(dataframe,log):
    df = dataframe.copy()
    for col in df.columns:
        assert not df[col].isnull().all()
    log.write_log("Dataframe integrity check passed", level="SUCCESS")

class Log:
    def __init__(self, save_file:bool = True, dir_path:str = "logs"):
        self.save_file = save_file
        self.dir_path = dir_path
        self.__level_list = ["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]
        self.__date = pendulum.now()
        self.__add_log_file(self.dir_path, level="DEBUG")
        logger.remove(0)
        
    def __add_log_file(self, file_name:str, level:str = "INFO"):
        if self.save_file:
            if not os.path.exists(self.dir_path):
                os.mkdir(self.dir_path)
            file_name = f"{self.__date.year}-{self.__date.month}-{self.__date.day}-{self.__date.hour}-{self.__date.minute}-{self.__date.second}.log"
            logger.add(f"{self.dir_path}/{file_name}", format="{time} {level} {message}", level=level)

    def write_log(self, message:str, level:str = "INFO"):
        if level not in self.__level_list:
            raise ValueError(f"Invalid level: {level}. Level must be one of the following values: {self.__level_list}")
        logger.log(level, message)
        
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
    def load_data(path="Data/processed_data.csv" ):
        try:
            df = pd.read_csv(path)
        except Exception as e:
            log.write_log(f"Program finished with error: {e}", level="CRITICAL")
            exit()
        return df
    
if __name__ == "__main__":
    start_program_time = time.time()
    log = Log()
    processor = VCFProcessor(log)
    cleaner = DataCleaner(log)
    log.write_log("Starting program", level="INFO")
    try:
        df = processor.get_dataframe()
        df = cleaner.clean_dataframe(df)
        check_integrity(df,log=log)
        log.write_log("Program finished successfully", level="SUCCESS")
        df.to_csv("Data/processed_data.csv", index=False)
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
    except KeyboardInterrupt:
        log.write_log("Program interrupted by user", level="CRITICAL")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
    except Exception as e:
        log.write_log(f"Program finished with error: {e}", level="CRITICAL")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
