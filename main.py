import warnings
import time
import pandas as pd
from src.Log import Log
from src.VCFProcessor import VCFProcessor, DataCleaner
from src.utils import check_integrity
from src.EnsemblAPI import EnsemblAPI

warnings.filterwarnings('ignore')

if __name__ == "__main__":
    start_program_time = time.time()
    log = Log()
    processor = VCFProcessor(log)
    cleaner = DataCleaner(log)
    api = EnsemblAPI(log)
    log.write_log("Starting program", level="INFO")
    try:
        df = pd.read_csv("data/processed_data.csv")
        api.get_api_info_from_df(df, path_json="memory.json")
        log.write_log("Program finished successfully", level="SUCCESS")
        df.to_csv("Data/processed_data.csv", index=False)
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
    except KeyboardInterrupt:
        log.write_log("Program interrupted by user", level="CRITICAL")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
