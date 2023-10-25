import warnings
import time
from src.Log import Log
from src.VCFProcessor import VCFProcessor, DataCleaner
from src.utils import check_integrity

warnings.filterwarnings('ignore')

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
