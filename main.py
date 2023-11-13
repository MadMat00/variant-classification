import warnings
import time
import pandas as pd
from src.Log import Log
from src.VCFProcessor import VCFProcessor, DataCleaner
from src.utils import check_integrity
from src.EnsemblAPI import EnsemblAPI
from src.MSPUpdater import MSPUpdater
from src.RisComparator import RisComparator

warnings.filterwarnings('ignore')

if __name__ == "__main__":
    start_program_time = time.time()
    log = Log()
    processor = VCFProcessor(log)
    cleaner = DataCleaner(log)
    api = EnsemblAPI(log)
    msp = MSPUpdater("data/Estrazione GMO 2018-2023.xls", log)
    ris = RisComparator("data/BRCA_completo_nuovo.xlsx", "data/Database HC.xlsx", log)
    log.write_log("Starting program", level="INFO")
    try:
        df = pd.read_csv("data/processed_data.csv")
        df = api.get_api_info_from_df(df, path_json="memory.json")
        log.write_log("Program finished successfully", level="SUCCESS")
        df = msp.aggiungi_colonna_msp(df)
        df = ris.compare_vcf_xlsx(df)
        df.to_csv("Data/msp.csv", index=False)
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
    except KeyboardInterrupt:
        log.write_log("Program interrupted by user", level="CRITICAL")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
