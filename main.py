import warnings
import time
import pandas as pd
import numpy as np
from src.Log import Log
from src.VCFProcessor import VCFProcessor, DataCleaner
from src.EnsemblAPI import EnsemblAPI
from src.MSPUpdater import MSPUpdater
from src.RisComparator import RisComparator
from ydata_profiling import ProfileReport

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
        df = processor.get_dataframe()
        df = cleaner.clean_dataframe(df)
        df = api.get_api_info_from_df(df, path_json="memory.json")
        df = msp.aggiungi_colonna_msp(df)
        df = ris.compare_vcf_xlsx(df)
        df = df[df["MSP"].notna()]
        df = df[df["RIS."].notna()]
        df = df[df["RIS."] != "nan"]
        df.to_csv("Data/dataset.csv", index=False)
        log.write_log("Program finished successfully", level="SUCCESS")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
    except KeyboardInterrupt:
        log.write_log("Program interrupted by user", level="CRITICAL")
        log.write_log(f"Total time: {time.time() - start_program_time} seconds", level="DEBUG")
        
    df = pd.read_csv("data/dataset.csv")
    profile = ProfileReport(df, title="Profiling Report")
    profile.to_file("your_report.html")