import pandas as pd
from tqdm import tqdm


class RisComparator: #TODO cambiare readcsv to readexcel
    def __init__(self, brca_path, hc_path, log):
        self.brca_df = self.__load_dataframe(brca_path)
        self.hc_df = self.__load_dataframe(hc_path)
        self._full_df = pd.concat([self.brca_df, self.hc_df])
        self.__log = log

    def __load_dataframe(self, path):
        df = pd.read_excel(path)
        df = df[df["ANNO"]>2017]
        return df

    def __missing_MSP(self,vcf_df, msp_df, identificativo):
        msp_hc_brca = set(vcf_df["MSP"])
        msp = set(msp_df["MSP"])
        for x in list(msp_hc_brca - msp):
            self.__log.write_log(level="DEBUG",message=f"Nel dataframe {identificativo}, abbiamo un MSP IN più con codice: {x}")

    def compare_vcf_xlsx(self, df:pd.DataFrame):
        for msp in tqdm(set(df["MSP"])):
            found = False
            msp_df = df[df["MSP"] == msp]
            genes_ris = self._full_df[self._full_df["MSP"] == msp]
            
            if genes_ris.empty:
                self.__log.write_log(level="ERROR",message=f"MSP: {msp}, non trovato")
                continue
            
            if "NEG" in genes_ris["RIS."].values[0]:
                df.loc[df["MSP"] == msp, "RIS."] = "NEG"
            else:
                
                df.loc[df["MSP"] == msp, "RIS."] = "NEG"
                hgvs = genes_ris["VARIANTE"].values[0]
                if isinstance(hgvs, float):
                    self.__log.write_log(level="ERROR",message=f"MSP: {msp}, HGVS sbagliato: {hgvs}")
                    continue
                hgvs = hgvs.replace(" ", "").replace("\n", "")
                for _, row in msp_df.iterrows():
                    if isinstance(row["hgvsc"], float):
                        continue

                    if row["hgvsc"].split(":")[1] in hgvs:
                        found = True
                        df.loc[(df["MSP"] == msp) & (df["hgvsc"] == row["hgvsc"]), "RIS."] = genes_ris["RIS."].values[0]
                
                if not found:
                    self.__log.write_log(level="ERROR",message=f"MSP: {msp}, HGVS sbagliato: {hgvs}")
        
        return df
        