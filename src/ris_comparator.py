import pandas as pd
import numpy as np


class Ris_comparator: #TODO cambiare readcsv to readexcel
    def __init__(self, brca_path, hc_path, log):
        self.brca_dataframe = self.__load_dataframe(brca_path)
        self.hc_dataframe = self.__load_dataframe(hc_path)
        self.__log = log

    def __load_dataframe(self, path):
        df = pd.read_csv(path)
        df = df[df["ANNO"]>2017]
        df_pulito = df[["RIS.", "MSP"]]
        return df_pulito

    def __missing_MSP(self,hc_brca_dataframe,dataframe,identificativo="brca" or "hc"):
        lista_MSP_hc_brca = set(hc_brca_dataframe["MSP"])
        lista_MSP = set(dataframe["MSP"])
        differenza = lista_MSP_hc_brca - lista_MSP
        for x in list(differenza):
            self.__log.write_log(level="ERROR",message=f"Nel dataframe {identificativo}, abbiamo un MSP IN più con codice: {x}")
        return

    def create_new_dataframe(self, dataframe):
        self.__missing_MSP(self.brca_dataframe,dataframe,identificativo="brca")
        self.__missing_MSP(self.hc_dataframe,dataframe,identificativo="hc")

        #self.__missing_MSP(self.hc_dataframe)
        lista_ris = []
        #c=0
        for index, row in dataframe.iterrows():
            ctype, msp, ris = row["CTYPE"], row["MSP"], row["CLEANVARPAT"]
            if ctype == "BRCA":
                genes_ris = self.brca_dataframe[self.brca_dataframe["MSP"] == msp]["RIS"].values
            else:
                genes_ris = self.hc_dataframe[self.hc_dataframe["MSP"] == msp]["RIS"].values
            if genes_ris.size > 0:
                genes_ris = genes_ris[0]
                if ris == genes_ris:
                    ris_finale = genes_ris
                elif ris == "VUS":
                    ris_finale = genes_ris
                elif genes_ris == "VUS":
                    ris_finale = ris
                else:
                    ris_finale = np.nan
                    self.__log.write_log(level="ERROR",message=f"Incongruenza fra i dati dei VCF Tipo di varianza: {ctype}, MSP: {msp}")

                lista_ris.append(ris_finale)
            else:
                lista_ris.append("NOT FOUND")
                self.__log.write_log(level="ERROR",
                                     message=f"MSP: {msp}, non prensente nel BRCA e neanche nel HC")

                #dataframe = dataframe.drop(c)
            #c+=1
        dataframe["RIS"] = lista_ris
        return dataframe

