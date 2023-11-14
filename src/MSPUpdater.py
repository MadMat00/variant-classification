import pandas as pd
import re
from Log import Log

class MSPUpdater:
    def __init__(self, gmo_path, log):
        self.gmo_path = gmo_path
        self.log = log
        self.nomi_non_trovati_registrati = set()

    @staticmethod
    def pulisci_stringa(s):
        # Pulisce la stringa: rimuove spazi bianchi, converte in maiuscolo
        s = s.strip()
        s = s.replace(" ", "")
        s = s.upper()
        return s

    def rimuovi_zero_dopo_brca_hc_dataset(self, df):
        # Applica la rimozione degli zeri alla colonna "NAME" del DataFrame
        #df["NAME"] = df["NAME"].apply(lambda s: re.sub(r'BRCA0(\d)', r'BRCA\1', s) if pd.notna(s) else s)
        #df["NAME"] = df["NAME"].apply(lambda s: re.sub(r'HC0(\d)', r'HC\1', s) if pd.notna(s) else s)
        for _,row in df.iterrows():
            nome = row["NAME"]
            if nome[2] == 0 and nome[:2] == "HC":
                nome[2].pop()
            elif nome[4] == 0 and nome[:4] == "BRCA":
                nome[4].pop()
            else:
                ...
            row["NAME"] = nome
        return df

    def aggiungi_colonna_msp(self, df):
        # Carica il file GMO come DataFrame
        origine = pd.read_excel(self.gmo_path, engine="xlrd")
        origine.iloc[:, 1] = origine.iloc[:, 1].apply(self.pulisci_stringa)

        # Carica il dataset come DataFrame
        df = pd.read_csv(dataset_file)

        # Applica la rimozione degli zeri direttamente al DataFrame del dataset
        df = self.rimuovi_zero_dopo_brca_hc_dataset(df)

        # Esegue il merge tra il DataFrame in input e il DataFrame originale
        risultati = pd.merge(df, origine, left_on="NAME", right_on="Informazione addizionale: GEN-MOL1", how="left")

        # Aggiunge la colonna "MSP" al DataFrame in input
        df["MSP"] = risultati["Codice Esterno"]
        df = df.dropna(subset=["MSP"])

        # Identifica i nomi non trovati e registra un avviso nel log
        nomi_non_trovati = risultati[pd.isna(risultati["Codice Esterno"])]
        for nome_non_trovato in nomi_non_trovati["NAME"]:
            if nome_non_trovato not in self.nomi_non_trovati_registrati:
                self.log.write_log(f"Nome non trovato: {nome_non_trovato}", level="WARNING")
                self.nomi_non_trovati_registrati.add(nome_non_trovato)

        return df  # Restituisci il DataFrame con la colonna "MSP"

# Esempio di utilizzo della classe
dataset_file = "src\processed_data.csv"
origine_file = "src\Estrazione GMO 2018-2023.xls"
log = Log()  # Assicurati che la classe Log sia disponibile

# Inizializza l'oggetto MSPUpdater
msp_updater = MSPUpdater(origine_file, log)

# Leggi il dataset come DataFrame
dataset = pd.read_csv(dataset_file)

# Aggiungi la colonna "MSP" al dataset
dataset = msp_updater.aggiungi_colonna_msp(dataset)  # Assegna il risultato al DataFrame

# Salva il dataset finale con la colonna "MSP" in un nuovo file CSV eliminando le righe con "MSP" vuoto
dataset = dataset.dropna(subset=["MSP"])
dataset.to_csv("dataset_con_msp.csv", index=False)
