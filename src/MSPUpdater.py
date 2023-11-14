import pandas as pd
import re

class MSPUpdater:
    def __init__(self, gmo_path, log):
        self.gmo_path = gmo_path
        self.log = log
        self.nomi_non_trovati_registrati = set()  # Set per tenere traccia dei nomi non trovati già registrati

    @staticmethod
    def pulisci_stringa(s):
        # Rimuovi spazi all'inizio e alla fine della stringa
        s = s.strip()
        # Rimuovi spazi extra all'interno della stringa
        s = s.replace(" ", "")
        # Trasforma tutto in maiuscolo
        s = s.upper()
        return s

    @staticmethod
    def rimuovi_zero_dopo_brca_hc(s):
        # Rimuovi "0" dopo "BRCA" e "HC"
        s = re.sub(r'BRCA0(\d+)', r'BRCA\1', s)
        s = re.sub(r'HC0(\d+)', r'HC\1', s)
        return s

    def aggiungi_colonna_msp(self, df):
        # Leggi il file Excel di origine in un DataFrame
        origine = pd.read_excel(self.gmo_path, engine="xlrd")

        # Applica la pulizia e la rimozione degli zeri
        origine.iloc[:, 1] = origine.iloc[:, 1].apply(self.pulisci_stringa)
        origine.iloc[:, 1] = origine.iloc[:, 1].apply(self.rimuovi_zero_dopo_brca_hc)

        # Effettua una fusione (merge) tra il dataset e il file di origine
        risultati = pd.merge(df, origine, left_on="NAME", right_on="Informazione addizionale: GEN-MOL1", how="left")

        # Crea una nuova colonna "MSP" e assegna i valori dalla colonna "Codice Esterno"
        df["MSP"] = risultati["Codice Esterno"]

        # Utilizza il logger personalizzato per registrare eventuali nomi non trovati, evitando duplicati
        nomi_non_trovati = risultati[pd.isna(risultati["Codice Esterno"])]
        for nome_non_trovato in nomi_non_trovati["NAME"]:
            if nome_non_trovato not in self.nomi_non_trovati_registrati:
                self.log.write_log(f"Nome non trovato: {nome_non_trovato}", level="WARNING")
                self.nomi_non_trovati_registrati.add(nome_non_trovato)

        return df

'''
# Esempio di utilizzo della classe
dataset_file = "src\processed_data.csv"
origine_file = "src\Estrazione GMO 2018-2023.xls"
log = Log()  # Assicurati che la classe Log sia disponibile

msp_updater = MSPUpdater(origine_file, log)
dataset = pd.read_csv(dataset_file)
msp_updater.aggiungi_colonna_msp(dataset)
'''