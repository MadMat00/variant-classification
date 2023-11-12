import pandas as pd
import re

class MSPUpdater:
    def __init__(self, dataset_file, origine_file, log):
        self.dataset_file = dataset_file
        self.origine_file = origine_file
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

    def aggiungi_colonna_msp(self):
        # Leggi il dataset Pandas
        dataset = pd.read_csv(self.dataset_file)

        # Leggi il file Excel di origine in un DataFrame
        origine = pd.read_excel(self.origine_file, engine="xlrd")

        # Applica la pulizia solo alla colonna [:, 1] (seconda colonna)
        colonna = origine.iloc[:, 1].apply(self.pulisci_stringa)
        
        origine.iloc[:, 1] = origine.iloc[:, 1].apply(self.pulisci_stringa)

        # Effettua una fusione (merge) tra il dataset e il file di origine
        risultati = pd.merge(dataset, origine, left_on="NAME", right_on="Informazione addizionale: GEN-MOL1", how="left")

        # Crea una nuova colonna "MSP" e assegna i valori dalla colonna "Codice Esterno"
        dataset["MSP"] = risultati["Codice Esterno"]

        # Utilizza il logger personalizzato per registrare eventuali nomi non trovati, evitando duplicati
        nomi_non_trovati = risultati[pd.isna(risultati["Codice Esterno"])]
        for nome_non_trovato in nomi_non_trovati["NAME"]:
            if nome_non_trovato not in self.nomi_non_trovati_registrati:
                self.log.write_log(f"Nome non trovato: {nome_non_trovato}", level="WARNING")
                self.nomi_non_trovati_registrati.add(nome_non_trovato)

        return dataset


'''
# Esempio di utilizzo della classe
dataset_file = "processed_data.csv"
origine_file = "Estrazione GMO 2018-2023.xls"
log = Log()  # Inizializza il logger personalizzato

msp_updater = MSPUpdater(dataset_file, origine_file, log)
msp_updater.aggiungi_colonna_msp()
'''