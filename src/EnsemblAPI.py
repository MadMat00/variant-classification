import httpx, urllib, json, os
import pandas as pd
from tqdm import tqdm
from src.Log import Log
import numpy as np

class EnsemblAPI:
    def __init__(self, log:Log) -> None:
        self.__log = log
        self.__headers = {
            "Host": "rest.ensembl.org",
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.7; rv:42.0) Gecko/20100101 Firefox/42.0",
            "Accept": "application/json, text/javascript, */*; q=0.01",
            "Accept-Language": "en-US,en;q=0.5",
            "Referer": "http://asia.ensembl.org/Tools/VEP?redirect=no",
            "Origin": "http://asia.ensembl.org",
            "Connection": "keep-alive",
        }
        
        # GeneInfo DataFrame
        self.__geneinfo_df = pd.read_csv("data/HCS_region_map.bed", sep="\t", comment="t", header=None)
        header = ["chrom", "chromStart", "chromEnd", "Exon", "GENEINFO"]
        self.__geneinfo_df.columns = header[:len(self.__geneinfo_df.columns)]
        # Raggruppamento per geneinfo prendendo il minimo dello start e il massimo dello stop
        self.__geneinfo_df = self.__geneinfo_df.groupby("GENEINFO").agg({"chrom": "mean", "chromStart": "min", "chromEnd": "max"})
        self.__geneinfo_df.reset_index(inplace=True)
        
        # DataFrame dei identificativi dei geni
        self.__transcript_id_df = pd.read_csv("data/GRCh38_genes_MANE_Select.txt", sep="\t", comment="t", header=None)
        header = ["Gene", "ENST", "NM"]
        self.__transcript_id_df.columns = header[:len(self.__geneinfo_df.columns)]
        
        self.__memory = {}
        
        # Chiavi delle api dentro transcript_consequence
        self.__transcript_consequences_info_column = json.load(open("key.json"))["transcript_consequences"]
        
        # Chiavi delle api dentro colocated_variants
        self.__colocated_variants_info_column = json.load(open("key.json"))["colocated_variants"]
    
    def __request_data(self, url:str)->dict:
        try_count = 0
        while True:
            try_count += 1
            try:
                if try_count > 10: return None
                r = httpx.get(url, headers=self.__headers)
                if r.status_code != 200 and r.status_code != 400: # Random error
                    self.__log.write_log(f"Could not connect to Ensembl API: {r.status_code}. Try: {try_count}", "ERROR")
                    continue
                elif r.status_code == 400: # Bad Request
                    self.__log.write_log(f"Could not connect to Ensembl API: {r.status_code}. Try: {try_count}", "ERROR")
                    return None
                return r.json()
            except Exception as e: # Errore nella chiamata
                self.__log.write_log(f"Could not connect to Ensembl API: {e}", "ERROR")
            
    
    def __get_variant_type(self, variant):
        '''
        Simple heuristics to get the variant type

        1 182712 182712 A/C 1  --> 1:182712-182712:1/C
        3 319781 319781 A/- 1  --> 3:319781-319781:1/-
        19 110748 110747 -/T 1 --> 19:110748-110747:1/T
        1 160283 471362 DUP 1 --> 1:160283-471362:1/DUP
        1 1385015 1387562 DEL 1 --> 1:1385015-1387562:1/DEL  

        1 182712 . A C . . . --> 1:182712-182712:1/C
        3 319780 . GA G . . . --> 3:319781-319781:1/- 
        3 319780 . GAA G . . . --> 3:319781-319782:1/- 
        19 110747 . G GT . . . --> 9:110748-110747:1/T 
        19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
        1 160283 sv1 . <DUP> . . SVTYPE=DUP;END=471362 . --> UNABLE TO GENERATE PREVIEW!
        1 1385015 sv2 . <DEL> . . SVTYPE=DEL;END=1387562 . --> UNABLE TO GENERATE PREVIEW!
        '''

        new_variant = variant
        variant_splitted = variant.split()
        if len(variant_splitted) == 5:
            # Ensembl default

            ref_alt = variant_splitted[3].split('/')
            if len(ref_alt) == 2:
                if ref_alt[1] == '-':
                    new_variant = ''.join([
                        variant_splitted[0], 
                        ':', 
                        variant_splitted[1], 
                        '-', 
                        variant_splitted[2],
                        ':',
                        variant_splitted[4],
                        '/',
                        '-',
                        ])
                else:
                    new_variant = ''.join([
                        variant_splitted[0], 
                        ':', 
                        variant_splitted[1], 
                        '-', 
                        variant_splitted[2],
                        ':',
                        variant_splitted[4],
                        '/',
                        ref_alt[1],
                        ])
            elif variant_splitted[3] in ['DUP', 'DEL']:
                new_variant = ''.join([
                    variant_splitted[0], 
                    ':', 
                    variant_splitted[1], 
                    '-', 
                    variant_splitted[2],
                    ':',
                    variant_splitted[4],
                    '/',
                    variant_splitted[3],
                    ])
            else:
                self.__log.write_log(f"Impossibile riconoscere questa variante {variant}", "ERROR")

            return 'region', new_variant

        if len(variant_splitted) == 8:
            if len(variant_splitted[3]) == 1 and len(variant_splitted[4]) == 1:
                #Simple SNP
                new_variant = ''.join([
                    variant_splitted[0],
                    ':',
                    variant_splitted[1],
                    '-',
                    variant_splitted[1],
                    ':',
                    '1', # Not sure what this is..
                    '/',
                    variant_splitted[4],
                    ])
            elif len(variant_splitted[3]) > len(variant_splitted[4]):
                # Deletion
                new_variant = ''.join([
                    variant_splitted[0],
                    ':',
                    str(int(variant_splitted[1]) + len(variant_splitted[4])),
                    '-',
                    str(int(variant_splitted[1]) + len(variant_splitted[3]) - len(variant_splitted[4])),
                    ':',
                    '1', # Not sure what this is..
                    '/',
                    '-',
                    ])
            elif len(variant_splitted[3]) < len(variant_splitted[4]):
                # Insertion
                # 19 110747 . G GT . . . --> 9:110748-110747:1/T 
                # 19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
                new_variant = ''.join([
                    variant_splitted[0],
                    ':',
                    str(int(variant_splitted[1]) + len(variant_splitted[3])),
                    '-',
                    variant_splitted[1],
                    ':',
                    '1', # Not sure what this is..
                    '/',
                    variant_splitted[4][len(variant_splitted[3]):],
                    ])
            else:
                self.__log.write_log(f"Impossibile riconoscere questa variante {variant}", "ERROR")

            return 'region', new_variant

        if variant[0:2] == 'rs':
            # Variant identifier
            return 'id', new_variant

        if ':' in variant:
            # HGVS notation 
            return 'hgvs', new_variant

        if len(variant_splitted) > 1:
            self.__log.write_log(f"Impossibile riconoscere questa variante {variant}", "ERROR")

        return 'id', new_variant # Hoping for the best
        
    def __get_transcript_id(self, gene:str):
        if self.__transcript_id_df[self.__transcript_id_df["Gene"] == gene].shape[0] == 0:
            return None
        transcript_id = self.__transcript_id_df[self.__transcript_id_df["Gene"] == gene]["ENST"].values[0]
        return transcript_id.split(".")[0]
       
    def __grch37_to_grch38(self, grch37:str, ref:str, alt:str, chrom:str):
        url = f"https://rest.ensembl.org/map/human/GRCh37/{grch37}/GRCh38?"
        r = self.__request_data(url)
        if r == None:
            return (None, None)
        
        start = r["mappings"][0]["mapped"]["start"]
        variant = f"{chrom} {start} . {ref} {alt} . . ."
        variant_type, variant = self.__get_variant_type(variant)
        variant = urllib.parse.quote(variant)
        return (variant_type, variant)
    
    def __get_geneinfo(self, chrom:str, pos:str):
        return self.__geneinfo_df[(self.__geneinfo_df["chrom"] == int(chrom)) & 
                                  (self.__geneinfo_df["chromStart"] <= int(pos)) & 
                                  (self.__geneinfo_df["chromEnd"] >= int(pos))
                                  ]["GENEINFO"].values[0]

    def __get_transcript_consequences_info(self, transcript_consequences:dict):
        api_dict = {}
        for key in self.__transcript_consequences_info_column:
            try:
                if "count" in key:
                    api_dict[key] = len(transcript_consequences[key.replace("_count", "")][0])
                else:
                    api_dict[key] = transcript_consequences[key]
            except KeyError:
                pass
            
        return api_dict

    def __get_colocated_variants_info(self, colocated_variants:dict):
        for dictionary in colocated_variants:
            for key in self.__colocated_variants_info_column:
                if key not in dictionary :
                    break
                
                if isinstance(dictionary, dict):
                    colocated_variants = dictionary
        
        api_dict = {}
        for key in self.__colocated_variants_info_column:
            try:
                if "count" in key:
                    api_dict[key] = len(colocated_variants[key.replace("_count", "")])
                else:
                    api_dict[key] = colocated_variants[key]
            except KeyError:
                pass
            except TypeError:
                try:
                    if "count" in key:
                        api_dict[key] = len(colocated_variants[0][key.replace("_count", "")])
                    else:
                        api_dict[key] = colocated_variants[0][key]
                except KeyError:
                    pass
            
        return api_dict
    
    def __get_clean_clin_sig_allele(self, clin_sig_allele, alt):
        weight = {
            "uncertain_significance": 0,
            "conflicting_interpretations_of_pathogenicity": 0,
            "not_provided": 0,
            "benign/likely_benign": 0.8,
            "benign": 1,
            "likely_benign": 0.5,
            "pathogenic": -1,
            "likely_pathogenic": -0.5,
            "pathogenic/likely_pathogenic": -0.8,
            "other": 0
        }
        
        
        allele_sig = []
        value = 0
        for allele in clin_sig_allele.split(";"):
            if allele[0] == alt[0]:
                allele_sig.append(allele)
        
        if len(allele_sig) == 1:
            value = weight[allele_sig[0].split(":")[1]]
        elif len(allele_sig) == 0:
            return None

        
        for allele in allele_sig:
            value += weight[allele.split(":")[1]]
            
        value /= len(allele_sig)
        
        if value >= 0.3: return "NEG"
        elif value <= -0.3: return "POS"
        else: return "VUS"
    
    def get_api_info(self, chrom:str, pos:str, ref:str, alt:str, path_json:str):
        first_variant = f"{chrom} {pos} . {ref} {alt} . . ."
        
        # Controllo in memoria
        if first_variant in self.__memory:
            api = self.__memory[first_variant]
            
            try:
                api_dict = self.__get_transcript_consequences_info(api["transcript_consequences"][0])
                api_dict = {**self.__get_colocated_variants_info(api["colocated_variants"]), **api_dict}
            except KeyError:
                return None
            
            # Pulizia del clin_sig_allele
            try:
                api_dict["clin_sig_allele"] = self.__get_clean_clin_sig_allele(api_dict["clin_sig_allele"], alt)
            except KeyError:
                pass

            for key in ["most_severe_consequence", "variant_class"]:
                try:
                    api_dict["variant_class"] = api[key]
                except KeyError:
                    pass
            
            api_dict["seq_region_name"] = api["seq_region_name"]
            return api_dict
        
        # Valore per la richiesta api di VEP
        result = self.__get_variant_type(first_variant)
        if result is None:
            self.__log.write_log(f"Risultato del calcolo per la chiamata VEP è {result} quindi non è possibile fare la chiamata", level="ERROR")
            return None
        
        # Conversione GRCh37 a GRCh38
        grch37 = result[1].split("/")[0]
        variant_type, variant = self.__grch37_to_grch38(grch37, ref, alt, chrom)

        # Ricerca delle GENEINFO
        try:
            geneinfo = self.__get_geneinfo(chrom, pos)
        except IndexError:
            self.__log.write_log(f"Could not find geneinfo for variant {variant}", "ERROR")
            return None
        
        # Chiamata VEP per informazioni utili
        transcript_id = self.__get_transcript_id(geneinfo)
        if transcript_id == None:
            return None
        url = f"http://rest.ensembl.org/vep/human/{variant_type}/{variant}?mane=1&transcript_id={transcript_id}&LoF=1&dbNSFP=ALL&variant_class=1&Geno2MP=1&domains=1&dbscSNV=1&hgvs=1"
        r = self.__request_data(url)
        if r == None:
            return None
        
        # Salvataggio in memoria
        self.__memory[first_variant] = r[0]
        self.__log.write_log(f"Added variant {first_variant} to memory", "DEBUG")
        with open(path_json, "w") as f:
            json.dump(self.__memory, f, indent=4)
        
        # Ricerca dei transcript_consequences e colocated_variants
        try:
            transcript_consequences = r[0]["transcript_consequences"][0]
            colocated_variants = r[0]["colocated_variants"]
        except KeyError:
            return None
         
        
        
        
        # Estrazione delle informazioni utili
        api_dict = self.__get_transcript_consequences_info(transcript_consequences)
        api_dict = {**self.__get_colocated_variants_info(colocated_variants), **api_dict}
        
        # Pulizia del clin_sig_allele
        try:
            api_dict["clin_sig_allele"] = self.__get_clean_clin_sig_allele(api_dict["clin_sig_allele"], alt)
        except KeyError:
            pass
        
        # Estrazione del variant_class
        for key in ["most_severe_consequence", "variant_class"]:
            try:
                api_dict["variant_class"] = r[0][key]
            except KeyError:
                pass
        api_dict["seq_region_name"] = r[0]["seq_region_name"]
        return api_dict
        
    def get_api_info_from_df(self, df:pd.DataFrame, path_json:str):
        # Creazione delle colonne
        df[self.__transcript_consequences_info_column] = np.nan
        df[self.__colocated_variants_info_column] = np.nan
        
        # Caricamento della memoria
        if os.path.exists(path_json):
            self.__memory = json.load(open(path_json, "r"))
            self.__log.write_log(f"Loaded {len(self.__memory)} variants from {path_json}", "DEBUG")
    
        for index, row in tqdm(df.iterrows(), total=df.shape[0]):
            chrom = row["CHROM"]
            pos = row["POS"]
            ref = row["REF"]
            alt = row["ALT"]
            api_dict = self.get_api_info(chrom, pos, ref, alt, path_json)
            
            if api_dict is None:
                continue
            
            # Assegnazione dei valori nelle colonne
            for key in api_dict:
                df.at[index, key] = api_dict[key]
            
            # Salvataggio dei dati in un CSV
            if index % 1000 == 0 and index > 0:
                df.head(index).to_csv("data/data_vep.csv", index=False)
                self.__log.write_log(f"Saved the first {index} in '{path_json}'", level="SUCCESS")
                
        df.to_csv("data/data_vep.csv", index=False)
        return df