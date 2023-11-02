import httpx, urllib, json, os
import pandas as pd
from tqdm import tqdm
from src.Log import Log
import time

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
        self.__geneinfo_df = pd.read_csv("data/HCS_region_map.bed", sep="\t", comment="t", header=None)
        header = ["chrom", "chromStart", "chromEnd", "Exon", "GENEINFO"]
        self.__geneinfo_df.columns = header[:len(self.__geneinfo_df.columns)]
        self.__geneinfo_df = self.__geneinfo_df.groupby("GENEINFO").agg({"chrom": "mean", "chromStart": "min", "chromEnd": "max"})
        self.__geneinfo_df.reset_index(inplace=True)
        
        self.__transcript_id_df = pd.read_csv("data/GRCh38_genes_MANE_Select.txt", sep="\t", comment="t", header=None)
        header = ["Gene", "ENST", "NM"]
        self.__transcript_id_df.columns = header[:len(self.__geneinfo_df.columns)]
        self.__memory = {}
    
    def __request_data(self, url:str):
        while True:
            try:
                r = httpx.get(url, headers=self.__headers)
                if r.status_code != 200 and r.status_code != 400:
                    self.__log.write_log(f"Could not connect to Ensembl API: {r.status_code}", "ERROR")
                    continue
                elif r.status_code == 400:
                    self.__log.write_log(f"Could not connect to Ensembl API: {r.status_code}", "ERROR")
                    return None
                break
            except Exception as e:
                self.__log.write_log(f"Could not connect to Ensembl API: {e}", "ERROR")
        return r.json()
    
    def __get_variant_type(self, variant:str):
        """
        Simple heuristics to get the variant type
        
        1 182712 . A C . . . --> 1:182712-182712:1/C
        3 319780 . GA G . . . --> 3:319781-319781:1/- 
        3 319780 . GAA G . . . --> 3:319781-319782:1/- 
        19 110747 . G GT . . . --> 9:110748-110747:1/T 
        19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
        """

        variant_splitted = variant.split()

        if len(variant_splitted) == 8:
            if len(variant_splitted[3]) == 1 and len(variant_splitted[4]) == 1:
                #Simple SNP
                new_variant = "".join([variant_splitted[0], ":",
                                       variant_splitted[1], "-",
                                       variant_splitted[1], ":", "1", "/", 
                                       variant_splitted[4]])
                
            elif len(variant_splitted[3]) > len(variant_splitted[4]):
                #Deletion
                new_variant = "".join([variant_splitted[0], ":",
                                       str(int(variant_splitted[1]) + len(variant_splitted[4])), "-",
                                       str(int(variant_splitted[1]) + len(variant_splitted[3]) - len(variant_splitted[4])), ":", "1", "/", "-",
                                    ])
            elif len(variant_splitted[3]) < len(variant_splitted[4]):
                #Insertion
                #19 110747 . G GT . . . --> 9:110748-110747:1/T 
                #19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
                new_variant = "".join([variant_splitted[0], ":",
                                       str(int(variant_splitted[1]) + len(variant_splitted[3])), "-",
                                       variant_splitted[1], ":", "1", "/",
                                       variant_splitted[4][len(variant_splitted[3]):]
                                       ])
            else:
                self.__log.write_log(f"Could not parse VCF variant: {variant}", "ERROR")
                return None

            return "region", new_variant
        
    def __get_transcript_id(self, gene:str):
        transcript_id = self.__transcript_id_df[self.__transcript_id_df["Gene"] == gene]["ENST"].values[0]
        return transcript_id.split(".")[0]
        
    def __grch37_to_grch38(self, grch37:str, ref:str, alt:str, chrom:str):
        url = f"https://rest.ensembl.org/map/human/GRCh37/{grch37}/GRCh38?"
        r = self.__request_data(url)
        if r == None:
            return None
        
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
        for key in ["sift_score", "polyphen_score"]:
            try:
                #TODO: Cercare tutte le informazioni utili
                ...
            except KeyError:
                pass

    def get_api_info(self, chrom:str, pos:str, ref:str, alt:str, path_json:str):
        first_variant = f"{chrom} {pos} . {ref} {alt} . . ."
        if first_variant in self.__memory:
            transcript_consequences = self.__memory[first_variant]
            self.__get_transcript_consequences_info(transcript_consequences)
            return None
            
        result = self.__get_variant_type(first_variant)
        if result is None:
            return None
        
        grch37 = result[1].split("/")[0]
        variant_type, variant = self.__grch37_to_grch38(grch37, ref, alt, chrom)

        try:
            geneinfo = self.__get_geneinfo(chrom, pos)
        except IndexError:
            self.__log.write_log(f"Could not find geneinfo for variant {variant}", "ERROR")
            return None
        transcript_id = self.__get_transcript_id(geneinfo)
        
        url = f"http://rest.ensembl.org/vep/human/{variant_type}/{variant}?mane=1&transcript_id={transcript_id}&LoF=1&dbNSFP=MutationTaster_pred,LRT_pred,Polyphen2_HDIV_score,1000Gp3_AF&variant_class=1&Geno2MP=1&domains=1&dbscSNV=1&hgvs=1"
        r = self.__request_data(url)
        if r == None:
            return None
        
        transcript_consequences = r[0]["transcript_consequences"][0]
        self.__memory[first_variant] = r[0]
        self.__log.write_log(f"Added variant {first_variant} to memory", "DEBUG")
        with open(path_json, "w") as f:
            json.dump(self.__memory, f, indent=4)
        self.__get_transcript_consequences_info(transcript_consequences)
        
    def get_api_info_from_df(self, df:pd.DataFrame, path_json:str):
        if os.path.exists(path_json):
            self.__memory = json.load(open(path_json, "r"))
            self.__log.write_log(f"Loaded {len(self.__memory)} variants from {path_json}", "DEBUG")
    
        for row in tqdm(df.itertuples(), total=df.shape[0]):
            chrom = row[1]
            pos = row[2]
            ref = row[3]
            alt = row[4]
            self.get_api_info(chrom, pos, ref, alt, path_json)