from Log import Log
import httpx, urllib, json
import pandas as pd
from tqdm import tqdm

class EnsemblAPI:
    def __init__(self, log) -> None:
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
        
        self.__transcript_id_df = pd.read_csv("data/GRCh38_genes_MANE_Select.txt", sep="\t", comment="t", header=None)
        header = ["Gene", "ENST", "NM"]
        self.__transcript_id_df.columns = header[:len(self.__geneinfo_df.columns)]
    
    def __get_variant_type(self, variant):
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
                # Deletion
                new_variant = "".join([variant_splitted[0], ":",
                                       str(int(variant_splitted[1]) + len(variant_splitted[4])), "-",
                                       str(int(variant_splitted[1]) + len(variant_splitted[3]) - len(variant_splitted[4])), ":", "1", "/", "-",
                                    ])
            elif len(variant_splitted[3]) < len(variant_splitted[4]):
                # Insertion
                # 19 110747 . G GT . . . --> 9:110748-110747:1/T 
                # 19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
                new_variant = "".join([variant_splitted[0], ":",
                                       str(int(variant_splitted[1]) + len(variant_splitted[3])), "-",
                                       variant_splitted[1], ":", "1", "/",
                                       variant_splitted[4][len(variant_splitted[3]):]
                                       ])
            else:
                self.__log.write_log(f"Could not parse VCF variant: {variant}", "ERROR")

            return "region", new_variant
        
    def __get_transcript_id(self, gene):
        transcript_id = self.__transcript_id_df[self.__transcript_id_df["Gene"] == gene]["ENST"].values[0]
        return transcript_id.split(".")[0]
        
    def __grch37_to_grch38(self, grch37, ref, alt, chrom):
        url = f"https://rest.ensembl.org/map/human/GRCh37/{grch37}/GRCh38?"
        r = httpx.get(url, headers={ "Content-Type" : "application/json"})
        r = r.json()
        
        start = r["mappings"][0]["mapped"]["start"]
        variant = f"{chrom} {start} . {ref} {alt} . . ."
        variant_type, variant = self.__get_variant_type(variant)
        variant = urllib.parse.quote(variant)
        return (variant_type, variant)
    
    def __get_geneinfo(self, chrom, pos):
        return self.__geneinfo_df[(self.__geneinfo_df["chrom"] == int(chrom)) & 
                                  (self.__geneinfo_df["chromStart"] <= int(pos)) & 
                                  (self.__geneinfo_df["chromEnd"] >= int(pos))
                                  ]["GENEINFO"].values[0]

    def get_api_info(self, chrom, pos, ref, alt):
        variant = f"{chrom} {pos} . {ref} {alt} . . ."
        result = self.__get_variant_type(variant)
        grch37 = result[1].split("/")[0]
        
        variant_type, variant = self.__grch37_to_grch38(grch37, ref, alt, chrom)
        geneinfo = self.__get_geneinfo(chrom, pos)
        transcript_id = self.__get_transcript_id(geneinfo)
        
        url = f"http://rest.ensembl.org/vep/human/{variant_type}/{variant}?mane=1&transcript_id={transcript_id}&LoF=1&dbNSFP=MutationTaster_pred,LRT_pred,Polyphen2_HDIV_score,1000Gp3_AF&variant_class=1&Geno2MP=1&domains=1&dbscSNV=1&hgvs=1"
        r = httpx.get(url, headers=self.__headers)
        r = r.json()
        hgvsc = r[0]["transcript_consequences"][0]["hgvsc"]
        
        url = f"http://rest.ensembl.org/vep/human/hgvs/{hgvsc}?mane=1&transcript_id={transcript_id}&LoF=1&&dbNSFP=ALL&variant_class=1&Geno2MP=1&domains=1&dbscSNV=1&transcript_match=1"
        r = httpx.get(url, headers={ "Content-Type" : "application/json"})
        r = r.json()[0]
        transcript_consequences = r["transcript_consequences"][0]
        try:
            ...
            #TODO Trovare le informazioni giuste
        except KeyError:
            pass
        
    def get_api_info_from_df(self, df):
        for row in tqdm(df.itertuples(), total=df.shape[0]):
            chrom = row[1].replace("chr", "")
            pos = row[2]
            ref = row[4]
            alt = row[5]
            self.get_api_info(chrom, pos, ref, alt)
    

api = EnsemblAPI(log=Log())
api.get_api_info_from_df(pd.read_csv("data/processed_data.csv"))