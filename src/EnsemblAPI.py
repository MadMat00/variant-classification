from Log import Log
import httpx, urllib, json

class EnsemblAPI:
    def __init__(self, log) -> None:
        self.__log = log
        self.__headers = {
            'Host': 'rest.ensembl.org',
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.7; rv:42.0) Gecko/20100101 Firefox/42.0',
            'Accept': 'application/json, text/javascript, */*; q=0.01',
            'Accept-Language': 'en-US,en;q=0.5',
            'Referer': 'http://asia.ensembl.org/Tools/VEP?redirect=no',
            'Origin': 'http://asia.ensembl.org',
            'Connection': 'keep-alive',
        }
    
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
                new_variant = ''.join([variant_splitted[0], ':',
                                       variant_splitted[1], '-',
                                       variant_splitted[1], ':', '1', '/', 
                                       variant_splitted[4]])
                
            elif len(variant_splitted[3]) > len(variant_splitted[4]):
                # Deletion
                new_variant = ''.join([variant_splitted[0], ':',
                                       str(int(variant_splitted[1]) + len(variant_splitted[4])), '-',
                                       str(int(variant_splitted[1]) + len(variant_splitted[3]) - len(variant_splitted[4])), ':', '1', '/', '-',
                                    ])
            elif len(variant_splitted[3]) < len(variant_splitted[4]):
                # Insertion
                # 19 110747 . G GT . . . --> 9:110748-110747:1/T 
                # 19 110747 . G GTT . . . --> 19:110748-110747:1/TT 
                new_variant = ''.join([variant_splitted[0], ':',
                                       str(int(variant_splitted[1]) + len(variant_splitted[3])), '-',
                                       variant_splitted[1], ':', '1', '/',
                                       variant_splitted[4][len(variant_splitted[3]):]
                                       ])
            else:
                self.__log.write_log(f'Could not parse VCF variant: {variant}', "ERROR")

            return 'region', new_variant

    def get_api_info(self, chrom, pos, ref, alt):
        variant = f"{chrom} {pos} . {ref} {alt} . . ."
        result = self.__get_variant_type(variant)
        grch37 = result[1].split("/")[0]
        url = f"https://rest.ensembl.org/map/human/GRCh37/{grch37}/GRCh38?"
        r = httpx.get(url, headers={ "Content-Type" : "application/json"})
        r = r.json()
        
        start = r["mappings"][0]["mapped"]["start"]
        variant = f"{chrom} {start} . {ref} {alt} . . ."
        variant_type, variant = self.__get_variant_type(variant)
        variant = urllib.parse.quote(variant)
        url = f"http://rest.ensembl.org/vep/Homo_sapiens/{variant_type}/{variant}?variant_class=True"
        r = httpx.get(url, headers=self.__headers)
        decoded = r.json()
        print(json.dumps(decoded, indent=4))
        
    def get_api_info_from_df(self, df):
        pass
    

api = EnsemblAPI(log=Log())
api.get_api_info("17", "41244936", "G", "A")