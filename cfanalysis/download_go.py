"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksandra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
Changes: Patryk Jarnot (2026)
"""

GO_ANNOTATIONS_FILE = "go_annotations.json"
GO_NAMES_FILE ="go_names.csv"
GO_MAX_PATH_FILE ="go_max_path.csv"


import logging
import os
import sys
import typing
import json
from pathlib import Path
from optparse import OptionParser
from cfanalysis.src.utils import http_get

from cfanalysis.src.cache import Cache


def fill_names(go_ids, cache_db: Cache, save_file="/tmp/tmp_go.csv"):
    for e, go_id in enumerate(list(go_ids)):
        URL = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=%s"
        print(URL, e, len(go_ids))
        value = http_get(URL, [go_id], "search", cache_db=cache_db)
        for result in json.loads(value)["results"]:
            if result["id"] == go_id:
                go_name = result['name']
                aspect = result["aspect"]
                print(go_id, (go_name, aspect))
                with open(save_file, "a") as f:
                    f.write(f"{go_id}\t{go_name}\t{aspect}\n")
                break


def get_GO(
        protein_list: iter,
        exclude: list,
        #save_go_file: str,
        aspect: str,
        #lack_goes: str,
        cache_db: Cache,
) -> (typing.Dict, set):
    aspect_dict = dict(F="molecular_function",
                       P="biological_process",
                       C="cellular_component")
    result = {}
    all_go = set()
    e = 0
    number_seq = "?"
    
    #protein_list = [i for i in protein_list]
    begining = len(protein_list)
    protein_go_dict = {}
    while protein_list:
        protein_run = protein_list[:1]
        print(protein_run[0])
        protein_list = protein_list[1:]
        tries = 0
        url = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=%s"
        header = dict(Accept='text/tsv')
        print(protein_run)
        text = http_get(url, [','.join(protein_run)], "geneproductid", header, cache_db=cache_db)
        logging.info(f"GO info downloaded for {protein_run} from {url} left {e}/{number_seq}")
        print(url, f"seq_no={e}", f"text={text}", f"tries={tries}")
        e += 1
        if text is None:
            print("err")
            break
        if not text.strip():
            print("lack of content", protein_run)
        for line in text.split("\n"):
            new_line = line.split("\t")
            if not line.startswith("GENE PRODUCT DB") and line.strip():
                if new_line[1] in protein_run:
                    protein_go = new_line[4]
                    aspect_go = new_line[5]
                    if aspect_go == aspect or aspect_go == aspect_dict.get(aspect):
                        annotation_type_go = new_line[7]
                        if annotation_type_go not in exclude:
                            all_go.add(protein_go)
                            if new_line[1] not in result:
                                result[new_line[1]] = [protein_go]
                            else:
                                result[new_line[1]].append(protein_go)
        for protein_acc in protein_run:
            if protein_acc in result.keys():
                #TO CLEAN
                # logging.info(f"Save go in {save_go_file} data:{str(set(result[protein_acc]))} ")
                #save_go(save_go_file, {protein_acc: set(result[protein_acc])}, mode="a")
                protein_go_dict[protein_acc] = result[protein_acc]
        #    else:
        #        with open(lack_goes, "a") as f:
        #            f.write(protein_acc + "\n")
            if e % 1000 == 0:
                logging.info(f"GO info taken from file for {protein_acc} left {e}/{number_seq}")
            if not result.get(protein_acc):
                result[protein_acc] = []
        print(len([1 for i, j in result.items() if len(j) > 0]), len(protein_list), begining,
              begining - len(protein_list))
    return result, all_go, protein_go_dict


def check_aspect(go, all_go, aspect, cache_db: Cache):
    aspect_dict = dict(F="molecular_function",
                       P="biological_process",
                       C="cellular_component")
    # print(go, all_go)
    if go is None:
        return False
    if go in all_go:
        return True
    else:
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/"
        text = http_get(url, [go], "terms", cache_db=cache_db)
        request_json = json.loads(text)
        if request_json.get("results", {}):
            aspect_go = [i for i in request_json["results"] if i["id"] == go][0]["aspect"]
            if aspect_go == aspect or aspect_go == aspect_dict.get(aspect):
                return True
    return False


def get_ancestors(
        go_list: set,
        #save_go_file: str,
        ancestors_old: dict,
        all_go: set,
        aspect: str,
        cache_db: Cache,
) -> (dict, set):
    ancestors = {}
    number_seq = len(go_list)
    for e, go in enumerate(go_list):
        if go not in ancestors_old:
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/ancestors?relations=is_a%%2Cpart_of%%2Coccurs_in%%2Cregulates"
            print(url)
            print(f"GO ancestor info downloaded for {go} from {url} left {e+1}/{number_seq}")
            tries = 0
            success = False
            while tries < 10 and not success:
                try:
                    #request = requests.get(url, timeout=10)
                    print(f"go: {go}")
                    text = http_get(url, [go.replace(':', '%3A')], "ancestors", cache_db=cache_db)
                    if text is not None:
                        request_json = json.loads(text)
                        if request_json.get("results", {}):
                            ancestors[go] = [i.get("ancestors") for i in request_json.get("results", {}) if
                                             i["id"] == go and i.get("ancestors") is not None]
                            ancestors[go] = [i for sublist in ancestors[go] for i in sublist if
                                             i != go and check_aspect(i, all_go, aspect, cache_db)]

                            all_go = all_go.union(set(ancestors[go]))
                            #save_go(save_go_file, {go: ancestors[go]}, "a")
                            success = True
                    else:
                        tries += 1
                except Exception as err:
                    print(err)
                    tries += 1
                    raise err

            if not success:
                logging.info(f"Lack of GO ancestor info for {go} left {e}/{number_seq}")
                ancestors[go] = []
    return ancestors, all_go


def save_go(file: str,
            go_set: dict,
            mode="a"
            ) -> None:
    stored_exception = None
    with open(file, mode) as f:
        for protein, goes in go_set.items():
            for go in goes:
                try:
                    f.write(f"{protein}\t{go}\n")
                except KeyboardInterrupt:
                    stored_exception = sys.exc_info()
    if stored_exception:
        raise (stored_exception[0], stored_exception[1], stored_exception[2])


def read_mapped_file(file, sign="\t"):
    result = {}
    if not os.path.exists(file):
        return result
    with open(file) as f:
        for line in f:
            if line.strip():
                line = line.split(sign)
                if line[0] not in result.keys():
                    result[line[0]] = [line[1].strip()]
                else:
                    if line[1].strip() not in result[line[0]]:
                        result[line[0]].append(line[1].strip())
    return result



def get_proteins(input_file_path):
    
    #for reading uniprot fasta file
    # with open(input_file_path) as f:
    #    for l in f.readlines():
    #        if l.startswith(">"):
    #            yield l.split("|")[1]

    proteins = []

    with open(input_file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip() != "":
                proteins.append(line.strip())
    return proteins

def get_max_path(child: str, main_GO: str, cache_db):
    url_path = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/paths/%s/"

    text = http_get(url_path, [child, main_GO], "paths", cache_db=cache_db)
    if text is not None:
        max_path_len = 0
        for result_path in json.loads(text)["results"]:
            if len(result_path) > max_path_len:
                max_path_len = len(result_path)
                return max_path_len
        return max_path_len


def get_paths(go: iter, path_path: str, aspect: str, cache_db: Cache):
    aspect_dict = dict(F="GO:0003674",
                       molecular_function="GO:0003674",
                       P="GO:0008150",
                       biological_process="GO:0008150",
                       C="GO:0005575",
                       cellular_component="GO:0005575")
    with open(path_path, "w") as f:
        for go_id in go:
            path_len = get_max_path(go_id, aspect_dict[aspect], cache_db)
            f.write(f"{go_id}\t{path_len}\n")

def add_ancestors(
        ancestors: dict,
        go_protein: dict
) -> dict:
    for e, protein in enumerate(go_protein.keys()):
        if e % 1000 == 0:
            logging.info(f"Add ancestors to {protein} {e}/{len(go_protein.keys())}")
        for old_go, anc in ancestors.items():
            if old_go in go_protein[protein]:
                go_protein[protein] += anc
                go_protein[protein] = list(set(go_protein[protein]))
    return go_protein


#old prepare_data function from go_analyse.py
def crate_annotation_file(all_go, ancestors, ouput_annotation_file):
    
    logging.info(f"Add ancestors info to protein GO")    
    all_go = add_ancestors(ancestors, all_go)


    with open(ouput_annotation_file, 'w', encoding='utf-8') as f:
        json.dump(all_go, f, indent=4, ensure_ascii=False)


def prepare_folders(input_file, exclude_IEA, ouptput_dir):

    #check if input file with proteins exists
    if not os.path.isfile(input_file):
        sys.exit("Input file with protein IDs does not exists. Exiting...")

    #create output directory    
    os.makedirs(ouptput_dir, exist_ok=True)

    #old files removal and creation of new empty files

    #go_max_path_file_path = os.path.join(ouptput_dir, GO_MAX_PATH_FILE)
    #if os.path.exists(go_max_path_file_path):
    #    os.remove(go_max_path_file_path)
    #path = Path(go_max_path_file_path)
    #path.touch()  

    go_names_file_path = os.path.join(ouptput_dir, GO_NAMES_FILE)
    if os.path.exists(go_names_file_path):
        os.remove(go_names_file_path)
    path = Path(go_names_file_path)
    path.touch()
    

    #go_proteins_file_path = os.path.join(ouptput_dir, GO_PROTEINS_FILE) 
    #if os.path.exists(go_proteins_file_path):
    #    os.remove(go_proteins_file_path)
    #path = Path(go_proteins_file_path)
    #path.touch()
    
    go_annotations_file_path = os.path.join(ouptput_dir, GO_ANNOTATIONS_FILE)

    if(exclude_IEA == "yes"):
        exclude_IEA = ["IEA"]
    else:
        exclude_IEA = []    


    return go_names_file_path, go_annotations_file_path, exclude_IEA

def main(options, cache_db):
    
    go_names_file_path, go_annotations_file_path, exclude_IEA = prepare_folders(options.protein_id_path, options.exclude_IEA, options.output_dir)
    
    proteins = get_proteins(options.protein_id_path)

    proteins_go, all_go, protein_go_dict = get_GO(protein_list=proteins,
                                 exclude=exclude_IEA,
                                 aspect=options.aspect, cache_db=cache_db)
                                 #lack_goes=lack_go_file_path)

    ancestors, all_go = get_ancestors(set([item for sublist in list(proteins_go.values()) for item in sublist]),
                                      ancestors_old={},
                                      all_go=all_go,
                                      aspect=options.aspect,
                                      cache_db=cache_db)

    #create file with names of GO terms
    fill_names(all_go, cache_db, save_file=go_names_file_path)

    #create file with max paths of GO terms
    #get_paths(all_go, path_path=go_max_path_file_path, aspect=options.aspect)

    #create json file with GO protein GO annotations including ancestors
    crate_annotation_file(protein_go_dict, ancestors, go_annotations_file_path)


def get_options():
    parser = OptionParser(description="desc")
    parser.add_option("-i", "--input", dest="protein_id_path", default=None,
                      help="List of proteins for annotations", metavar="FASTA")
    parser.add_option("-e", "--exclude-iea", dest="exclude_IEA", default="no",
                      help="Exclude GO terms with IEA? yes/no", metavar="STRING")
    parser.add_option("-s", "--aspect", dest="aspect", default="F",
                      help="Aspect of GO", metavar="STRING")
    parser.add_option('-o', '--output-dir', dest="output_dir", default='./gbsc_functional_results/',
                      help='Project directory')
    parser.add_option('-c', '--cache-file', dest="cache_file", default='cache.sqlite',
                      help='Cache directory for data downloaded')
    options, args = parser.parse_args()

    return options, args

class DownloadGo:
    class Params:
        def __init__(self):
            self.protein_id_path = None
            self.exclude_IEA = "no"
            self.aspect = "F"
            self.output_dir = "./gbsc_functional_results/"
            self.cache_file = "cache.sqlite"

    def __init__(self):
        self.params = DownloadGo.Params()
        self.cache = Cache()

    def run(self):
        main(self.params, self.cache)


if __name__ == "__main__":
    cache = Cache()
    try:
        options, args = get_options()
        main(options, cache)
    except KeyboardInterrupt:
        print("Shutdown requested...exiting")
        cache.cursor.close()
        cache.connection.close()
    except Exception as e:
        cache.cursor.close()
        cache.connection.close()
        raise e

