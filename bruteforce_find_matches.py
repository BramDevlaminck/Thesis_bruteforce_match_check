import json
import multiprocessing
from multiprocessing.shared_memory import ShareableList

import regex as re
from tqdm import tqdm


def read_tsv(file: str) -> list[tuple[str, str]]:
    proteins = []
    with open(file) as fp:
        for line in fp:
            (accession_nr, _taxon_id, sequence) = line.rstrip("\n").split("\t")
            proteins.append((accession_nr, sequence))

    return proteins

def search_peptides(process_num: int, proteins: ShareableList, equalize_i_and_l: bool, peptides: list[str], return_dict):
    for peptide in peptides:
        peptide = peptide.rstrip("\n").upper()
        search_peptide = peptide
        if equalize_i_and_l:
            search_peptide = re.sub("[IL]", "[IL]", search_peptide)

        accessions = []
        for accession_nr, sequence in proteins:
            for _ in re.findall(search_peptide, sequence, overlapped=True):
                accessions.append(accession_nr)

        return_dict[process_num].append({
            "sequence": peptide,
            "uniprot_accessions": accessions
        })

def search_file(proteins: list[tuple[str, str]], peptide_file: str, output_file: str, equalize_i_and_l: bool = False):

    results = {"result": []}

    num_cpu = multiprocessing.cpu_count()

    peptides = []
    with open(peptide_file) as fp:
        for peptide in fp:
            peptides.append(peptide)

    block_size = len(peptides) // num_cpu

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for i in range(num_cpu - 1):
        return_dict[i] = manager.list()
        p = multiprocessing.Process(target=search_peptides, args=(i, proteins, equalize_i_and_l, peptides[i*block_size:(i+1)*block_size], return_dict))
        jobs.append(p)
        p.start()

    # handle last iteration separately to ensure that we don't skip some data at the end
    return_dict[num_cpu - 1] = manager.list()
    p = multiprocessing.Process(target=search_peptides, args=(num_cpu - 1, proteins, equalize_i_and_l, peptides[(num_cpu - 1) * block_size:], return_dict))
    jobs.append(p)
    p.start()

    for proc in jobs:
        proc.join()

    # flatten output from subprocesses and save result
    results["result"] = [j for sub in return_dict.values() for j in sub]

    with open(output_file, "w") as output_fp:
        json.dump(results, output_fp)


if __name__ == '__main__':
    data = read_tsv(
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database_minimal.tsv")
    search_file(data,
                "small_swissprot_no_mch.tsv",
                "swissprot_no_mch_inexact_bruteforce.json",
                True)
