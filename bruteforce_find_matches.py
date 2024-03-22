import json
import regex as re
from tqdm import tqdm


def read_tsv(file: str) -> list[tuple[str, str]]:
    proteins = []
    with open(file) as fp:
        for line in fp:
            (accession_nr, _taxon_id, sequence) = line.rstrip("\n").split("\t")
            proteins.append((accession_nr, sequence))

    return proteins


def search_file(proteins: list[tuple[str, str]], peptide_file: str, output_file: str, equalize_i_and_l: bool = False):
    results = {"result": []}
    with open(peptide_file) as fp:
        for peptide in tqdm(fp):
            peptide = peptide.rstrip("\n").upper()
            search_peptide = peptide
            if equalize_i_and_l:
                search_peptide = re.sub("[IL]", "[IL]", search_peptide)

            accessions = []
            for accession_nr, sequence in proteins:
                for _ in re.findall(search_peptide, sequence, overlapped=True):
                    accessions.append(accession_nr)

            results["result"].append({
                "sequence": peptide,
                "uniprot_accessions": accessions
            })

    with open(output_file, "w") as output_fp:
        json.dump(results, output_fp)


if __name__ == '__main__':
    data = read_tsv(
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database_minimal.tsv")
    search_file(data,
                "small_swissprot_no_mch.tsv",
                "swissprot_no_mch_inexact_bruteforce.json",
                True)
