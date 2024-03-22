import json
import sys


def compare(bruteforce_output_file: str, suffixarray_output_file: str):
    with open(bruteforce_output_file, 'r') as fp_bruteforce:
        bruteforce_data = json.load(fp_bruteforce)
        with open(suffixarray_output_file, 'r') as fp_suffixarray:
            suffixarray_data = json.load(fp_suffixarray)

            for entry in bruteforce_data["result"]:
                sequence = entry["sequence"]
                accessions = entry["uniprot_accessions"]
                if accessions:  # if no accessions found, we don't return it in the rust code, so skip it
                    match_found = False
                    for sa_entry in suffixarray_data["result"]:
                        if sequence == sa_entry["sequence"] and ("visited" not in sa_entry.keys() or not sa_entry["visited"]):
                            match_found = True
                            if sorted(accessions) == sorted(sa_entry["uniprot_accessions"]):
                                sa_entry["visited"] = True  # add a boolean to indicate if we visited that entry
                                break
                    else:
                        if match_found:
                            print(f"The bruteforce and suffixarray results differ for peptide {sequence}",
                                  file=sys.stderr)
                        else:
                            print(f"Sequence {sequence} not found in suffixarray data", file=sys.stderr)
                        assert False

            # check if we actually visited all the nodes from the suffix array results
            for entry in suffixarray_data["result"]:
                assert "visited" in entry.keys() and sa_entry["visited"]


if __name__ == "__main__":
    compare("swissprot_no_mch_inexact_bruteforce.json", "swissprot_no_mch_inexact_suffixarray.json")
