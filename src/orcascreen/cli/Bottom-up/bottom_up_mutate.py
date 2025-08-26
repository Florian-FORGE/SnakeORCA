import argparse
import textwrap
import os
import sys

import orcascreen.mutation as m
import orcascreen.mutate_with_rdm as mm

from pysam import FastaFile
from Bio import SeqIO


"""
This script applies mutations to a reference sequence and generates random mutations if specified.
It processes a given region of a reference sequence by applying mutations from a BED file or a TSV file.
It also generates a specified number of random mutations and applies them to the sequence.
The mutated sequences and their corresponding traces are saved to the specified output path.

Usage:
    python bottom_up_mutate.py --abs_rdm_log_path /path/to/rdm.log 
                               --annot_f_traces_log_path /path/to/annot_f_traces.log
                               --mut_path /path/to/output

Dependencies:
    pandas
    biopython
    pysam
    orcascreen (custom library)

Notes:
    - The function ensures that random mutations do not overlap with existing mutations.
    - The random seed for generating mutations is incremented for each random mutation set.
"""



def mutate(abs_rdm_log_path: str, annot_f_traces_log_path:str, mut_path: str):
    """
    Applies mutations to a reference sequence and generates random mutations if specified.
    This function processes a given region of a reference sequence by applying mutations
    from a BED file or a TSV file. It also generates a specified number of random mutations
    and applies them to the sequence. The mutated sequences and their corresponding traces
    are saved to the specified output path.
    
    Args:
        abs_rdm_log_path (str): File path leading to the file containing data about the 
        fasta and mutations to apply.
        mut_path (str): Directory path where the mutated sequences and traces will be saved.
        muttype (str): Type of mutation to apply (e.g., SNP, insertion, deletion).
        sequence (str): In case the mutation type is 'insertion', the seaquence to insert.
    Raises:
        FileNotFoundError: If the reference sequence file does not exist.
        ValueError: If the region format is invalid or if mutation data is missing.
    Outputs:
        - Mutated sequences are saved in FASTA format under `mut_path/{mutation_name}/sequence.fa`.
        - Mutation traces are saved as CSV files under `mut_path/{mutation_name}/trace_{mutation_name}.csv`.
    Notes:
        - The function ensures that random mutations do not overlap with existing mutations.
        - The random seed for generating mutations is incremented for each random mutation set.
    """
    with open(abs_rdm_log_path, "r") as ffa, open(annot_f_traces_log_path, "r") as ftsv :
        lines = [line.strip().split("\t") for file in (ffa, ftsv) for line in file if not line.startswith("#")]
    
    filepaths = {lines[i][0].lower() : lines[i][1] for i in range(len(lines))}

    relative_tsv = filepaths["relative_tsv"] if "relative_tsv" in filepaths.keys() else None
    relative_fasta = filepaths["relative_fasta"] if "relative_fasta" in filepaths.keys() else None


    mutations = mm.read_mutations_from_tsv(relative_tsv)
    
    fasta_handle = FastaFile(relative_fasta)

    mutator = m.Mutator(fasta_handle, mutations)
    
    mutator.mutate()
    seq_records = mutator.get_mutated_chromosome_records()

    # name = os.path.basename(os.path.dirname(mut_path))
    name = os.path.basename(mut_path)
        
    if not os.path.exists(f"{mut_path}"):
        os.makedirs(f"{mut_path}")
        
    output_path = f"{mut_path}/sequence.fa"
    SeqIO.write(seq_records, output_path, "fasta")
    
    trace = mutator.get_trace()
    trace.to_csv(f"{mut_path}/trace_{name}.csv", 
            sep="\t", index=False, header=True)
    
    global_log_path = "/".join(mut_path.split("/")[:-1]) if mut_path.split("/")[-1] == "genome" else mut_path
    global_log_path += "/mutate.log"
    if os.path.exists(global_log_path) :
        os.remove(global_log_path)
    
    with open(global_log_path, "w") as fglog :
        fglog.write(f"# Relative fasta file and trace:\n")
        fglog.write(f"{mut_path}/sequence.fa\t{mut_path}/trace_{name}.csv\n")
        


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    
    parser.add_argument('--abs_rdm_log_path',
                        required=True, help="The path to a .log file in which the path to the genome as follow: 'name\tpath'.")
    parser.add_argument('--annot_f_traces_log_path',
                        required=True, help="The path to a .log file in which the path to the mutation file to use are given as follow: 'name\tpath'.")
    parser.add_argument("--mut_path",
                        required=True, help="The path to the output directory to store the relative sequences and traces.")
    
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    
    mutate(abs_rdm_log_path=args.abs_rdm_log_path, 
           annot_f_traces_log_path=args.annot_f_traces_log_path, 
           mut_path=args.mut_path)


if __name__ == '__main__':
    main()
