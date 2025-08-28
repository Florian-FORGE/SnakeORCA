import argparse
import textwrap
import random
import os
from Bio import SeqIO

import logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format="%(asctime)s - %(levelname)s - %(message)s"
)

import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
This script reads a trace file (in tsv format) describing mutations and generates a bed file with relative positions of mutations.
It handles overlapping intervals by merging them, in case there are many files.
It saves the resulting bed file and a log file summarizing the operations performed.
It is designed to be used in a pipeline for genome mutation analysis.

Usage:
    python annot_from_traces.py --trace trace1.tsv,trace2.tsv 
                                --mut_path /path/to/output 
                                --chrom_name fake_chr

Dependencies:
    pandas
    biopython
"""


def check_overlapping(df: pd.DataFrame):
    row_modif = False
    
    df = df.sort_values(by='start').reset_index(drop=True)
    # df['length'] = df['end'] - df['start']
    # df = df.sort_values(by=['start', 'length'], ascending=[True, False]).reset_index(drop=True)

    # Initialize empty list to act as a stack
    stack = []

    # Iterate through rows
    for idx, row in df.iterrows():
        if not stack:
            stack.append(row)
        else:
            last = stack[-1]
            diff = last['end'] - row['start']

            # Only add if current start is after last end
            if diff < 0:
                stack.append(row)
            
            # If it is not the case the interval should be enlarged
            else :
                if last['end'] < row['end'] :
                    stack[-1]['end'] = row['end']
                    stack[-1]['sequence'] += row['sequence'][diff:]


    # Final filtered DataFrame
    df_filtered = pd.DataFrame(stack)
    df_filtered = df_filtered.sort_values(by='start').reset_index(drop=True)

    return df_filtered, row_modif


def read_annotations_from_trace(trace: str, chrom_name: str):
    """
    Read a database describing mutations and return the needed information in a database.
    """
    df = pd.read_csv(trace, sep="\t")
    df = df[['chrom', 'start', 'end', 'name', 'strand', 'operation', 'ref_seq']]  # Keep only these columns

    df = df.rename(columns={'ref_seq': 'sequence'})
    df["operation"] = "insertion"
    df.loc[df["chrom"] == chrom_name, "chrom"] = "fake_chr"
    df['start'] = pd.to_numeric(df['start']) - 1                                  # Convert to 0-based for processing (end is exclusive in bed format)
    
    return df


def bed_from_traces(traces: list, mut_path: str, chrom_name: str):
    """
    Create a bed file with relative positions of mutations from a list of traces.
    Handle overlapping intervals by merging them.
    Save the resulting bed file and a log file summarizing the operations performed.

    Parameters
    ----------
        traces (list) : List of paths to trace files (in tsv format).
        mut_path (str) : Path to the output directory to store the relative sequence and bed.
        chrom_name (str) : Name of the chromosome in the trace file (most likely 'fake_chr$').
    
    Returns
    -------
        None
    
    Side effects
    ------------
        Creates output directory if it does not exist.
        Saves the resulting bed file and a log file in the specified output directory.
    
    Notes
    -----
        If multiple traces are provided, overlapping intervals are merged.
        The output bed file is named 'relative_position_mutations.csv'.
        A log file named 'traces_resume.log' is created to summarize the operations performed.
    """
    if not os.path.exists(mut_path):
        os.makedirs(mut_path)

    new_bed_path = f"{mut_path}/relative_position_mutations.csv"
    if os.path.exists(new_bed_path) :
        os.remove(new_bed_path)
    
    raw_bed_path = f"{mut_path}/raw_intervals.csv"
    if os.path.exists(raw_bed_path) :
        os.remove(raw_bed_path)
    
    # Get the annotations
    create = True
    df_combined = pd.DataFrame()
    for trace in traces :
        df = read_annotations_from_trace(trace=trace, chrom_name=chrom_name)
        
        if create :
            if os.path.exists(raw_bed_path) :
                df.to_csv(raw_bed_path, sep="\t", index=False, mode='w', header=True)
            else :
                df.to_csv(raw_bed_path, sep="\t", index=False, mode='x', header=True)
            create = False
        else :
            df.to_csv(raw_bed_path, sep="\t", index=False, mode='a', header=False)
                
        df_combined = pd.concat([df_combined, df], ignore_index=True)

    df_final, row_modif = check_overlapping(df=df_combined)

    if row_modif:
        logging.info("Since there was overlapping intervals, merging was applied...Proceeding")

    if os.path.exists(new_bed_path) :
        df_final.to_csv(new_bed_path, sep="\t", index=False, mode='w', header=True)
    else :
        df_final.to_csv(new_bed_path, sep="\t", index=False, mode='x', header=True)

    mutation_line = f"Relative_tsv\t{new_bed_path}"

    # Save data about what has been done
    log_path = f"{mut_path}/traces_resume.log"
    trace_data = "\n".join([f"- {trace}" for trace in traces])
    with open(log_path, "w") as flog :
        flog.write("The trace(s) used for the annotations was/were stored in : "
                   f"{trace_data}"
                   "The relevant mutations (with relative positions) to apply "
                   f"were stored in : {new_bed_path}\n")
    
    global_log_path = "/".join(mut_path.split("/")[:-1]) if mut_path.split("/")[-1] == "genome" else mut_path
    global_log_path += "/annot_f_traces.log"
    if os.path.exists(global_log_path) :
        os.remove(global_log_path)
    
    with open(global_log_path, "w") as fglog :
        fglog.write(f"{mutation_line}\n"
                    f"Chrom_name\tfake_chr")

 

def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    
    parser.add_argument("--trace",
                        required=True, help="the trace in tsv format.")
    parser.add_argument("--mut_path",
                        required=True, help="the path to the output directory to store the relative sequence and bed.")
    parser.add_argument("--chrom_name",
                        required=True, help="the name (in the trace file -- most likely 'fake_chr$' -- of the chomosome which is studied.")
    
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()
    
    traces = args.trace.split(",")
    
    bed_from_traces(traces=traces,  
                    mut_path=args.mut_path,
                    chrom_name=args.chrom_name)


if __name__ == '__main__':
    main()
