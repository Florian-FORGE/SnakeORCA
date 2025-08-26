import argparse
import textwrap
import os
import sys

import orcascreen.process_sequence as ps
import orcascreen.matrices as mat

import pandas as pd




def predict_and_orcarun_descript(resol_model: str,
                                 mutate_log_path: str, 
                                 mpos: int, 
                                 cool_resol: int, 
                                 strict: bool, 
                                 padding_chr: str, 
                                 use_cuda: bool, 
                                 use_memmapgenome: bool, 
                                 pred_path: str, 
                                 abs_rdm_log_path: str, 
                                 builder_path: str) :
    """
    Generates predictions and associated OrcaRun descriptions from a given sequence.
    This function processes a sequence to generate predictions for a set of mutations 
    (wild-type and random mutations) and creates associated OrcaRun descriptions. It 
    also generates a reference OrcaRun description for the wild-type sequence.

    Args:
        resol_model (str): Resolution model to be used for predictions.
        mutate_log_path (str): Path to a file storing paths to the logs of the mutations (and reference).
        mpos (int): Mutation position.
        cool_resol (int): Resolution for the cool file.
        strict (bool): Whether to enforce strict processing.
        padding_chr (str): Padding character for chromosome sequences.
        use_cuda (bool): Whether to use CUDA for computations.
        use_memmapgenome (bool): Whether to use memory-mapped genome files.
        pred_path (str): Path to store prediction outputs.
        abs_rdm_log_path (str): Path to a file in which data abot the reference FASTA file and studied chromosome are stored.
        builder_path (str): Path to store the generated OrcaRun description files.
    Raises:
        Exception: If there is an error creating the builder directory.
    Outputs:
        - Generates prediction files for each mutation and the reference sequence.
        - Creates OrcaRun description files (`orcarun.csv` and `ref_orcarun.csv`) 
          containing metadata for the predictions.
    """
    data = []
    
    head = ["name", "list_resol", "path", "gtype", "genome", "obj_type", "refgenome", "mtype", "trace_path"]

    with open(abs_rdm_log_path, "r") as ffa:
        lines = [line.strip().split("\t") for line in ffa if not line.startswith("#")]
        data_rel = {lines[i][0].lower() : lines[i][1] for i in range(len(lines))}
        ref_fasta = data_rel["relative_fasta"] if "relative_fasta" in data_rel.keys() else None
        chrom = data_rel["chrom_name"] if "chrom_name" in data_rel.keys() else None

    
    ref_path = f"{pred_path}/reference"
    if os.path.exists(ref_path) and os.path.isdir(ref_path):
        os.rmdir(ref_path)


    ps.main(chrom=chrom,
            output_prefix="reference",
            mutation="wt",
            resol_model=resol_model,
            mpos=mpos,
            fasta=ref_fasta,
            cool_resol=cool_resol,
            strict=strict,
            padding_chr=padding_chr,
            use_cuda=use_cuda,
            use_memmapgenome=use_memmapgenome,
            pred_path=pred_path)

    l_resol = mat.extract_resol_asc(ref_path)

    data_ref = [f"Reference", 
                f"{l_resol}", 
                f"{ref_path}", 
                f"wt", 
                f"{ref_fasta}", 
                f"MatrixView", 
                f"{ref_fasta}", 
                f"OrcaMatrix",
                f"{None}"]
    
    data.append(data_ref)
    
    pred_log = f"{ref_path}/pred_log.csv"

    df_ref = pd.DataFrame([data_ref],columns=head)
    df_ref.to_csv(pred_log, sep="\t", index=False, mode='x', header=True)

    
    with open(mutate_log_path, "r") as fmut:
        lines = [line.strip().split("\t") for line in fmut if not line.startswith("#")]

    for fasta_path, trace_path in lines:
        name = fasta_path.split("/")[-2]
        expe_path = f"{pred_path}/{name}"
        if os.path.exists(expe_path) and os.path.isdir(expe_path):
            os.rmdir(expe_path)

        ps.main(chrom=chrom,
                output_prefix=name,
                mutation=name,
                resol_model=resol_model,
                mpos=mpos,
                fasta=fasta_path,
                cool_resol=cool_resol,
                strict=strict,
                padding_chr=padding_chr,
                use_cuda=use_cuda,
                use_memmapgenome=use_memmapgenome,
                pred_path=pred_path)
        
        
        l_resol = mat.extract_resol_asc(expe_path)

        expe_data = [f"{name}", 
                     f"{l_resol}", 
                     f"{expe_path}", 
                     f"{name}", 
                     f"{fasta_path}", 
                     f"MatrixView", 
                     f"{ref_fasta}", 
                     f"OrcaMatrix", 
                     f"{trace_path}"]
                     
        data.append(expe_data)
        
        pred_log = f"{expe_path}/pred_log.csv"

        df_expe = pd.DataFrame([expe_data],columns=head)
        df_expe.to_csv(pred_log, sep="\t", index=False, mode='x', header=True)
    
    
    
    if builder_path:
        try:
            builder_path = os.path.abspath(builder_path)
            os.makedirs(builder_path, exist_ok=True)
            print(f"Directory ensured: {builder_path}")
            
            if os.path.exists(builder_path):
                print(f"Directory exists: {builder_path}")
            else:
                print(f"Directory does not exist: {builder_path}")
        except Exception as e:
            print(f"Error creating directory: {e}")
    else:
        print("builder_path is not defined or is empty.")

    df = pd.DataFrame(data,columns=head)

    df.to_csv(f"{builder_path}/runs.csv", sep="\t", index=False, mode='x', header=True)

    
    global_log_path = "/".join(mutate_log_path.split("/")[:-1]) if mutate_log_path.split("/")[-1] == "mutate.log" else mutate_log_path
    global_log_path += "/prediction.log"
    if os.path.exists(global_log_path) :
        os.remove(global_log_path)
    
    with open(global_log_path, "w") as fglog :
        fglog.write(f"# The file needed to build the matrices:\n{builder_path}")



def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    
    # parser.add_argument("--pred_prefix",
    #                     required=True, help="The the prediction prefix that is used to differentiate these runs from others.")
    parser.add_argument("--resol_model", 
                        required=False, help="The resolution model to use (either '32Mb' or '256Mb').")
    parser.add_argument("--mutate_log_path", 
                        required=True, help="The path to the directory in which all the genome directories are stored.")
    parser.add_argument("--mpos", 
                        required=True, type=int, help="The coordinate to zoom into for multiscale prediction.")
    parser.add_argument("--pred_path",
                        required=True, help="Path to the directory where the predictions will be saved. Defaults to None. If None, the predictions will be saved in the same directory as the genomes.")
    parser.add_argument("--abs_rdm_log_path", 
                        required=True, 
                        help="The full path to the .log file of the relative genome.")
    parser.add_argument("--builder_path", 
                        required=True, help="Path to the directory where the data to build the matrices will be saved. Defaults to None. If None, the data will be saved in the same directory as the genomes.")
    parser.add_argument("--cool_resol", 
                        required=False, type=int, default=128_000, 
                        help="The resolution of a .mcool file that could be used for comparison. Used to achieve good alignment of bins (the start poition should be divisible by cool_resol). Defaults to 128_000.")
    parser.add_argument("--strict", 
                        required=False, type=bool, default=False, 
                        help="Whether the start position should be used directly or not. Defaults to False.")
    parser.add_argument("--padding_chr", 
                        required=False, help="If resol_model is '256Mb', padding is generally needed to fill the sequence to 256Mb. The padding sequence will be extracted from the padding_chr.")
    parser.add_argument("--no_cuda", 
                        required=False,
                        help="Whether to use CUDA for GPU acceleration. If False, CUDA is used. Defaults to False.")
    parser.add_argument("--use_memmapgenome", 
                        required=False, type=bool, default=True, 
                        help=" Whether to use memmory-mapped genome. Defaults to True.")
    
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()
    
    no_cuda = bool(args.no_cuda.lower() == "true") if args.no_cuda is not None else True
    use_cuda = not no_cuda
    
    predict_and_orcarun_descript(resol_model=args.resol_model, 
                                 mutate_log_path=args.mutate_log_path, 
                                 mpos=args.mpos, 
                                 cool_resol=args.cool_resol, 
                                 strict=args.strict, 
                                 padding_chr=args.padding_chr, 
                                 use_cuda=use_cuda, 
                                 use_memmapgenome=args.use_memmapgenome, 
                                 pred_path=args.pred_path, 
                                 abs_rdm_log_path=args.abs_rdm_log_path, 
                                 builder_path=args.builder_path)



if __name__ == '__main__':
    main()
