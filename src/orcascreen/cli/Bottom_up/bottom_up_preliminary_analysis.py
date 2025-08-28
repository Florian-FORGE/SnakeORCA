import argparse
import textwrap
import os
import sys

import orcascreen.matrices as mat
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd

import random as rd


"""
This script produces a pdf file containing the analysis plots for a given run.
It can be used as a command line tool or as a function in a python script.

Usage :
    python bottom_up_preliminary_analysis.py --expe_descrip "Description of the experiment" \
                                             --run_path "path/to/run_file.tsv" 
                                             --analysis_path "path/to/analysis_directory" 
                                             --l_score_types "insulation_count,PC1" 
                                             --l_resol "1Mb,2Mb,4Mb"
                                             --show_compartments True

Dependencies :
    -matplotlib
    -pandas
    -orcascreen (custom library)

Notes :
    -The run file must be a tsv file with a header.
    -The analysis directory will be created if it does not exist.
"""


def analysis_slide(expe_descrip: str, run_path: str, analysis_path:str, 
                   l_score_types: list, l_resol: list = None, 
                   show_compartments: bool = False):
    """
    """
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)
    
    df = pd.read_csv(run_path, header=0, sep='\t')

    row = next(df.itertuples(index=False))

    mat_view = mat._build_MatrixView_(row=row)

    mutation = True
    if mat_view.list_mutations is None :
        mutation = False

    if l_resol is None or l_resol == ['None'] :
        l_resol = []
        if "insulation_count" in l_score_types or "insulation_correl" in l_score_types: 
            l_resol += ["1Mb", "2Mb", "4Mb"]
        if "PC1" in l_score_types : 
            l_resol += ["8Mb", "16Mb", "32Mb"]

    if mutation :
        nb_mut_max = mat_view.di[l_resol[-1]].nb_mutated_pb
        if nb_mut_max > 1_000_000 :
            nb_mut_max = f"{nb_mut_max/1_000_000:.1f}Mb"
        elif nb_mut_max > 1_000 :
            nb_mut_max = f"{nb_mut_max/1_000:.1f}kb"
        else :
            nb_mut_max = f"{nb_mut_max}pb"
    
    
    log_info = ""
    for resol in l_resol :
        plots_path = f"{analysis_path}/plots_{resol}.pdf"
        saddle_path = f"{analysis_path}/saddle_{resol}.pdf" \
                                    if resol in ["8Mb", "16Mb", "32Mb"] else None

        if mutation :
            nb_mut_resol = mat_view.di[resol].nb_mutated_pb
            if nb_mut_resol > 1_000_000 :
                nb_mut_resol = f"{nb_mut_resol/1_000_000:.1f}Mb"
            elif nb_mut_resol > 1_000 :
                nb_mut_resol = f"{nb_mut_resol/1_000:.1f}kb"
            else :
                nb_mut_resol = f"{nb_mut_resol}pb"
        
        
        if os.path.exists(plots_path) :
            os.remove(plots_path)
        
        with PdfPages(plots_path, keep_empty=False) as pdf:
        
            nb_scores = len(l_score_types)
            height_ratios = [0.9] + [0.75/nb_scores for _ in range(nb_scores)]
            width_ratios = [45, 5, 45, 5, 45]
            gs = GridSpec(nrows= 1 + nb_scores, ncols=5, height_ratios=height_ratios, width_ratios=width_ratios)
            f = plt.figure(clear=True, figsize=(60, 33.75))
            f.subplots_adjust(hspace=.5, top=0.92, bottom=0.05)
            if mutation : 
                f.suptitle(f"{expe_descrip} ({nb_mut_max} mutated - {nb_mut_resol} locally)--{resol}.", fontsize=48)
            else :
                f.suptitle(f"{expe_descrip} --{resol}.", fontsize=48)

            
            log_info += f"To produce the plots in plots_{resol}.pdf the following method and arguments were used :\n"

            ax = mat_view.di[resol].heatmap_plot(gs=gs, f=f, j=2, compartment=show_compartments, mutation=mutation)
            log_info += f"mat_view.di[{resol}].heatmap_plot(gs=gs, f=f, compartment=True, mutations={mutation})\n"
            ax.tick_params(axis='both', labelsize=32)

            mat_view.di[resol].hist_mutations(gs=gs, f=f, j=3, title=f"Mutation repartition\n")
            log_info += f"mat_view.di[{resol}].hist_mutations(gs=gs, f=f, j=1, title='Mutation repartition')\n"
            
            k=1
            for score_type in l_score_types:
                title = "Insultion Score (IS)" if score_type == "insultation_count" else score_type
                ax = f.add_subplot(gs[k, :])
                mat_view.di[resol]._score_plot(gs=gs, f=f, ax=ax, title=title, score_type=score_type)
                ax.tick_params(axis='both', labelsize=32)
                log_info += f"mat_view.di[{resol}]._score_plot(gs=gs, f=f, ax=ax, title={title}, score_type={score_type})\n"
                k+=1
            
        pdf.savefig()
        pdf.close()
        log_info += f"\n"
    
        if saddle_path is not None :
            if os.path.exists(saddle_path) :
                os.remove(saddle_path)
            
            with PdfPages(saddle_path, keep_empty=False) as pdf:
                mat_view.di[resol].saddle_plot(mutation=mutation)
                pdf.savefig()
                pdf.close()

    log_path = f"{analysis_path}/plots.log"
    with open(log_path, "w") as fout:
        fout.write(log_info)

    
    log_path = f"{analysis_path}/analysis.log"
    if os.path.exists(log_path) :
        os.remove(log_path)
    
    with open(log_path, "w") as fglog :
        fglog.write(f"# The following function has successfully been executed:\n")
        fglog.write(f"analysis_slide(expe_descrip={expe_descrip}, analysis_path={analysis_path}, l_score_types={l_score_types}, "
                    f"l_resol={l_resol}, show_compartments={show_compartments})\n")





def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    
    parser.add_argument("--expe_descrip",
                        required=True, help="The name of the experiment. It will be used in the title inside of the output pdf.")
    parser.add_argument("--run_path",
                       required=True, help="The path to the file in which the run's data is stored.")
    parser.add_argument("--analysis_path",
                        required=True, help='The path to the directory in which the analysis plots will be saved.')
    parser.add_argument("--l_score_types",
                        required=True, help='The list of scores for which plots should specifically be done.')
    parser.add_argument("--l_resol",
                        required=False, help="The list of resolutions to study. If not specified, the more " \
                                                        "representative resolutions for the given score types will " \
                                                        "automatically be selected.")
    parser.add_argument("--show_compartments",
                        required=False, help="Whether to show the compartments in the plots. If not specified, it will be set to False.")
    
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()

    l_score_types = args.l_score_types.split(",")
    l_resol = args.l_resol.split(",") if args.l_resol is not None else None
    show_compartments = bool(args.show_compartments.lower() == "true") if args.show_compartments is not None else False
   
    analysis_slide(expe_descrip=args.expe_descrip, 
                   run_path=args.run_path, 
                   analysis_path=args.analysis_path, 
                   l_score_types=l_score_types,
                   l_resol=l_resol,
                   show_compartments=show_compartments)


if __name__ == '__main__':
    main()

   

