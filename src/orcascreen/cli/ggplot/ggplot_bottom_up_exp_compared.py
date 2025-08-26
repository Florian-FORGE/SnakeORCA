import argparse
import textwrap
import os
import threading
import logging
import sys

import orcascreen.matrices as mat

import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from plotnine import (
                      ggplot, aes, geom_boxplot, geom_jitter, scale_fill_brewer, 
                      scale_x_discrete, theme_minimal, theme, xlab, ylab, coord_cartesian, 
                      element_text, facet_wrap, element_blank
                      )
import re



# Helper to get input with a timeout
def prompt_with_timeout(prompt, timeout, result_container):
    try:
        result_container.append(input(prompt))
    except Exception as e:
        logging.error(f"Error during input: {e}")



def natural_key(s):
    # Split the string into a list of strings and integers
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def extract_scores_data(resol: str, create: bool, data_file: str, rename:bool = True, wdir:str = None, expe_name_dict: dict = None):
    """
    Extracts insulation scores and PC1 scores from multiple repositories and saves them to a specified data file.

    Parameters:
    ----------
    resol : str
        The resolution to filter the data.
    create : bool
        If True, creates a new data file; if False, appends to the existing file.
    data_file : str
        The path to the file where the extracted data will be saved.
    rename : bool
        If True, prompts the user to enter experiment names corresponding to each repository.
    
    Returns:
    -------
    None

    Side Effects:
    ----------
    - Creates or appends to a data file with extracted scores.
    - Prompts the user for experiment names if `rename` is True.
    """
    l_repo = sorted([repo 
                     for repo in os.listdir(wdir) 
                     if os.path.isdir( f"{wdir}/{repo}") 
                     and os.path.exists(f"{wdir}/{repo}/pred_log.csv")], key=natural_key)
    
    unique_parts = set()

    for item in l_repo:
        parts = item.split('_') 
        unique_parts.update(parts)

    l_pseudo = sorted(list(unique_parts))
    
    if expe_name_dict is None :
        if rename and len(l_pseudo) < 5 :
            expe_name_dict = {repo : input(f"Enter the experiment name corresponding to {repo}\t") for repo in l_pseudo}
        elif rename :
            l_expe_name = input(f"Enter the experiment name for each repository (comma-separated) - repository order being {l_pseudo}\t").split(",")
            if len(l_expe_name) != len(l_pseudo) :
                raise ValueError(f"There should be as many experiment names as repositories for which data is extracted ({len(l_pseudo)})...Exiting.")
            expe_name_dict = {l_pseudo[i]: l_expe_name[i] for i in range(len(l_pseudo))}
        else :
            expe_name_dict = {pseudo: pseudo for pseudo in l_pseudo}

    print(expe_name_dict)

    mat_comparisons = mat.build_CompareMatrices(filepathref=f"{wdir}/predictions/reference/pred_log.csv", 
                                                filepathcomp=f"{wdir}/matrices_builder/runs.csv")

    df_IS = mat_comparisons.extract_data(data_type="score", score_type="insulation_count", standard_dev=True, 
                                         ref_name="Reference", wanted_pattern="Reference")

    
    df_IS["run_name"] = df_IS["name"]
    
    df_IS["nb_annot"] = df_IS["name"].apply(lambda x: len(str(x).split("_")))
    
    df_IS['name'] = df_IS['name'].apply(lambda x: '_'.join([expe_name_dict.get(part, part) 
                                                            for part in str(x).split('_')]))


    df_PC1 = mat_comparisons.extract_data(data_type="score", score_type="PC1", standard_dev=True, 
                                          ref_name="Reference", wanted_pattern="Reference")
    
    df_PC1["run_name"] = df_PC1["name"]
    
    df_PC1["nb_annot"] = df_PC1["name"].apply(lambda x: len(str(x).split("_")))
    
    df_PC1["name"] = df_PC1["name"].apply(lambda x: '_'.join([expe_name_dict.get(part, part) 
                                                              for part in str(x).split('_')]))

    # Concatenate both DataFrames
    df_combined = pd.concat([df_IS, df_PC1], ignore_index=True)
    df_combined = df_combined[df_combined["resolution"] == resol]

    df_combined.loc[df_combined['name'] == 'Reference', 'nb_annot'] = 0

    df_combined = df_combined.sort_values(by=["run_name", "nb_annot"], 
                                          ascending=[True, True]).reset_index(drop=True)

    # Save to a tab-separated CSV file
    if create :
        if os.path.exists(data_file) :
            df_combined.to_csv(data_file, sep="\t", index=False, mode='w', header=["name", "resolution", "data_type", "values", "reference", "score_type", "mutation_distance", "run_name", "nb_annot"])
        else :
            df_combined.to_csv(data_file, sep="\t", index=False, mode='x', header=["name", "resolution", "data_type", "values", "reference", "score_type", "mutation_distance", "run_name", "nb_annot"])
        create = False
    else :
        df_combined.to_csv(data_file, sep="\t", index=False, mode='a', header=False)


def wrap_text(text: str, width: int, sep: str = " ", new_sep: str = ""):
    # Protecting parts needed for proper printing
    protected_segments = re.findall(r"\$[^$]*\$+", text) 
    for segment in protected_segments:
        text = text.replace(segment, segment.replace(sep, "\uFFFF"))  # Use a rare separator
    
    words = text.split(sep=sep)
    lines = []
    current_line = ""
    for word in words:
        word = word.replace("\uFFFF", sep)  # Restore original spaces inside $...$
        if len(current_line) + len(word) + (1 if current_line else 0) <= width:
            current_line += (new_sep if current_line else "") + word
        else:
            # Check if adding the word would split it more than halfway
            if len(word) > width:
                # There is an exception for a specific syntax
                if word.startswith("$") and word.endswith("$"):
                    if len(current_line) + len(word) + (1 if current_line else 0) <= width:
                        current_line += (new_sep if current_line else "") + word
                    else:
                        if current_line:
                            lines.append(current_line)
                        current_line = word
                    continue
                
                else :
                # For very long words, split at width
                    while len(word) > width:
                        lines.append(word[:width])
                        word = word[width:]
                    current_line = word
            else:
                lines.append(current_line)
                current_line = word
    if current_line:
        lines.append(current_line)
    return "\n".join(lines)


def boxplot_slide(data_file: str, analysis_path:str, output_file: str, score_types: list, 
                  create_data: bool, rename: bool , resol: str = None, wdir:str = None, 
                  expe_name_dict: dict = None, names: list = None):
    """
    Generates a boxplot slide for multiple mutation analyses.

    Parameters:
    ----------
    data_file : str
        The path to the file in which the data is stored.
    analysis_path : str
        The path to the directory in which the analysis plots will be saved.
    output_file : str
        The name of the output PDF file.
    score_types : list
        The list of scores to study (e.g., ['insulation_count', 'PC1']).
    create_data : bool, optional
        If True, the script will create a new data file; if False, it will append to the existing file.
    rename : bool, optional
        If True, prompts the user to enter experiment names corresponding to each repository.
    resol : str, optional
        The resolution to filter the data. If None, defaults to "32Mb".
    wdir : str, optional
        The working directory in which the different experiments are saved. If None, defaults to the current directory.
    expe_name_dict : dict, optional
        A dictionary mapping experiment names to their corresponding display names.
    names : list, optional
        The names (from the 'run_name' column) of the studied experiments to study. If None, all experiments will be included.
        
    Returns:
    -------
    None

    Side Effects:
    ----------
    - Creates a PDF file with boxplots for each score type.
    - Creates a directory for the analysis if it does not exist.
    - Extracts scores data from multiple repositories and saves it to a specified data file.
    """
    wdir = f"{os.path.abspath(os.curdir)}/{wdir}" if wdir is not None else "."
    if create_data :
        if os.path.exists(data_file):
            logging.warning("The data file already exists...Waiting for confirmation.")

            result = []
            thread = threading.Thread(target=prompt_with_timeout, args=("Do you wish to overwrite the already existing file? (yes/no): ", 300, result))
            thread.start()
            thread.join(timeout=300)  # 5 minutes timeout

            if thread.is_alive():
                logging.error("User took too long to respond.")
                raise TimeoutError("Validation timed out after 5 minutes.")

            overwrite = result[0].strip().lower()

            if overwrite not in ['yes', 'y']:
                raise FileExistsError("The data file already exists and overwrite was not confirmed.")
        
        resol = resol if resol else "32Mb"
        extract_scores_data(resol=resol, create=create_data, data_file=data_file, rename=rename, 
                            wdir=wdir, expe_name_dict=expe_name_dict)
    else :
        if not os.path.exists(data_file):
            print(f"Data file {data_file} does not exist. Please create it first.")
            return

    data = pd.read_csv(data_file, sep="\t")
    # data["values"] = data["values"] ** 2 # Uncomment if you want the square deviation

    data= data[data["run_name"].isin(names)] if names is not None else data
           
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)
    
    plots_path = f"{analysis_path}/{output_file}"

    if os.path.exists(plots_path) :
        os.remove(plots_path)
       
    with PdfPages(plots_path, keep_empty=False) as pdf:
        l_steps = sorted(data["nb_annot"].unique(), reverse=True)
        nb_steps = len(l_steps)
        name_width = {step: max(3, 15-data[data["nb_annot"] == step]["run_name"].nunique()//1.5)
                      for step in l_steps}

        for score in score_types:
            _data = data[data["score_type"] == score]
            if _data.empty:
                continue

            _data["wrapped_name"] = _data.apply(
                                                lambda row: wrap_text(row["run_name"], 
                                                                      width=int(name_width[row["nb_annot"]]), 
                                                                      sep="_", 
                                                                      new_sep="-")
                                                if row["run_name"] != "Reference" else row["run_name"],
                                                axis=1
                                                )
            _data["nb_annot"] = pd.Categorical(_data["nb_annot"])

            ymin = _data["values"].min()
            ymin = 1.05 * ymin if ymin < 0 else .95 * ymin
            ymax = _data["values"].max()
            ymax = 1.05 * ymax if ymax > 0 else .95 * ymax

            score_label = "IS" if score == "insulation_count" else score
            add_fig_size = nb_steps*5
            add_font_size = nb_steps*8

            plot = (
                ggplot(_data, aes(x='wrapped_name', y='values', fill='nb_annot')) +
                geom_boxplot() +
                geom_jitter(color="black", size=0.4, alpha=0.9) +
                scale_fill_brewer(type='div', palette='Spectral') +
                scale_x_discrete() +
                theme_minimal() +
                theme(
                    panel_spacing_y=.5*(nb_steps - 1),                              # Adjusts the spacing between facets
                    figure_size=(40+add_fig_size, 55+2*add_fig_size),               # Adjust as needed
                    legend_position='none',
                    axis_text_x=element_text(size=24+add_font_size, weight="bold"),
                    axis_text_y=element_text(size=24+add_font_size, weight="bold"),
                    axis_title_y=element_text(size=40+add_font_size, weight="bold"),
                    axis_title_x=element_text(size=40+add_font_size, weight="bold"),
                ) +
                xlab("") +
                ylab(f"{score_label} deviation") +
                coord_cartesian(ylim=(ymin, ymax)) +
                facet_wrap('~nb_annot', ncol=1, scales='free_x')  # ncol can be adjusted
            )
            
            plot += theme(
                          strip_text = element_blank(),        # Removes the text labels
                          strip_background = element_blank(),  # Removes the background box
                          plot_margin=0.005
                          )
            
            correspondence = "\t".join(
                                       f"{row['run_name']} → {row['name']}"
                                       for _, row in data.drop_duplicates(subset=["name", "run_name"]).iterrows()
                                       if (row["run_name"] != "Reference" and row["nb_annot"] == 1)
                                       )
            correspondence = wrap_text(correspondence, width=66+nb_steps*3, sep="\t", new_sep=" • ")
            
            fig = plot.draw()
            fig.subplots_adjust(top=0.999, bottom=0.01, left=0.07)
            fig.text(.5, .99, correspondence, ha='center', va='top', 
                     fontsize=24+add_font_size, fontweight='bold', style='italic', 
                     bbox=dict(boxstyle='round,pad=0.6', facecolor='white', alpha=0.8))
            pdf.savefig(fig)
                
        pdf.close()

        log_path = f"{wdir}/comparison.log"
        if os.path.exists(log_path) :
            os.remove(log_path)
        
        with open(log_path, "w") as fglog :
            fglog.write(f"# The comparison graphs are consultable in:\n{plots_path}")




def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    
    parser.add_argument("--data_file",
                       required=True, help="The path to the file in which the data is stored")
    parser.add_argument("--analysis_path",
                        required=True, help='The path to the directory in which the analysis plots will be saved.')
    parser.add_argument("--output_file",
                        required=True, help='The name of the output PDF file.')
    parser.add_argument("--score_types",
                        required=True, help="The list of scores to study (format : 'score1,score2,...')")
    parser.add_argument("--create_data",
                        required=False,
                        action='store_true', help="If True, the script will create a new data file; if False, it will append to the existing file.")
    parser.add_argument("--rename",
                        required=False,
                        action='store_true', help="If True, prompts the user to enter experiment names corresponding to each repository.")
    parser.add_argument("--resol",
                        required=False, help="The resolution to filter the data. If None, defaults to '32Mb'.")
    parser.add_argument("--wdir",
                        required=False, help="The working directory in which the different experiments are saved, in case it not the current directory.")
    parser.add_argument("--expe_names",
                        required=False, help="The dict of experiments names (format : 'expe1:name1,expe2:name2,...'). Note : if there are LateX elements '$...$' use '\$...\$'.")
    parser.add_argument("--studied_names",
                        required=False, help="The names (from the 'run_name' column) of the studied experiments to study (format : 'expe1,expe2,...'). If None, all experiments will be included. Note : if there are LateX elements '$...$' use '\$...\$'.")
    
    
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()

    score_types = args.score_types.split(",")
    
    expe_names = None
    if args.expe_names is not None :
        expe_names = [elt.split(":") for elt in args.expe_names.split(",")]
        expe_names = {elt[0]: elt[1] for elt in expe_names}
    
    studied_names = args.studied_names.split(",") if args.studied_names is not None else None
    
    boxplot_slide(data_file=args.data_file, 
                  analysis_path=args.analysis_path, 
                  output_file=args.output_file, 
                  score_types=score_types, 
                  create_data=args.create_data, 
                  rename=args.rename, 
                  resol=args.resol,
                  wdir=args.wdir,
                  expe_name_dict=expe_names,
                  names=studied_names)


if __name__ == '__main__':
    main()
