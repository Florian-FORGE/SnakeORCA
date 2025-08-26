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

from seaborn import boxenplot, color_palette
import re


"""
This script generates boxplot slides comparing multiple mutation analyses based on insulation scores and PC1 scores.
It extracts scores from multiple repositories, processes the data, and creates visualizations in a PDF format.
It is designed to be used in a pipeline for genome mutation analysis.

Usage:
    python bottom_up_expe_compared.py --descrip "Description of the analysis"
                                      --data_file /path/to/data_file.tsv
                                      --analysis_path /path/to/analysis
                                      --output_file output.pdf
                                      --score_types insulation_count,PC1
                                      [--create_data]
                                      [--rename]
                                      [--resol 32Mb]
                                      [--wdir /path/to/working_directory]
                                      [--expe_names expe1:name1,expe2:name2,...]

Dependencies:
    pandas
    matplotlib
    seaborn
    orcascreen

Note:
    - If the data file already exists and --create_data is specified, the user will be prompted to confirm overwriting.
    - If --rename is specified, the user will be prompted to enter experiment names corresponding to each repository.
    - The script assumes that the repositories follow a specific naming convention for extracting experiment names.
"""



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


def boxplot_slide(descrip: str, data_file: str, analysis_path:str, output_file: str, 
                  score_types: list, create_data: bool, rename: bool , resol: str = None, 
                  wdir:str = None, expe_name_dict: dict = None):
    """
    Generates a boxplot slide for multiple mutation analyses.

    Parameters:
    ----------
    descrip : str
        The description that will be used in the title inside of the output pdf.
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
    
    # Retrieving the names of the experiences
    names = data["run_name"].astype(str).unique()

    # Retrieving the name of the 'reference'
    ## Step 1: Get the annot names
    cut_names = [name.split("_") for name in names]
    all_annot = max(cut_names, key=len)
    
    ## Step 2: Find the element not used in any composite
    ref_candidates = [item for item in names if item not in all_annot and '_' not in item]

    ## Step 3: Assume there's only one true 'ref'
    ref = ref_candidates[0] if ref_candidates else None
        
    # Getting the number of steps and choosing a maximmum number of columns in the GridSpec
    nb_steps = data["nb_annot"].dropna().nunique()
    ncols_max = 11 if nb_steps <= 7 else int((((11 * (nb_steps - 6)) // 2) * 2) + 1) # Increasing the maximum number of columns and ensuring it is odd (for centering purposes) 

    # Getting the number of rows and columns for the GridSpec
    nrows = data["nb_annot"].dropna().nunique()
    grouped = data.groupby("nb_annot")["run_name"].nunique() # Group by 'nb_annot' and count distinct 'name' values in each group
    exceeding = [nb_expe for nb_expe in grouped if nb_expe > ncols_max]
    for nb_expe in exceeding :
        nrows += nb_expe // ncols_max

    ncols = data['nb_annot'].value_counts().max()
    ncols = ncols if ncols < ncols_max else ncols_max
    
            
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)
    
    plots_path = f"{analysis_path}/{output_file}"

    if os.path.exists(plots_path) :
        os.remove(plots_path)
       
    with PdfPages(plots_path, keep_empty=False) as pdf:
        gs = GridSpec(nrows=nrows, ncols=ncols)
        f = plt.figure(clear=True, figsize=(60*(ncols_max//11), nrows*10))
        descrip = wrap_text(descrip, width=(55 - nrows // 3), new_sep=" ")
        f.suptitle(f"\n{descrip}", fontsize=64+(nrows*3))
        f.subplots_adjust(hspace=.5, top=0.92, bottom=0.05)
                
        leg, new_page = False, False
        name_wdth = max(5, 10 - nrows // 3)

        for score in score_types :
                    
            _data = data[data["score_type"] == score]

            ymin = _data["values"].min()
            ymin = 1.1 * ymin if  ymin < 0 else .9 * ymin
            ymax = _data["values"].max()
            ymax = 1.1 * ymax if  ymax > 0 else .9 * ymax
            

            if score == "insulation_count":
                palette_str = "ch:s=2.75, rot=.1, dark=.7, light=.3, hue=.5, gamma=.9"
            elif score == "PC1" :
                palette_str = "ch:s=1.75, rot=.1, dark=.7, light=.3, hue=.5, gamma=.9"
            else :
                palette_str = "ch:s=0, rot=3, dark=.6, light=.4, hue=2, gamma=1.75"

            l_steps = sorted(_data["nb_annot"].unique(), reverse=True)
            n_colors = len(l_steps)
            palette = color_palette(palette=palette_str, n_colors=n_colors)

            color_map = dict(zip(l_steps, palette)) # Keys are integers
            
            if new_page :
                leg = False
                f = plt.figure(clear=True, figsize=(60, nrows*10))
                f.suptitle(f"\n{descrip}", fontsize=64+(nrows*3))
                f.subplots_adjust(hspace=.5, top=0.92, bottom=0.05)

            ax_i =0
            for i in range(nb_steps) :

                if i == 0:
                    data_i = _data[_data["run_name"] == ref]
                    names_i = ["Reference"]
                    pos = ncols // 2
                    ax = f.add_subplot(gs[ax_i, pos])
                    ax_i += 1
                
                else :
                    data_i = _data[_data["nb_annot"] == i]

                    names_i = data_i["run_name"].astype(str).unique()
                    k = len(names_i)
                    
                    if k > ncols_max :
                        for u in range(0, k, ncols_max) :
                            names_i_u = names_i[u : u+ncols_max]
                            data_i_u = data_i[data_i["run_name"].isin(names_i_u)]
                            nb_iter = k // ncols_max
                            curr_iter = 0
                            
                            if u+ncols_max < k :
                                ax = f.add_subplot(gs[ax_i, 0:ncols])
                                ax_i +=1
                                names_i_u = [wrap_text(name, width=name_wdth, sep="_", new_sep="-") for name in names_i_u]
                                delim = [v +.5 for v in range(0, len(names_i_u))]
                                pos = ax.get_position()
                                y_space = .01 * (curr_iter - (nb_iter // 2))
                                new_pos = [pos.x0, pos.y0 + y_space, pos.width, pos.height]
                                ax.set_position(new_pos)

                                color = color_map[i]

                                boxenplot(data=data_i_u, x="run_name", y="values", color=color, ax=ax, width_method="linear")

                                for d in delim:
                                    ax.axvline(x=d, color='gray', linestyle='--', linewidth=2)
                                ax.tick_params(axis='both', labelsize=30+(nrows*3))
                                ticks = ax.get_xticks()
                                ax.set_xticks(ticks=ticks, labels=names_i_u)
                                ax.set_xlabel("")
                                score = "IS" if score == "insulation_count" else score
                                ax.set_ylabel(f"{score} deviation", fontsize=32+(nrows*3))
                                ax.set_ylim(ymin, ymax)
                                curr_iter += 1
                            
                            elif u+ncols_max == k :
                                names_i = names_i[u : u+ncols_max]
                                data_i = data_i[data_i["run_name"].isin(names_i)]
                                ax = f.add_subplot(gs[ax_i, 0:ncols])
                                ax_i += 1
                                pos = ax.get_position()
                                y_space = .01 * (curr_iter - (nb_iter // 2))
                                new_pos = [pos.x0, pos.y0 + y_space, pos.width, pos.height]
                                ax.set_position(new_pos)
                                

                            else :
                                names_i = names_i[u : u+ncols_max]
                                data_i = data_i[data_i["run_name"].isin(names_i)]
                                l_offset = (ncols - len(names_i)) // 2
                                r_offset = l_offset + len(names_i)
                                ax = f.add_subplot(gs[ax_i, l_offset:r_offset])
                                ax_i += 1
                                x_space = .033 if len(names_i) % 2 == 0 else 0
                                pos = ax.get_position()
                                coef = (curr_iter - (nb_iter // 2))
                                y_space = .01 if coef == 0 else .01 * coef
                                new_pos = [pos.x0 + x_space, pos.y0 + y_space, pos.width, pos.height]
                                ax.set_position(new_pos)
                    
                    elif k == ncols_max :
                        ax = f.add_subplot(gs[ax_i, 0:ncols])
                        ax_i += 1
                    
                    elif k == 0 :
                        logging.info(f"Apparently there is no experimental data ({names_i}) for this " \
                                     f"number of ({i}) annotations")
                        continue

                    else :
                        l_offset = (ncols - k) // 2
                        r_offset = l_offset + k
                        ax = f.add_subplot(gs[ax_i, l_offset:r_offset])
                        ax_i += 1
                        if ((len(names_i) % 2 == 0) and (ncols % 2 != 0)
                            or
                            (len(names_i) % 2 != 0) and (ncols % 2 == 0)) :
                            pos = ax.get_position()
                            x_space = .033 * (ncols / 11)
                            new_pos = [pos.x0 + x_space, pos.y0, pos.width, pos.height]
                            ax.set_position(new_pos)

                names_i = [wrap_text(name, width=name_wdth, sep="_", new_sep="-") for name in names_i] if i!=0 else names_i
                delim = [v +.5 for v in range(0, len(names_i))]

                color = color_map[i]

                boxenplot(data=data_i, x="run_name", y="values", color=color, ax=ax, width_method="linear")

                for d in delim:
                    ax.axvline(x=d, color='gray', linestyle='--', linewidth=2)
                ax.tick_params(axis='both', labelsize=30+(nrows*3))
                ticks = ax.get_xticks()
                ax.set_xticks(ticks=ticks, labels=names_i)
                ax.set_xlabel("")
                score = "IS" if score == "insulation_count" else score
                ax.set_ylabel(f"{score} deviation", fontsize=32+(nrows*3))
                ax.set_ylim(ymin, ymax)

                
                if not leg :
                    if expe_name_dict is None :
                        expe_name_dict = {}
                        for pseudo in _data[_data["nb_annot"] == 1]["run_name"].astype(str).unique() :
                            expe_name_dict[pseudo] = _data[_data["run_name"] == pseudo]["name"].astype(str).unique()[0]

                    legend = "\n".join([f"{pseudo} : {name}" for pseudo, name in expe_name_dict.items()])
                    if ncols > 4 :
                        ax = f.add_subplot(gs[0, ncols-2])
                        x = .5
                    elif ncols > 2 :
                        ax = f.add_subplot(gs[0, ncols-1])
                        x = - .5
                    
                    ax.text(x, .5, legend, ha="left", va="center",
                    transform=ax.transAxes, fontsize=35+(nrows*3), verticalalignment='top',
                    bbox=dict(boxstyle="round", facecolor="white"))
                    ax.axis("off")

                    leg = True
            
            pdf.savefig()
            new_page = True
        
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
    
    parser.add_argument("--descrip",
                        required=True, help="The description that will be used in the title inside of the output pdf.")
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
    
    
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()

    score_types = args.score_types.split(",")
    
    expe_names = None
    if args.expe_names is not None :
        expe_names = [elt.split(":") for elt in args.expe_names.split(",")]
        expe_names = {elt[0]: elt[1] for elt in expe_names}
    
    
    boxplot_slide(descrip=args.descrip, 
                  data_file=args.data_file, 
                  analysis_path=args.analysis_path, 
                  output_file=args.output_file, 
                  score_types=score_types, 
                  create_data=args.create_data, 
                  rename=args.rename, 
                  resol=args.resol,
                  wdir=args.wdir,
                  expe_name_dict=expe_names)


if __name__ == '__main__':
    main()
