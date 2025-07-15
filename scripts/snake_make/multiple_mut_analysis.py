import argparse
import textwrap
import os
import sys
c_path = "/home/fforge/Stage-IA3D/scripts/"
sys.path.append(f"{c_path}/orcanalyse")

import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from seaborn import boxenplot
import re

def natural_key(s):
    # Split the string into a list of strings and integers
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def wrap_text(text: str, width: int, sep: str = " "):
    words = text.split(sep=sep)
    lines = []
    current_line = ""
    for word in words:
        if len(current_line) + len(word) + (1 if current_line else 0) <= width:
            current_line += (" " if current_line else "") + word
        else:
            # Check if adding the word would split it more than halfway
            if len(word) > width:
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


def boxplot_slide(descrip: str, data_file: str, analysis_path:str, output_file: str, score_types: list):
    """
    """
    data = pd.read_csv(data_file, sep="\t")
    # data["values"] = data["values"] ** 2
    
    # Remove '_rdm' suffix if present, using pandas vectorized string operations
    base_names = data["name"].astype(str).str.replace(r"_rdm$", "", regex=True).unique()
    names = sorted(base_names, key=natural_key)

    # Build all_names using a list comprehension
    all_names = [item for name in names for item in (name, f"{name}_rdm")]
            
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)
    
    plots_path = f"{analysis_path}/{output_file}"

    if os.path.exists(plots_path) :
        os.remove(plots_path)
       
    with PdfPages(plots_path, keep_empty=False) as pdf:
        if len(names) >= 11 :
            rest = (len(names) % 10) / 10
            width_ratios = [(1-rest)/2, rest, (1-rest)/2]
        gs = GridSpec(nrows= len(score_types), ncols=1) if len(names) < 11 else GridSpec(nrows=math.ceil(len(names)/10), ncols=3, width_ratios=width_ratios)
        f = plt.figure(clear=True, figsize=(60, 33.75))
        descrip = wrap_text(descrip, width=91)
        f.suptitle(f"\n{descrip}", fontsize=48)

        ticks = [i +.5 for i in range(0, 2*len(names), 2)]
        delim = [i +.5 for i in range(1, 2*len(names), 2)]

        if len(names) < 11 :
            name_wdth = 230//len(names) if len(names) < 11 else 30
            names = [wrap_text(name, width=name_wdth, sep="_") for name in names]

            handles, labels, leg = None, None, False
            for i, score in enumerate(score_types) :
                _data = data[data["score_type"] == score]
                ax = f.add_subplot(gs[i, 0])

                boxenplot(data=_data, x="name", y="values", hue="hue_name", ax=ax, width_method="linear")

                for d in delim:
                    ax.axvline(x=d, color='gray', linestyle='--', linewidth=2)
                ax.tick_params(axis='both', labelsize=24)
                ax.set_xticks(ticks=ticks, labels=names)
                ax.set_xlabel("")
                ax.set_ylabel("Insulation score deviation", fontsize=26)

                # Extract legend info
                if handles is None and labels is None :
                    handles, labels = ax.get_legend_handles_labels()

                # Remove the default legend
                if ax.get_legend() is not None:
                    ax.get_legend().remove()

                # Place the legend outside the plot
                if not leg :
                    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(-.1, -.1), title="", fontsize=26)
                    leg = True

            pdf.savefig()
            pdf.close()
        
        else :
            handles, labels, leg, new_page = None, None, False, False
            name_wdth = 20

            for score in score_types :
                _data = data[data["score_type"] == score]
                
                if new_page :
                    handles, labels, leg = None, None, False
                    f = plt.figure(clear=True, figsize=(60, 33.75))
                    f.suptitle(f"\n{descrip}", fontsize=48)

                for ax_i, i in enumerate(range(0, len(names), 10)) :
                    ymin = _data["values"].min()
                    ymin = 1.1 * ymin if  ymin < 0 else .9 * ymin
                    ymax = _data["values"].max()
                    ymax = 1.1 * ymax if  ymax > 0 else .9 * ymax
                    
                    all_names_i = all_names[i*2 : i*2+20]
                    names_i = names[i : i+10]
                    data_i = _data[_data["name"].isin(all_names_i)]
                    ax = f.add_subplot(gs[ax_i, 0:3]) if len(names_i) == 10 else f.add_subplot(gs[ax_i, 1])

                    names_i = [wrap_text(name, width=name_wdth, sep="_") for name in names_i]
                    ticks = [i +.5 for i in range(0, 2*len(names_i), 2)]
                    delim = [i +.5 for i in range(1, 2*len(names_i), 2)]

                    boxenplot(data=data_i, x="name", y="values", hue="hue_name", ax=ax, width_method="linear")

                    for d in delim:
                        ax.axvline(x=d, color='gray', linestyle='--', linewidth=2)
                    ax.tick_params(axis='both', labelsize=24)
                    ax.set_xticks(ticks=ticks, labels=names_i)
                    ax.set_xlabel("")
                    ax.set_ylabel("Insulation score deviation", fontsize=26)
                    ax.set_ylim(ymin, ymax)

                    # Extract legend info
                    if handles is None and labels is None :
                        handles, labels = ax.get_legend_handles_labels()
                        if labels == ["orcarun_Wtd_mut", "merged_rdm"] :
                            labels = ["Wanted mutation", "Merged random"]


                    # Remove the default legend
                    if ax.get_legend() is not None:
                        ax.get_legend().remove()

                    # Place the legend outside the plot
                    if not leg :
                        ax.legend(handles, labels, loc='center left', bbox_to_anchor=(-.14, -.66), title="", fontsize=26)
                        leg =True

                pdf.savefig()
                new_page = True
            
            pdf.close()




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
    parser.add_argument("--score_types",
                        required=True, help="The list of scores to study (format : 'score1,score2,...')")
    
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_arguments()

    score_types = args.score_types.split(",")

    boxplot_slide(descrip=args.descrip, 
                   data_file=args.data_file, 
                   analysis_path=args.analysis_path, 
                   score_types=score_types)