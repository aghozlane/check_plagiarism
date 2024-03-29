#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Check plagiarism."""

# pycode_similar
import argparse
import sys
import numpy as np
import subprocess
import re
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
import pickle
import pandas as pd
from pathlib import Path


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def isdir(path: str) -> Path:  # pragma: no cover
    """Check if path can be valid directory.

    :param path: Path to the directory

    :raises ArgumentTypeError: If directory does not exist

    :return: (str) Path object of the directory
    """
    mydir = Path(path)
    # if not mydir.is_dir():
    if mydir.is_file():
        msg = f"{mydir.name} is a file."
        raise argparse.ArgumentTypeError(msg)
        # else:
        #     msg = f"{mydir.name} does not exist."
    return mydir


def get_arguments():
    """Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-p",
        dest="project",
        type=str,
        choices=["debruijn", "agc", "gpred"],
        nargs="+",
        default="debruijn",
        help="Which program to check.",
    ),
    parser.add_argument(
        "-i",
        dest="input_dir",
        type=isdir,
        default=Path.cwd(),
        help="Input directory containing all students files (default current directory).",
    ),
    parser.add_argument(
        "-c",
        dest="cluster",
        action="store_true",
        default=False,
        help="Clusterize the heatmap.",
    )
    parser.add_argument(
        "-s",
        dest="saveobj",
        type=isfile,
        help="Saved comparison" "debruijn/debruijn.py for instance",
    )
    parser.add_argument(
        "-o", dest="output_dir", type=Path, default="plagirism", help="Output directory"
    )
    return parser.parse_args()


def get_debruijn(pattern, directory):
    """Find all code"""
    to_compare = []
    students_list = []
    for direct in directory.iterdir():
        debruijn = list(direct.glob(f"**/{pattern}"))
        if len(debruijn) == 1:
            assert debruijn[0].exists()
            to_compare += [[direct, debruijn[0]]]
            students_list += [direct]
        else:
            print(f"Directory ignored: {direct}")
    return to_compare, students_list


def run_command(cmd):
    """Run command
    Arguments:
        cmd: Command to run
    """
    res = {}
    regex_global = re.compile(rb"([0-9]+\.[0-9]+)\s\%.+")
    regex_func = re.compile(rb"([0-9]+\.[0-9]+)\s*\: ref (\S+)<.+")
    # regex_func = re.compile(r"[0-9]+\.[0-9]+ \: ref (\S+)\s[0-9]+:[0-9]+\s\, candidate isfile\s[0-9]+:[0-9]+\s")
    try:
        p = subprocess.Popen(
            cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )
        for line in iter(p.stdout.readline, b""):
            global_match = regex_global.match(line)
            func_match = regex_func.match(line)
            if global_match:
                res["global"] = float(global_match.group(1))
            elif func_match:
                res[func_match.group(2).decode("ascii")] = (
                    float(func_match.group(1)) * 100
                )
            # else:
            #    print("no match: {}".format(line))
        # Cas erreur
        assert len(p.stderr.readlines()) == 0
    except OSError as e:
        sys.exit("Execution failed: {0}".format(e))
    # except AssertionError:
    #     sys.exit("The command failed: {}".format(cmd))
    # except:
    #    sys.exit("There is something wrong with the command: {0}".format(cmd))
    # print(res)
    return res


def get_matrix(matrix_dat, category):
    # print(matrix_dat)
    matrix_dist = np.zeros((len(matrix_dat[0]), len(matrix_dat[0])))
    for i in range(len(matrix_dat[0])):
        # print(matrix_dat[i])
        for j in range(len(matrix_dat[0])):
            if category in matrix_dat[i + 1][j]:
                matrix_dist[i, j] = matrix_dat[i + 1][j][category]
            else:
                print(matrix_dat[0][i])
                print(matrix_dat[0][j])
                print(matrix_dat[i + 1][j])
                sys.exit()
    return pd.DataFrame(matrix_dist, columns=matrix_dat[0], index=matrix_dat[0])


def plot_heatmap(matrix_dist, output_file):
    """Plot half heatmap"""
    fig, ax = plt.subplots()
    sns_plot = sns.heatmap(
        matrix_dist, xticklabels=True, yticklabels=True, linewidths=0.5, cmap="YlGnBu"
    )
    figure = sns_plot.get_figure()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(10)
    plt.tight_layout()
    figure.savefig(output_file, dpi=900)


def plot_clustered_heatmap(matrix_dist, output_file):
    """Plot half heatmap"""
    sns_plot = sns.clustermap(
        matrix_dist, col_cluster=True, linewidths=0.5, cmap="YlGnBu"
    )
    sns_plot.savefig(output_file, dpi=600)


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    args = get_arguments()
    if "debruijn" in args.project:
        pattern = "debruijn/debruijn.py"
        totest = [
            "global",
            "cut_kmer",
            "read_fastq",
            "build_graph",
            "build_kmer_dict",
            "simplify_bubbles",
            "solve_bubble",
            "get_starting_nodes",
            "get_sink_nodes",
            "get_contigs",
            "save_contigs",
            "path_average_weight",
            "remove_paths",
        ]
    elif "agc" in args.project:
        pattern = "agc/agc.py"
        totest = [
            "global",
            "read_fasta",
            "dereplication_fulllength",
            "get_identity",
            "abundance_greedy_clustering",
            "write_OTU",
        ]
        # "get_chunks", "cut_kmer", "get_unique_kmer", "search_mates", "detect_chimera", "chimera_removal",
    elif "gpred" in args.project:
        pattern = "gpred/gpred.py"
        totest = [
            "global",
            "read_fasta",
            "find_start",
            "find_stop",
            "predict_genes",
            "has_shine_dalgarno",
            "predict_genes",
        ]
    else:
        sys.exit("I don't know this program")

    if not args.saveobj:
        to_compare, student_list = get_debruijn(pattern, args.input_dir)
        nb_student = len(student_list)
        matrix_dat = [student_list]
        cmd = "pycode_similar -p 0 {} {}"
        for i in range(nb_student):
            one_compair = []
            for j in range(nb_student):
                one_compair += [
                    run_command(cmd.format(to_compare[i][1], to_compare[j][1]))
                ]
            matrix_dat += [one_compair]
        # with open(os.getcwd() + "/student_compare.pkl", "wb") as compair:
        #     pickle.dump(matrix_dat, compair)
    else:
        with args.saveobj.open("rb") as compair:
            matrix_dat = pickle.load(compair)
    # Function
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # "read_fasta",
    if args.cluster:
        print("Plot clustered heatmaps")
        output = args.output_dir / "plagiarism_clust_heatmap"
        hip = plot_clustered_heatmap
    else:
        print("Plot heatmaps")
        output = args.output_dir / "plagiarism_heatmap"
        hip = plot_heatmap
    for func in totest:
        print("Function: " + func)
        matrix_dist = get_matrix(matrix_dat, func)
        print(matrix_dist.mean().sort_values(ascending=False))
        # plot global
        hip(matrix_dist, f"{output}_{func}.png")
    print("done")
    # "global"


if __name__ == "__main__":
    main()
