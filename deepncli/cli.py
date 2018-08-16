# -*- coding: utf-8 -*-

"""Console script for deepncli."""
import os
import sys
import imp
import math
import click
import requests
import zipfile
import platform
import stat
from tqdm import tqdm
from functools import partial
# project imports
from .utils.io import check_and_create_folders
from .junction.main import junction_search, blast_search, parse_blast_results

green_fg = partial(click.style, fg='green')
yellow_fg = partial(click.style, fg='yellow')
magenta_fg = partial(click.style, fg='magenta')
cyan_fg = partial(click.style, fg='cyan')
red_fg = partial(click.style, fg='red')
deepn_option = partial(click.option, show_default=True)


class Config(object):
    def __init__(self):
        self.verbose = False
        self.debug = False


pass_config = click.make_pass_decorator(Config, ensure=True)

blast_dbs = {'hg38': os.path.join("data", "blastdb", "hg38NMgenes.db"),
             'saccer3': os.path.join("data", "blastdb", "sacCer3.db"),
             'mm10': os.path.join("data", "blastdb", "mm10mRNAnoSuffix.fa"),
             'hg38_pGAD': os.path.join("data", "blastdb", "hg38NMgenes.db"),
             'saccer3_pGAD': os.path.join("data", "blastdb", "sacCer3.db")
             }
gene_lists = {'hg38': os.path.join("data", "lists", "hg38GeneList.prn"),
              'saccer3': os.path.join("data", "lists", "sacCer3GeneList.prn"),
              'mm10': os.path.join("data", "lists", "mm10GeneList.prn"),
              'hg38_pGAD': os.path.join("data", "lists", "hg38GeneList.prn"),
              'saccer3_pGAD': os.path.join("data", "lists", "sacCer3GeneList.prn")
              }
junction_sequences = {'hg38': "CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC",
                      'saccer3': "CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC",
                      'mm10': "AATTCCACCCAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGG",
                      'saccer3_pGAD': "ATGATGAAGATACCCCACCAAACCCAAAAAAAGAGATCGAATTCCCGGGG",
                      'hg38_pGAD': "AATTCCACCCAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGG"}


def download_url(url, output_path):
    file_path = os.path.join(output_path, os.path.basename(url))
    output_dir = os.path.join(output_path, os.path.splitext(os.path.basename(url))[0].split("_")[0])
    if not os.path.exists(output_dir):
        # Streaming, so we can iterate over the response.
        r = requests.get(url, stream=True)
        # Total size in bytes.
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024
        wrote = 0
        click.echo(green_fg(">>> Downloading %s" % os.path.basename(url)))
        with open(file_path, 'wb') as f:
            for data in tqdm(r.iter_content(block_size), total=math.ceil(total_size//block_size), unit='B',
                             unit_scale=True):
                wrote = wrote + len(data)
                f.write(data)
        zip_ref = zipfile.ZipFile(file_path, 'r')
        zip_ref.extractall(output_dir)
        zip_ref.close()
        os.remove(file_path)
        if total_size != 0 and wrote != total_size:
            click.echo(red_fg(">>> ERROR, something went wrong"))


def download_data():
    download_list = ["https://github.com/emptyewer/deepncli/releases/download/support/lists.zip",
                     "https://github.com/emptyewer/deepncli/releases/download/support/blastdb.zip"]
    if platform.system() == "Windows":
        download_list.append("https://github.com/emptyewer/deepncli/releases/download/support/blast_win.zip")
    elif platform.system() == "Linux":
        download_list.append("https://github.com/emptyewer/deepncli/releases/download/support/blast_linux.zip")
    elif platform.system() == "Darwin":
        download_list.append("https://github.com/emptyewer/deepncli/releases/download/support/blast_osx.zip")
    output_dir = os.path.join(imp.find_module("deepncli")[1], "data")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for link in download_list:
        if not os.path.exists(os.path.join(output_dir, os.path.splitext(os.path.basename(link))[0])):
            download_url(link, output_dir)
    if platform.system() != "Windows":
        for dirpath, dirnames, filenames in os.walk(os.path.join(output_dir, "blast")):
            for filename in filenames:
                path = os.path.join(dirpath, filename)
                os.chmod(path, stat.S_IXUSR)


def verify_options(*args, **kwargs):
    if not kwargs['genome'] in blast_dbs.keys():
        click.echo(red_fg(">>> ERROR: Specified option for genome selection (%s) not available" % kwargs['genome']))
        sys.exit(1)
    if not os.path.exists(kwargs['dir']):
        click.echo(red_fg(">>> ERROR: Specified work folder (%s) does not exist." % kwargs['dir']))
        sys.exit(1)


@click.group()
@deepn_option("--verbose", is_flag=True, default=False)
@deepn_option("--debug", is_flag=True, default=False)
@pass_config
def main(config, verbose, debug):
    """Console script for deepncli."""
    click.echo(magenta_fg("\n{}  DEEPN  {}\n".format("*"*10, "*"*10)))
    config.verbose = verbose  # pragma: no cover
    config.debug = debug  # pragma: no cover


@main.command()
@deepn_option("--dir", required=True, help="path to work folder")
@pass_config
def initialize(config, *args, **kwargs):
    if not os.path.exists(kwargs['dir']):
        os.makedirs(kwargs['dir'])
    else:
        click.echo(red_fg(">>> ERROR: Specified work folder (%s) already exists." % kwargs['dir']))
        sys.exit(1)
    check_and_create_folders(kwargs['dir'], ['mapped_sam_files', 'unmapped_sam_files', 'sam_files'])
    click.echo(green_fg(">>> Sucessfully initiated work folder."))


@main.command()
@deepn_option("--dir", required=True, help="path to work folder")
@deepn_option("--genome", required=True, help="name of the reference organism. "
                                              "options: mm10/hg38/saccer3/hg38_pGAD,saccer3_pGAD")
@deepn_option("--seq", required=False, default="", help="junction sequence used in sequencing. "
                                                        "If provided default junction sequence will be ignored")
@deepn_option("--exclude_seq", required=False, default="", help="sequence to exclude from junction matching")
@deepn_option("--unmapped", is_flag=True, default=False, help="if flag is enabled, .sam files will be read from "
                                                              "unmapped_sam_files folder")
@deepn_option("--non_interactive", is_flag=True, default=False, help="if enabled interactive session will be turned off"
                                                                     ". (careful, files maybe be overwritten")
@pass_config
def junction_make(config, *args, **kwargs):
    click.echo(green_fg("\n{}  Junction Make  {}\n".format(">" * 10, "<" * 10)))
    input_data_folder = 'unmapped_sam_files'  # Manage name of input folder here
    junction_folder = 'junction_files'  # Manage name of junction reads output folder here
    blast_results_folder = 'blast_results'  # Manage name of blast results output folder here
    blast_results_query = 'blast_results_query'  # Manage name of blast results dictionary output folder here
    junction_sequence = junction_sequences[kwargs['genome']].replace(" ", "").split(",")
    if kwargs['seq']:
        junction_sequence = kwargs['seq'].replace(" ", "").split(",")
    exclusion_sequence = kwargs['exclude_seq'].replace(" ", "")
    blast_db = blast_dbs[kwargs['genome']]
    gene_list_file = gene_lists[kwargs['genome']]
    # verify if the options provided are valid
    verify_options(*args, **kwargs)
    # create folders for junction make
    check_and_create_folders(kwargs['dir'], ['junction_files', 'blast_results', 'blast_results_query'],
                             non_interactive=kwargs['non_interactive'])
    # search for junctions
    junction_search(kwargs['dir'], junction_folder, input_data_folder, blast_results_folder,
                    junction_sequence, exclusion_sequence)
    # blast the junctions
    blast_search(kwargs['dir'], blast_db, blast_results_folder)
    # parse blast results
    parse_blast_results(kwargs['dir'], blast_results_folder, blast_results_query, gene_list_file)


@main.command()
@pass_config
def gene_count(config, *args, **kwargs):
    pass


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover