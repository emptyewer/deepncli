# -*- coding: utf-8 -*-

"""Console script for deepncli."""
# project imports
from .utils.io import check_and_create_folders
from .utils.download import download_data
import joblib.parallel as parallel
from .junction.main import junction_search, blast_search, parse_blast_results
# from .genecount.main import count_genes
# Library imports
import os
import sys
import click
from functools import partial
import warnings
warnings.filterwarnings("ignore")

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

exon_dictionary = {'hg38': os.path.join("data", "exons", "hg38exonDict.p"),
                   'saccer3': os.path.join("data", "exons", "sacCer3GeneDict.p"),
                   'mm10': os.path.join("data", "exons", "mm10GeneDict.p"),
                   'hg38_pGAD': os.path.join("data", "exons", "hg38exonDict.p"),
                   'saccer3_pGAD': os.path.join("data", "liexonssts", "sacCer3GeneDict.p")
                   }

junction_sequences = {'hg38': "CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC",
                      'saccer3': "CCTCTGCGAGTGGTGGCAACTCTGTGGCCGGCCCAGCCGGCCATGTCAGC",
                      'mm10': "AATTCCACCCAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGG",
                      'saccer3_pGAD': "ATGATGAAGATACCCCACCAAACCCAAAAAAAGAGATCGAATTCCCGGGG",
                      'hg38_pGAD': "AATTCCACCCAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGG"}


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
    download_data()
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
                                              "options: mm10/hg38/saccer3/hg38_pGAD/saccer3_pGAD")
@deepn_option("--seq", required=False, default="", help="junction sequence used in sequencing. "
                                                        "If provided default junction sequence will be ignored")
@deepn_option("--threads", required=False, help="Number of threads to use for processing the files. "
                                                "Defaults to the number of processors.")
@deepn_option("--exclude_seq", required=False, default="", help="sequence to exclude from junction matching")
@deepn_option("--unmapped", is_flag=True, help="if flag is enabled, .sam files will "
                                               "be read from unmapped_sam_files folder")
@deepn_option("--interactive", is_flag=True, help="if enabled interactive session will be turned on.")
@pass_config
def junction_make(config, *args, **kwargs):
    click.echo(green_fg("\n{}  Junction Make  {}\n".format(">" * 10, "<" * 10)))
    threads = kwargs['threads'] if kwargs['threads'] else parallel.cpu_count()
    input_data_folder = 'unmapped_sam_files' if kwargs['unmapped'] else 'sam_files'
    junction_folder = 'junction_files'  # Manage name of junction reads output folder here
    blast_results_folder = 'blast_results'  # Manage name of blast results output folder here
    blast_results_query = 'blast_results_query'  # Manage name of blast results dictionary output folder here
    junction_sequence = junction_sequences[kwargs['genome']].replace(" ", "").split(",")
    if kwargs['seq'] != "":
        junction_sequence = kwargs['seq'].replace(" ", "").split(",")
    exclusion_sequence = kwargs['exclude_seq'].replace(" ", "")
    blast_db = blast_dbs[kwargs['genome']]
    gene_list_file = gene_lists[kwargs['genome']]
    # verify if the options provided are valid
    verify_options(*args, **kwargs)
    # create folders for junction make
    check_and_create_folders(kwargs['dir'], ['junction_files', 'blast_results', 'blast_results_query'],
                             interactive=kwargs['interactive'])
    if kwargs['interactive']:
        if not click.confirm(magenta_fg('\nDo you want to search junctions and blast?')):
            click.echo(red_fg("...Skipping search junctions and blast..."))
        else:
            # search for junctions
            junction_search(kwargs['dir'], junction_folder, input_data_folder, blast_results_folder,
                            junction_sequence, exclusion_sequence, threads)
            # blast the junctions
            blast_search(kwargs['dir'], blast_db, blast_results_folder)

        if not click.confirm(magenta_fg('\nDo you want to parse blast results')):
            click.echo(red_fg("ABORTING..."))
            sys.exit(1)
        else:
            # parse blast results
            parse_blast_results(kwargs['dir'], blast_results_folder, blast_results_query, gene_list_file, threads)
    else:
        # search for junctions
        junction_search(kwargs['dir'], junction_folder, input_data_folder, blast_results_folder,
                        junction_sequence, exclusion_sequence, threads)
        # blast the junctions
        blast_search(kwargs['dir'], blast_db, blast_results_folder)
        # parse blast results
        parse_blast_results(kwargs['dir'], blast_results_folder, blast_results_query, gene_list_file, threads)


# @main.command()
# @deepn_option("--dir", required=True, help="path to work folder")
# @deepn_option("--genome", required=True, help="name of the reference organism. "
#                                               "options: mm10/hg38/saccer3/hg38_pGAD/saccer3_pGAD")
# @deepn_option("--unmapped", is_flag=True, help="if flag is enabled, .sam files will "
#                                                "be read from unmapped_sam_files folder")
# @pass_config
# def gene_count(config, *args, **kwargs):
#     click.echo(green_fg("\n{}  Gene Count  {}\n".format(">" * 10, "<" * 10)))
#     input_data_folder = 'unmapped_sam_files' if kwargs['unmapped'] else 'sam_files'
#     exon_file = exon_dictionary[kwargs['genome']]
#     summary_folder = 'gene_count_summary'
#     check_and_create_folders(kwargs['dir'], ['chromosome_files', 'gene_count_summary', 'gene_count_indices'],
#                              interactive=kwargs['interactive'])
#     # Count genes
#     count_genes(kwargs['dir'], input_data_folder, summary_folder, exon_file)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
