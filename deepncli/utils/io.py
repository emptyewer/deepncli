# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import click
from functools import partial

yellow_fg = partial(click.style, fg='yellow')
magenta_fg = partial(click.style, fg='magenta')
red_fg = partial(click.style, fg='red')


def check_and_create_folders(directory, folder_list, non_interactive=False):
    for folder in folder_list:
        if os.path.exists(os.path.join(directory, folder)):
            click.echo(red_fg(">>> WARNING: Folder (%s) already exists in path (%s). "
                              "Existing files will be overwritten and garbled!" % (folder.upper(), directory)))
            if not non_interactive:
                if not click.confirm(magenta_fg('Do you want to continue? If you do, '
                                                'the existing files will be overwritten')):
                    click.echo(red_fg("ABORTING..."))
                    sys.exit(1)
        else:
            os.makedirs(os.path.join(directory, folder))


def get_file_list(directory, folder, suffix):
    file_list = os.listdir(os.path.join(directory, folder))
    return_list = []
    for fi in file_list:
        if os.path.splitext(fi)[1] == suffix:
            return_list.append(fi)
    return return_list


def get_sam_filelist(directory, input_folder):
    file_list = []
    if os.path.exists(os.path.join(directory, input_folder)):
        file_list = get_file_list(directory, input_folder, '.sam')
    return file_list


def make_fasta_file(junction_file, output_file):
    counter = 0
    input_file = open(junction_file, 'r')
    output_file = open(output_file, 'w')
    for line in input_file:
        line.strip()
        split = line.split()
        output_file.write(">%s\n%s\n" % (str(split[0]), str(split[5])))
        counter += 1
    input_file.close()
    output_file.close()
    click.echo(magenta_fg('>>> Converted  %d junctions '
                          'in %s to a FASTA file.' % (counter, os.path.split(junction_file)[-1])))


def concatenate_dicts(list):
    output_dict = {}
    for d in list:
        output_dict = dict(output_dict.items() + dict(d).items())
    return output_dict
