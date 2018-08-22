# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import imp
import math
import stat
import click
import zipfile
import requests
import platform
from tqdm import tqdm
from functools import partial

green_fg = partial(click.style, fg='green')
yellow_fg = partial(click.style, fg='yellow')
magenta_fg = partial(click.style, fg='magenta')
cyan_fg = partial(click.style, fg='cyan')
red_fg = partial(click.style, fg='red')


def check_and_create_folders(directory, folder_list, interactive=False):
    for folder in folder_list:
        if os.path.exists(os.path.join(directory, folder)):
            click.echo(red_fg(">>> WARNING: Folder (%s) already exists in path (%s). "
                              "Existing files will be overwritten and garbled!" % (folder.upper(), directory)))
            if interactive:
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
    click.echo(magenta_fg(">>> Attempting to download data..."))
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
