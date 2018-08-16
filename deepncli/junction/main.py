import os
import sys
import imp
import time
import click
import subprocess
from tqdm import tqdm
from functools import partial
from sys import platform as _platform
import joblib.parallel as parallel
from collections import Counter
import warnings
warnings.filterwarnings("ignore")
# project imports
from ..db.junctiondb import JunctionsDatabase, Gene, Junction, Stats
from ..utils.io import get_sam_filelist, get_file_list, make_fasta_file, concatenate_dicts
from proteinprocessor import ProteinProcessor as PProcessor

green_fg = partial(click.style, fg='green')
yellow_fg = partial(click.style, fg='yellow')
magenta_fg = partial(click.style, fg='magenta')
cyan_fg = partial(click.style, fg='cyan')
red_fg = partial(click.style, fg='red')

file_read_progress = {}


def make_search_junctions(junctions_array):
    junction_sequences = []
    for junc in junctions_array:
        junction_sequences.append(junc[30:50])
        junction_sequences.append(junc[26:46])
        junction_sequences.append(junc[22:42])
    return junction_sequences


def junctions_in_read(read, junction_sequences):
    match_index = -1
    junction_index = -1
    for i, j in enumerate(junction_sequences):
        if j in read:
            match_index = read.index(j)
            junction_index = i
    return junction_index, match_index


def search_for_junctions(filepath, jseqs, exclusion_sequence, output_filehandle, f, input_file_size):
    hits_count = 0
    processor = PProcessor()

    def check_matching_criteria(l, indexes, jseqs, read):
        value = 0
        if indexes[0] != -1:
            junction = jseqs[indexes[0]]
            downstream_rf = read[len(junction) + indexes[1] + (indexes[0] % 3) * 4:]
            if len(downstream_rf) > 25:
                if exclusion_sequence not in downstream_rf or exclusion_sequence == '':
                    value = 1
                    protein_sequence = processor.translate_orf(downstream_rf)
                    output_filehandle.write(" ".join(l[:4]) + " " + l[9] + " " + downstream_rf
                                            + " " + protein_sequence + "\n")
        return value

    input_filehandle = open(filepath)
    for line in tqdm(input_filehandle.readlines(), unit=' lines'):
        line_split = line.strip().split()
        if line_split[0][0] != "@" and line_split[2] == "*":
            sequence_read = line_split[9]
            rev_sequence_read = processor.reverse_complement(sequence_read)
            fwd_indexes = junctions_in_read(sequence_read, jseqs)
            hit = check_matching_criteria(line_split, fwd_indexes, jseqs, sequence_read)
            if hit == 0:
                rev_indexes = junctions_in_read(rev_sequence_read, jseqs)
                hit = check_matching_criteria(line_split, rev_indexes, jseqs, rev_sequence_read)
            hits_count += hit
    input_filehandle.close()


def multi_convert(directory, infolder, outfolder):
    file_list = get_file_list(directory, infolder, ".txt")
    for f in file_list:
        make_fasta_file(os.path.join(directory, infolder, f), os.path.join(directory, outfolder, f[:-4] + ".fa"))


def jsearch(directory, filename, input_data_folder, junction_folder, junction_sequence, exclusion_sequence):
    exclusion_sequence = exclusion_sequence.upper() if exclusion_sequence else ""
    click.echo(green_fg('>>> Searching junctions in file: %s' % filename))
    start = time.time()
    filepath = os.path.join(directory, input_data_folder, filename)

    input_file_size = os.path.getsize(filepath)
    output_file_handle = open(os.path.join(directory, junction_folder, filename.replace(".sam", '.junctions.txt')), 'w')
    search_for_junctions(filepath, junction_sequence, exclusion_sequence,
                         output_file_handle, filename, input_file_size)
    output_file_handle.close()
    finish = time.time()
    elapsed = finish - start
    click.echo(cyan_fg("Finished searching junctions in %s in time %f seconds" % (filename, elapsed)))


def junction_search(directory, junction_folder, input_data_folder, blast_results_folder,
                    junction_sequence, exclusion_sequence):
    unmap_files = get_sam_filelist(directory, input_data_folder)
    if not len(unmap_files):
        click.echo(red_fg(">>> ERROR: No .sam files found in directory %s." % directory))
        sys.exit(1)
    junction_seqs = make_search_junctions(junction_sequence)
    click.echo(cyan_fg(">>> The primary, secondary, and tertiary sequences searched are:"))
    for j in junction_seqs:
        click.echo(yellow_fg("    %s" % j))
    cores = parallel.cpu_count()
    click.echo(cyan_fg('>>> Starting junction search on %s cores.' % cores))
    parallel.Parallel(n_jobs=cores)(parallel.delayed(jsearch)(directory, f, input_data_folder, junction_folder,
                                                              junction_seqs, exclusion_sequence) for f in unmap_files)
    multi_convert(directory, junction_folder, blast_results_folder)


def blast_search(directory, db_name, blast_results_folder):
    suffix = ''
    if _platform.startswith('win'):
        suffix = '.exe'
    blast_path = os.path.join(imp.find_module("deepncli")[1], "data", "blast")
    db_path = os.path.join(imp.find_module("deepncli")[1], db_name)
    click.echo(green_fg(">>> Selected Blast DB: %s" % db_name))
    file_list = get_file_list(directory, blast_results_folder, ".fa")
    for file_name in file_list:
        start = time.time()
        output_file = os.path.join(directory, blast_results_folder, file_name.replace(".junctions.fa", '.blast.txt'))
        click.echo(cyan_fg(">>> Running BLAST search for file: " + file_name))
        blast_command_list = [os.path.join(blast_path, 'blastn' + suffix),
                              '-query', os.path.join(directory, blast_results_folder, file_name), '-db', db_path,
                              '-task', 'blastn', '-dust', 'no', '-num_threads', str(parallel.cpu_count()),
                              '-outfmt', '7', '-out', output_file, '-evalue', '0.2', '-max_target_seqs', '10']
        blast_pipe = subprocess.Popen(blast_command_list, shell=False)
        blast_pipe.wait()
        finish = time.time()
        elapsed = finish - start
        click.echo(cyan_fg("Finished blasting file %s in time %f seconds" % (file_name, elapsed)))


def create_gene_list(gene_list_path):
    fh = open(gene_list_path, "r")
    for line in fh.readlines():
        split = line.split()
        if not Gene.select().where(Gene.gene_name == split[1], Gene.orf_start == int(split[6]) + 1,
                                   Gene.nm_number == split[0], Gene.orf_stop == int(split[7])):
            Gene.insert(gene_name=split[1], orf_start=int(split[6]) + 1, orf_stop=int(split[7]),
                        mrna=split[9], intron=split[8], chromosome=split[2], nm_number=split[0]).execute()


def generate_stats(blast_count):
    gene_query = Gene.select()
    for gene in gene_query:
        junc_query = Junction.select().where(Junction.gene == gene)
        frames = []
        orfs = []
        inframe_inorfs = []
        if len(junc_query) >= 1:
            for j in junc_query.dicts():
                frames.append(j['frame'])
                orfs.append(j['orf'])
                if j['inframe_inorf']:
                    inframe_inorfs.append('inframe_inorf')
            insert_dictionary = concatenate_dicts([Counter(frames), Counter(orfs), Counter(inframe_inorfs)])
            insert_dictionary['gene'] = gene
            insert_dictionary['total'] = len(junc_query)
            Stats.insert(**insert_dictionary).execute()
            for j in junc_query:
                j.ppm = j.count * 1000000.0 / blast_count
                j.save()


def _parse_blast_results(directory, blast_results_folder, blasttxt, blast_results_query_folder, gene_list_file):
    start = time.time()
    click.echo(magenta_fg(">>> Parsing blast results for file %s" % blasttxt))
    blast_parsed_results_filepath = os.path.join(directory, blast_results_query_folder,
                                                 blasttxt.replace(".blast.txt", ".db"))
    gene_list_path = os.path.join(imp.find_module("deepncli")[1], gene_list_file)
    # Initialize database for storage
    jdb = JunctionsDatabase(blast_parsed_results_filepath)
    jdb.create_tables()
    # Populate gene table
    create_gene_list(gene_list_path)
    blast_results_handle = open(os.path.join(directory, blast_results_folder, blasttxt), 'r')
    previous_bitscore = 0
    blast_count = 1
    collect_results = True
    click.echo(yellow_fg(">>> Counting blast hits for file %s ..." % blasttxt))
    for line in tqdm(blast_results_handle.readlines(), unit=' lines'):
        line.strip()
        split = line.split()
        if "BLASTN" in line:
            previous_bitscore = 0
            collect_results = True
            blast_count += 1
        if "hits" in line and int(split[1]) > 100:
            collect_results = False

        if split[0] != '#' and collect_results and float(split[2]) > 98 and float(split[11]) > 50.0 and \
            float(split[11]) > previous_bitscore:
            previous_bitscore = float(split[11]) * 0.98
            nm_number = split[1]
            gene_record = Gene.select().where(Gene.nm_number == nm_number)
            if not len(gene_record):
                click.echo(red_fg(">>> ERROR: NM number %s not found in the gene list." % nm_number))
                raise KeyError("%s not found." % nm_number)
            gene = gene_record[0]
            position = int(split[8])
            query_start = int(split[6])
            fudge_factor = query_start - 1
            _frame = position - gene.orf_start - fudge_factor
            # Frame Calculation
            frame = "not_in_frame"
            if _frame % 3 == 0 or _frame == 0:
                frame = "in_frame"
            if gene.intron == "INTRON":
                frame = "intron"
            if int(split[9]) - position < 0:
                frame = "backwards"
            # Orf calculation
            orf = "in_orf"
            if position < gene.orf_start:
                orf = "upstream"
            if position > gene.orf_stop:
                orf = "downstream"
            # Frame Orf calculation
            inframe_inorf = False
            if frame == 'in_frame' and orf == 'in_orf':
                inframe_inorf = True

            junction_record = Junction.select().where(Junction.gene == gene, Junction.position == position,
                                                      Junction.query_start == query_start)
            if not len(junction_record):
                Junction.insert(gene=gene, position=position, query_start=query_start,
                                frame=frame, orf=orf, ppm=0.0, inframe_inorf=inframe_inorf, count=1).execute()
            else:
                junction_record[0].count += 1
                junction_record[0].save()
    click.echo(green_fg(">>> Generating stats for file %s ..." % blasttxt))
    generate_stats(blast_count)
    finish = time.time()
    elapsed = finish - start
    click.echo(cyan_fg("Finished parsing file %s in time %f seconds" % (blasttxt, elapsed)))
    jdb.close_db()


def parse_blast_results(directory, blast_results_folder, blast_results_query_folder, gene_list_file):
    blast_results_list = get_file_list(directory, blast_results_folder, ".txt")
    cores = parallel.cpu_count()
    click.echo(cyan_fg('>>> Starting parsing blast results on %s cores.' % cores))
    parallel.Parallel(n_jobs=cores)(parallel.delayed(_parse_blast_results)(directory,
                                                                           blast_results_folder, f,
                                                                           blast_results_query_folder,
                                                                           gene_list_file) for f in blast_results_list)

