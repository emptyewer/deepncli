# project imports
from ..db.junctiondb import JunctionsDatabase, Gene, Junction, Stats
from ..utils.io import get_sam_filelist, get_file_list, make_fasta_file, concatenate_dicts, count_lines
from ..utils.time import elapsed_time
from proteinprocessor import ProteinProcessor as PProcessor
# Other imports
import os
import sys
import time
import click
import subprocess
from tqdm import tqdm
from functools import partial
from sys import platform as _platform
import joblib.parallel as parallel
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")


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


def search_for_junctions(filepath, jseqs, exclusion_sequence, output_filehandle):
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
    bar = tqdm(total=count_lines(filepath), unit=' lines')
    for line in input_filehandle:
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
        bar.update(1)
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

    output_file_handle = open(os.path.join(directory, junction_folder, filename.replace(".sam", '.junctions.txt')), 'w')
    search_for_junctions(filepath, junction_sequence, exclusion_sequence,
                         output_file_handle)
    output_file_handle.close()
    finish = time.time()
    hr, min, sec = elapsed_time(start, finish)
    click.echo(cyan_fg("Finished searching junctions in %s in time %d hr, %d min, %d sec" % (filename, hr, min, sec)))


def junction_search(directory, junction_folder, input_data_folder, blast_results_folder,
                    junction_sequence, exclusion_sequence, threads):
    unmap_files = get_sam_filelist(directory, input_data_folder)
    if not len(unmap_files):
        click.echo(red_fg(">>> ERROR: No .sam files found in directory %s." % directory))
        sys.exit(1)
    junction_seqs = make_search_junctions(junction_sequence)
    click.echo(cyan_fg(">>> The primary, secondary, and tertiary sequences searched are:"))
    for j in junction_seqs:
        click.echo(yellow_fg("    %s" % j))
    click.echo(cyan_fg('>>> Starting junction search on %s cores.' % threads))
    parallel.Parallel(n_jobs=threads)(parallel.delayed(jsearch)(directory, f, input_data_folder, junction_folder,
                                                                junction_seqs, exclusion_sequence) for f in unmap_files)
    multi_convert(directory, junction_folder, blast_results_folder)


def blast_search(directory, db_name, blast_results_folder):
    suffix = ''
    if _platform.startswith('win'):
        suffix = '.exe'
    blast_path = os.path.join(os.path.expanduser('~'), ".deepn", "data", "blast")
    db_path = os.path.join(os.path.expanduser('~'), ".deepn", db_name)
    click.echo(green_fg(">>> Selected Blast DB: %s" % db_name))
    file_list = get_file_list(directory, blast_results_folder, ".fa")
    for file_name in file_list:
        if not os.path.getsize(os.path.join(directory, blast_results_folder, file_name)) == 0:
            start = time.time()
            output_file = os.path.join(directory, blast_results_folder, file_name.replace(".junctions.fa", '.blast.txt'))
            click.echo(yellow_fg(">>> Running BLAST search for file: " + file_name))
            blast_command_list = [os.path.join(blast_path, 'blastn' + suffix),
                                  '-query', os.path.join(directory, blast_results_folder, file_name), '-db', db_path,
                                  '-task', 'blastn', '-dust', 'no', '-num_threads', str(parallel.cpu_count()),
                                  '-outfmt', '7', '-out', output_file, '-evalue', '0.2', '-max_target_seqs', '10']
            blast_pipe = subprocess.Popen(blast_command_list, shell=False)
            blast_pipe.wait()
            finish = time.time()
            hr, min, sec = elapsed_time(start, finish)
            click.echo(cyan_fg("Finished blasting file %s in time %d hr, %d min, %d sec" % (file_name, hr, min, sec)))
        else:
            click.echo(red_fg(">>> ERROR: File %s does not have any junctions, "
                              "please check if they right genome was chosen." % file_name))
            sys.exit(1)


def create_gene_list(gene_list_path):
    fh = open(gene_list_path, "r")
    nm_gene_dictionary = defaultdict()
    for line in fh:
        split = line.split()
        nm_gene_dictionary[split[0]] = [split[1], int(split[6]) + 1, int(split[7]), split[8]]
        Gene.create(gene_name=split[1], orf_start=int(split[6]) + 1, orf_stop=int(split[7]),
                    mrna=split[9].upper(), intron=split[8], chromosome=split[2], nm_number=split[0])
    return nm_gene_dictionary


def generate_stats(blast_count):
    gene_query = Gene.select()
    for gene in tqdm(gene_query, unit=" genes"):
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
    click.echo(magenta_fg(">>> Reading blast output for file %s" % blasttxt))
    blast_parsed_results_filepath = os.path.join(directory, blast_results_query_folder,
                                                 blasttxt.replace(".blast.txt", ".db"))
    gene_list_path = os.path.join(os.path.expanduser('~'), ".deepn", gene_list_file)
    # Initialize database for storage
    jdb = JunctionsDatabase(blast_parsed_results_filepath)
    jdb.create_tables()
    # Populate gene table
    nm_gene_dictionary = create_gene_list(gene_list_path)
    blast_results_handle = open(os.path.join(directory, blast_results_folder, blasttxt), 'r')
    previous_bitscore = 0
    blast_count = 0
    rejected_count = 0
    accepted_count = 0
    collect_results = True
    click.echo(yellow_fg(">>> Consolidating blast hits for file %s ..." % blasttxt))
    parsed_results = defaultdict(int)
    bar = tqdm(total=count_lines(os.path.join(directory, blast_results_folder, blasttxt)), unit=' lines')
    for line in blast_results_handle:
        line.strip()
        split = line.split()
        if "BLASTN" in line:
            previous_bitscore = 0
            collect_results = True
            blast_count += 1

        elif "hits" in line and int(split[1]) > 100:
            collect_results = False

        elif split[0] != '#' and collect_results and float(split[2]) > 98 and float(split[11]) > 50.0 and float(split[11]) > previous_bitscore:
            accepted_count += 1
            previous_bitscore = float(split[11]) * 0.98
            nm_number = split[1]
            gene = nm_gene_dictionary[nm_number][0]
            position = int(split[8])
            query_start = int(split[6])
            fudge_factor = query_start - 1
            _frame = position - nm_gene_dictionary[nm_number][1] - fudge_factor
            # Frame Calculation
            frame = "not_in_frame"
            if _frame % 3 == 0 or _frame == 0:
                frame = "in_frame"
            if nm_gene_dictionary[nm_number][3] == "INTRON":
                frame = "intron"
            if int(split[9]) - position < 0:
                frame = "backwards"
            # Orf calculation
            orf = "in_orf"
            if position < nm_gene_dictionary[nm_number][1]:
                orf = "upstream"
            if position > nm_gene_dictionary[nm_number][2]:
                orf = "downstream"
            # Frame Orf calculation
            inframe_inorf = False
            if frame == 'in_frame' and orf == 'in_orf':
                inframe_inorf = True
            parsed_results["|".join([gene, nm_number, frame, orf, str(int(inframe_inorf)), str(position),
                                     str(query_start)])] += 1
        else:
            rejected_count += 1
        bar.update(1)
    click.echo(red_fg(">>> Accepted %d and rejected %d blast hits for file %s ..." % (accepted_count,
                                                                                      rejected_count, blasttxt)))
    click.echo(magenta_fg(">>> Inserting junctions into "
                          "database %s ..." % os.path.basename(blast_parsed_results_filepath)))
    for key in tqdm(parsed_results.keys(), unit=" junctions"):
        gene_name, nm_number, frame, orf, inframe_inorf, position, query_start = key.split("|")
        count = parsed_results[key]
        gene = Gene.select().where(Gene.gene_name == gene_name)
        Junction.insert(gene=gene, position=position, query_start=query_start,
                        frame=frame, orf=orf, ppm=0.0, inframe_inorf=inframe_inorf, count=count).execute()

    click.echo(green_fg(">>> Generating gene stats for database %s ..." % os.path.basename(blast_parsed_results_filepath)))
    generate_stats(blast_count)
    finish = time.time()
    hr, min, sec = elapsed_time(start, finish)
    click.echo(cyan_fg("Finished parsing blast file %s in time %d hr, %d min, %d sec" % (blasttxt, hr, min, sec)))
    jdb.close_db()


def parse_blast_results(directory, blast_results_folder, blast_results_query_folder, gene_list_file, threads):
    blast_results_list = get_file_list(directory, blast_results_folder, ".txt")
    click.echo(cyan_fg('>>> Parsing blast results on %s cores.' % threads))
    parallel.Parallel(n_jobs=threads)(parallel.delayed(_parse_blast_results)(directory,
                                                                             blast_results_folder, f,
                                                                             blast_results_query_folder,
                                                                             gene_list_file) for f in blast_results_list)

