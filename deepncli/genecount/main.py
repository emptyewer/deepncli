import joblib.parallel as parallel
from ..utils.io import get_sam_filelist


# def lets_count(directory, input_data_folder, sam_file, summary_folder, exon_file)


def count_genes(directory, input_data_folder, summary_folder, exon_file):
    sam_files_list = get_sam_filelist(directory, input_data_folder)
    num_cores = parallel.cpu_count()
    if len(sam_files_list) > 0:
        parallel(n_jobs=num_cores - 1)(
            parallel.delayed(lets_count)(directory, input_data_folder, sam_file, summary_folder, exon_file) for sam_file
            in
            sam_files_list)
