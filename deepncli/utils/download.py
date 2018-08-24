import os
import imp
import stat
import requests
import zipfile
import platform


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
        print(">>> Downloading %s" % os.path.basename(url))
        with open(file_path, 'wb') as f:
            for data in r.iter_content(block_size):
                wrote = wrote + len(data)
                f.write(data)
        zip_ref = zipfile.ZipFile(file_path, 'r')
        zip_ref.extractall(output_dir)
        zip_ref.close()
        os.remove(file_path)
        if total_size != 0 and wrote != total_size:
            print(">>> ERROR, something went wrong")


def download_data():
    print(">>> Attempting to download data...")
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
