import os
import sys

class FileSystem:
    def __init__(self):
        pass

    def makedir_orchange(path): ##This function is the equivalent of the mkdir -p bash function
        """Returns the base directory, the base name (without extension) and the base directory of the base directory"""
        try:
            os.mkdir(path) ##equivalent of bash command mkdir for specified path 
            return path
        except FileExistsError:
            return path

    def get_base_dir(path):
        """Returns the base directory, the base name (without extension) and the base directory of the base directory"""
        base, ext = os.path.splitext(path) ##Split the path by /: if the path is /etc/nice7path/file.fasta you get base=/etc/nice7path/file, ext=.fasta
        asedir = base.split("/") ##Split the base by /: if the base is /etc/nice/path/file you get ["etc", "nice", "path", "file"]
        b = "/"
        basedir = b.join(asedir[:len(asedir)-1]) ##Get basedir: if the base is /etc/nice/path/file you get /etc/nice/path
        base_basedir = b.join(asedir[:len(asedir)-2]) ##Get the parent of basedir: if the base is /etc/nice/path/file you get /etc/nice/
        basename = asedir[len(asedir)-1] ##Get the basename: if the base is /etc/nice/path/file you get file
        return basedir, basename, base_basedir


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        try:
            y = gzip.open(x, "rt", encoding="latin-1")
        except FileNotFoundError:
            print(f"File {x} does not exist, proceeding with the analysis...", file=sys.stderr)
            return []
    else:
        try:
            y = open(x, "rt", encoding="latin-1")
        except FileNotFoundError:
            print(f"File {x} does not exist, proceeding with the analysis...", file=sys.stderr)
            return []
    return y
