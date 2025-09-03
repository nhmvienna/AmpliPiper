from pandas import read_csv
import os
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to the file with sample names", required=True)
argparse.add_argument(
    "-d", "--directory", help="Path to the directory where fastq files are stored", required=True)

args = argparse.parse_args()

inf = args.infile
fastqdir = args.directory


if __name__=="__main__":
    exts = ["fastq","fq",]
    names = list(read_csv(inf)["SAMPLENAME"])
    dirfiles = [os.path.join(fastqdir, f) for f in os.listdir(fastqdir) if (os.path.isfile(os.path.join(fastqdir, f)) and f.split(".")[1] in exts and f.split(".")[0] in names)]
    samples = {}
    for name in names:
        for dirfile in dirfiles:
            if os.path.basename(dirfile).split(".")[0] == name:
                samples.update({name: dirfile})
                break
            else:
                continue
    with open(os.path.join(fastqdir, "samples.csv"),"w") as c:
        for i in list(samples.keys()):
            c.write(f"{i},{samples[i]}\n")
    c.close()

