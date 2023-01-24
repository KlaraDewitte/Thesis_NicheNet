#!/usr/bin/env python
import glob
import os

if __name__ == "__main__":
    peak_files = glob.glob("narrowPeaks/*.bed")
    for file in peak_files:
        filename, file_extension = os.path.splitext(file)
        f_annStats = filename + '.annStats'
        f_annotated = filename + '.annotated'

        CMD_line = "annotatePeaks.pl " + file + " hg38 " + "-annStats " + f_annStats + " > " + f_annotated
        print(CMD_line)
        os.system(CMD_line)
