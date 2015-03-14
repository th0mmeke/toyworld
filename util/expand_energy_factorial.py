import config

import csv
import sys
import os
import string

"""Experiment, Repeat, Partition Start, Partition End, Dimensionality, Energy, Number of cycles, Length of longest cycle, Count of most common cycle
to:
Experiment, Repeat, Partition Start, Partition End, Energy, Number of cycles, Length of longest cycle, Count of most common cycle"""

if __name__ == "__main__":

    fname = sys.argv[1]
    filename = os.path.join(config.DataDir, fname)
    new_filename = filename + "-expanded"

    with open(filename, 'rb') as csvfile:
        with open(new_filename, 'wb') as outfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if row[5] == "+1":
                    row[5] = "300"
                else:
                    row[5] = "100"
                del(row[4])

                outfile.write("{}\n".format(string.join(row, ",")))
