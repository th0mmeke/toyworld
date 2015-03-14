import config

import csv
import sys
import os
import string

"""Experiment, Repeat, Partition Start, Partition End, Number of cycles, Length of longest cycle, Count of most common cycle
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
                experiment_name = row[0]
                energy = experiment_name.split('-')[-1]
                row.insert(4, energy)
                # print(string.join(row,","))
                outfile.write("{}\n".format(string.join(row, ",")))
