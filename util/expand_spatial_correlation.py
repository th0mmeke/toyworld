import config

import csv
import sys
import os
import json

if __name__ == "__main__":

    fname = sys.argv[1]
    filename = os.path.join(config.DataDir, fname)
    new_filename = filename + "-expanded"

    with open(filename, 'rb') as csvfile:
        with open(new_filename, 'wb') as outfile:
            spamreader = csv.reader(csvfile, delimiter='{')
            for row in spamreader:
                for x in row:
                    if x[-1] == ",":
                        x = x[:-1]
                    if x[-1] != "}":
                        header = x
                    else:
                        x = "{" + x
                        x = str.replace(x, "'", '"')
                        val = json.loads(x)
                        outfile.write("{},{},{},{}\n".format(header, val['t'], val['g'], val['n']))
