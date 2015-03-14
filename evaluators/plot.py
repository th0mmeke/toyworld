"""
Created on 6/03/2013

@author: thom
"""

from evaluator import Evaluator

import matplotlib.pyplot as plt
import logging


class Plot(Evaluator):

    def get_result_titles(self):
        return []

    def evaluate(self, results_filename, **kwargs):

        f1 = plt.figure()
        self.draw_figure(f1, results_filename, **kwargs)
        logging.info("Saving plot to {}".format(kwargs["output_filename"]))
        f1.savefig(kwargs["output_filename"], format='eps')
