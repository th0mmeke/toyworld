#!/usr/bin/env python
"""
Run a simulation or evaluate the output of an earlier simulation run.

-e, --experiment  Filename of Experiment design (xml)
--evaluator  Module.Class for evaluator - use instead of evaluators in experiment design (optional)
-l, --log_level  Set the logging level (optional)
-f, --log_filename  Filename for logging (relative to location of evaluation design file) (optional)
-a, --partition If true, partition the input data into four windows equally spaced in the input stream, each of 250 reactions long

Examples:

* python toyworld/main.py -e experiment-new/experiment_design.xml

* python toyworld/evaluate.py -e experiment-new/experiment_design.xml -l DEBUG -f new.log
* python toyworld/evaluate.py -e experiment-new/experiment_design.xml --evaluator evaluators.evaluator_iteration.IterationEvaluator
"""

import os
import logging
import argparse
import importlib
import string
import sys

import xml.etree.cElementTree as ElementTree

import config
import runner


def evaluate_experiments(experiments, evaluations, dirname, output):
    """Run a series of evaluations against a series of experiments.
    All evaluators are run against an experiment before moving on to the next experiment."""

    partition_length = 250  # number of iterations to include in each input partition, if partitioning requested
    output_filehandles = {}

    for evaluation in evaluations:
        if len(evaluation.get_result_titles()) > 0:

            # Configure output relative to dirname
            evaluator_class = evaluation.__class__.__name__
            output_filename = os.path.join(dirname, "{}.out".format(evaluator_class))
            logging.info("Will write output from evaluator {} to {}".format(evaluator_class, output_filename))
            fh = open(output_filename, "w")
            output_filehandles[evaluator_class] = fh
            factor_titles = experiments[0].get_factor_titles()

            header_string = "Experiment, Repeat, Partition Start, Partition End, "
            if len(factor_titles) > 0:
                fh.write(header_string + "{}, {}\n".format(string.join(factor_titles, ', '),
                                                           string.join(evaluation.get_result_titles(), ', ')))
            else:
                fh.write(header_string + "{}\n".format(string.join(evaluation.get_result_titles(), ', ')))

    for experiment in experiments:

        if experiment.end_iteration > partition_length:
            partition_generator = range(experiment.end_iteration / 4, experiment.end_iteration + 1,
                                        experiment.end_iteration / 4)
            experiment_partitions = [[{'start': i - partition_length, 'end': i}] for i in partition_generator if
                                     i > partition_length]
        else:
            experiment_partitions = [None]

        logging.info(
            "Evaluating {} with {} repeats and results from {}".format(experiment.name, experiment.repeats, dirname))

        for repeat in range(experiment.repeats):
            data_filename = experiment.get_results_filename(repeat)
            if not os.path.exists(data_filename):
                logging.warn("Data file for this repeat is missing ({})".format(data_filename))
            else:
                for evaluation in evaluations:

                    evaluator_class = evaluation.__class__.__name__
                    logging.info("Running evaluator: {}.{}".format(evaluation.__class__.__module__, evaluator_class))

                    selected_partitions = (experiment_partitions if evaluation.is_partitioned() else [None])

                    for partition in selected_partitions:
                        experiment_details = [experiment.name, str(repeat)]

                        if partition is None:
                            experiment_details.extend(["0", str(experiment.end_iteration)])
                            output_filename = os.path.join(dirname, "{}-{}-{}.eps".format(evaluator_class, experiment.name, repeat))
                        else:
                            output_filename = os.path.join(dirname, "{}-{}-{}-{}.eps".format(evaluator_class, experiment.name, repeat, partition[0]['start']))
                            experiment_details.extend([str(partition[0]['start']), str(partition[0]['end'])])

                        if len(experiment.get_factors()) > 0:
                            experiment_details.extend([{'High': '+1', 'Low': '-1'}[raw_factor] for raw_factor in experiment.get_factors()])

                        results = evaluation.evaluate(results_filename=data_filename,
                                                      states_filename=experiment.get_states_filename(repeat),
                                                      selected_partitions=partition, output_filename=output_filename,
                                                      experiment=experiment)

                        if type(results) != list:
                            results = [results]

                        for result in results:

                            if len(output_filehandles) > 0:
                                if evaluation.__class__.__name__ in output_filehandles.keys():
                                    results_string = string.join([str(x) for x in result], ',')
                                    string_output = string.join(experiment_details, ',') + "," + results_string
                                    logging.info("Writing {} to output file".format(string_output))
                                    output_filehandles[evaluation.__class__.__name__].write(string_output + "\n")
                                    # output_filehandles[evaluation.__class__.__name__].flush()  # force python file write
                                    # os.fsync(output_filehandles[evaluation.__class__.__name__].fileno())  # force OS file write

                        del results

    for f in output_filehandles.values():
        f.close()


def initialise_logging(args, basedir):
    level = getattr(logging, args.log_level.upper())
    logger = logging.getLogger()
    logger.setLevel(level)

    if not args.log_filename:
        args.log_filename = os.path.basename(basedir) + ".log"
    fh = logging.FileHandler(os.path.join(basedir, args.log_filename))
    fh.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)

    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)


def get_args():
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-e', '--experiment', required=True, help="Filename of Experiment design (xml)")
    parent_parser.add_argument('-m', '--evaluator', default=None,
                               help="Module.Class for evaluator (e.g., evaluators.evaluator_summary.SummaryEvaluator): use instead of evaluators in experiment design")
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level")
    parent_parser.add_argument('-f', '--log_filename',
                               help="Filename for logging (relative to location of evaluation design file) (optional)")

    if os.path.basename(sys.argv[0]) == "evaluate.py":
        parser.add_argument('-a', '--partition', action='store_true', help='Partition input data')
        parser.add_argument('-p', '--preview', action='store_true', help='Do not write evaluator output to file')

    parser = argparse.ArgumentParser(parents=[parent_parser])

    return parser.parse_args()


if __name__ == "__main__":

    args = get_args()
    basedir = os.path.join(config.DataDir, os.path.dirname(args.experiment))
    initialise_logging(args, basedir)

    basename = os.path.basename(args.experiment)
    if os.path.splitext(args.experiment)[1] == '':
        basedir = os.path.join(basedir, basename)
        basename = 'experiment_design.xml'

    experiment_filename = os.path.join(basedir, basename)
    if not os.path.isfile(experiment_filename):
        sys.exit("Couldn't find an experiment file at {0}".format(experiment_filename))

    logging.info("Reading experiment design from file {}".format(experiment_filename))
    runner = runner.Runner(ElementTree.parse(experiment_filename), basedir)

    if os.path.basename(sys.argv[0]) != "evaluate.py":
        runner.run()
    else:
        experiments = runner.get_experiments()

        if args.evaluator is None:
            evaluations = runner.get_evaluations()
        else:
            evaluator_module, evaluator_class = args.evaluator.rsplit(".", 1)
            evaluator_module = evaluator_module.rstrip()
            evaluator_class = evaluator_class.rstrip()
            evaluations = [
                getattr(importlib.import_module(evaluator_module), evaluator_class)(partition=args.partition)]

        evaluate_experiments(experiments, evaluations, basedir, not args.preview)

        logging.info("Finished evaluation")
