# ToyWorld

An artificial chemistry in Python. Makes use of [RDKit](https://github.com/rdkit/rdkit) for base chemical operations. RDKit provides a number of useful capabilities including format conversions to and from SMILES [4] and graphical forms of molecules; standard sanity checks for molecular structure, and molecular manipulations.

Molecules are modelled as an extension of standard RDKit Mol objects, constructed from RDKit Atoms connected with Bonds. Standard Lewis dot structures built on the inherited atomic properties are used to identify possible bonds, and a formal charge model is used to record the charge changes associated with modi􏰀cations to the molecular structure caused by reactions.

## Installation

Requires Python 2.7 with python packages python-networkx, python-numpy, python-rdkit, pygame and pymunk.

### Docker

This is the simplest way to use Toyworld.

`docker run -v=<path to data directory>:/toyworld/data th0mmeke/toyworld`

The data directory must contain `experiment_design.xml` and any referenced files (such as a population file). Results will also be returned in this directory.

### Ubuntu
Install pre-reqs:

    apt-get update
    apt-get -y install python-networkx python-numpy python-rdkit git-core
    apt-get install wget zip
    wget https://pymunk.googlecode.com/files/pymunk-4.0.0.zip -O /tmp/pymunk.zip
    cd /tmp; unzip pymunk.zip; cd pymunk-*; python setup.py install

Set environment variables:

    export TOYWORLD=<TOYWORLD_BASE>
    export TOYWORLDDATA=<full path to your folder or directory holding your experiment files> # Optional - default is $TOYWORLD/data
    export PYTHONPATH=$PYTHONPATH:$TOYWORLD

Finally, in directory ../$TOYWORLD: `git clone https://github.com/th0mmeke/toyworld.git`

### MacOS

Install RDKit through Homebrew:

    ruby -e “$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)”
    brew doctor
    brew tap rdkit/rdkit
    brew install rdkit

Install pre-reqs:

    sudo easy_install numpy
    sudo easy_install networkx
    sudo easy_install pymunk
    sudo easy_install pygame

Then for PyGame:

    brew install sdl sdl_image sdl_mixer sdl_ttf portmidi
    sudo pip install hg+http://bitbucket.org/pygame/pygame

Set environment variables:

    export RDBASE=/usr/local/share/RDKit
    export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/site-packages
    export TOYWORLD=<TOYWORLD_BASE>
    export TOYWORLDDATA=<full path to your folder or directory holding your experiment files> # Optional - default is $TOYWORLD/data
    export PYTHONPATH=$PYTHONPATH:$TOYWORLD

Finally, in directory ../$TOYWORLD: `git clone https://github.com/th0mmeke/toyworld.git`

## Usage

  1. Describe an experiment (create an experiment definition in XML) or sequence of experiments (`python toyworld/doe.py`)
  1. Run the experiment (`python toyworld/main.py -e <experiment file>`)
  1. Evaluation (`python toyworld/evaluate.py -e <experiment file>`)

## Configuration

Experiments are defined in an XML experiment definition. An example can be found in `docs/experiment_design.xml`. The experiment definition must include a reference to a population description, detailing the initial molecular population (see `docs\mixed_population.xml` for a sample.)


