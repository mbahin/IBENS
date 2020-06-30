#!/usr/bin/python

__author__ = "bahin"
""" <description> """

import argparse
import sys

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="input", required=True, help="<description>")
parser.add_argument("-f", dest="flag", action="store_true", help="<description>")
parser.add_argument("-o", dest="output", required=True, help="<description>")
options = parser.parse_args()

# Importing the "common_func" module
sys.path.insert(1, "/data/biocomp/bahin/biocomp-mb-scripts/Scripts")
from common_func import init_dict
sys.path.remove("/data/biocomp/bahin/biocomp-mb-scripts/Scripts")

# Parsing the input file
with open(options.input, "r") as input_file:
    for line in input_file:
