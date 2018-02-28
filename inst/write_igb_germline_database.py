#!/usr/bin/env python

# Script to process the extras.csv files from partis' germline directories
# Assumes the fasta files (e.g. ighv.fa) have been deduplicated

import sys

partis_dir = "/home/bolson2/Software/partis"
sys.path.insert(1, partis_dir + '/python')

import glutils

igb_path = "/home/bolson2/Software/igblast"
igb_database_path = igb_path + "/bin_deduplicated"

glfo = glutils.read_glfo(igb_database_path, locus='igh', debug=True)
glutils.write_glfo(igb_path + '/partis_friendly_bin', glfo, debug=True)

glfo = glutils.read_glfo(igb_database_path, locus='igk', debug=True)
glutils.write_glfo(igb_path + '/partis_friendly_bin', glfo, debug=True)

glfo = glutils.read_glfo(igb_database_path, locus='igl', debug=True)
glutils.write_glfo(igb_path + '/partis_friendly_bin', glfo, debug=True)
