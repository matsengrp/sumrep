#!/usr/bin/env python
import csv
import sys

partis_path = "/home/bolson2/Software/partis"  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')

import utils
import glutils

# read default germline info
glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus='igh')

arg_count = len(sys.argv)

if arg_count < 3:
    output_filename = "new_output.csv"
    if arg_count == 1:
        input_filename = "partis_output.csv"
else:
    input_filename, output_filename = sys.argv[1:3]

with open(input_filename) as infile:
    with open(output_filename, 'w') as outfile:
        reader = csv.DictReader(infile)
        writer_fieldnames = utils.annotation_headers + [name for name in utils.implicit_linekeys] 
        writer = csv.DictWriter(outfile, writer_fieldnames)
        writer.writeheader()
        for line in reader:
            if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
                continue
            utils.process_input_line(line)
            utils.add_implicit_info(glfo, line)
            line_to_write = {key:value for key, value in utils.get_line_for_output(line).items() \
                    if key in writer_fieldnames} 
            writer.writerow(line_to_write)

