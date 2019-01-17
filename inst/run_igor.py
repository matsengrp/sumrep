import sys, os
import pandas as pd
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np
import random
import re

# May need to move this to /fh/fast/
sys.path.insert(0, '/home/bolson2/Software/IGoR/pygor')

import pygor
from pygor import counters, aligns
from pygor.counters import bestscenarios
from pygor.counters.bestscenarios import bestscenarios
from pygor.aligns import CDR3_tools


def run_igor(input_file, \
             igor_wd, \
             chain, \
             num_scenarios=1, \
             eval_batch_name="igor", \
             gen_batch_name="igor", \
             species="human" \
             ):
    df = pd.read_csv(input_file)
    df.columns = ['Sequence']
    generated_seq_count = df.shape[0]
    
    igor_command = "sh inst/run_igor.sh" + \
            " -w " + igor_wd + \
            " -i " + str(input_file) + \
            " -n " + str(num_scenarios) + \
            " -g " + str(generated_seq_count) + \
            " -e " + eval_batch_name + \
            " -b " + gen_batch_name + \
            " -c " + chain + \
            " -s " + species 
    
    os.system(igor_command)
    
# Get and write annotation dataframe via pygor subpackage routine
def get_annotations(igor_wd, 
                    chain,
                    scenarios_filename,
                    output_filename):
    scenarios_file = os.path.join(igor_wd, scenarios_filename)
    model_parms_filename = "igor_evaluate/final_parms.txt"
    model_parms_file = os.path.join(igor_wd, model_parms_filename)
    annotation_dat = bestscenarios.read_bestscenarios_values(scenarios_file, model_parms_file)

    annotation_dat['v_call'] = annotation_dat['v_choice'].apply( \
            lambda x: extract_gene_from_full_string(x, chain, "V") )
    annotation_dat['j_call'] = annotation_dat['j_choice'].apply( \
            lambda x: extract_gene_from_full_string(x, chain, "J") )

    if chain == "beta":
        annotation_dat['d_call'] = annotation_dat['d_gene']
        annotation_dat = annotation_dat.rename(index=str, \
                columns={"v_3_del": "v_3p_del", \
                         "d_5_del": "d_5p_del", \
                         "d_3_del": "d_3p_del", \
                         "j_5_del": "j_5p_del", \
                         "vd_ins": "np1_length", \
                         "dj_ins": "np2_length", \
                         "vd_dinucl": "vd_insertion", \
                         "dj_dinucl": "dj_insertion" })
    else: 
        annotation_dat = annotation_dat.rename(index=str, \
                columns={"v_3_del": "v_3p_del", \
                         "j_5_del": "j_5p_del", \
                         "vj_ins": "np1_length", \
                         "vj_dinucl": "vj_insertion" })

    # Write annotations dataset to file
    annotation_dat.to_csv(os.path.join(igor_wd, output_filename))

def extract_gene_from_full_string(s, chain, gene):
    chain_letter = {"beta": "B",
                    "alpha": "A"}
    regex = "TR" + chain_letter[chain] + gene + "[0-9]+([-][0-9]+)?[*][0-9]+"
    regex_result = re.search(regex, s)
    gene_string = regex_result.group(0) if regex_result is not None else None
    return(gene_string)

def run_igor_analysis(input_file,
                      igor_dir_name,
                      chain,
                      output_filename
                     ):
    cwd = os.getcwd()
    igor_directory = os.path.join(cwd, igor_dir_name)
    if not os.path.exists(igor_directory):
        os.makedirs(igor_directory)
    run_igor(input_file=input_file, \
            igor_wd=igor_directory, \
            chain=chain)
    obs_scenarios_filename = "igor_output/best_scenarios_counts.csv"
    get_annotations(igor_directory, chain, obs_scenarios_filename, output_filename)
    sim_scenarios_filename = "igor_generated/generated_realizations_werr.csv"
    get_annotations(igor_directory, chain, sim_scenarios_filename, "sim.csv")

if __name__ == "__main__":
    cmv_annotations = run_igor_analysis(
        input_file=sys.argv[1],
        igor_dir_name=sys.argv[2],
        chain=sys.argv[3],
        output_filename=sys.argv[4]
    )
