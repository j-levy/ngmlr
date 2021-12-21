#/usr/bin/env python3

import os
import re # regular expressions
import time
import subprocess
import math

utf8 =  'utf-8'



from enum import Enum
class SvType(Enum):
    NONE = "NONE"
    DUP = "DUPLICATION"
    INDEL = "INDEL"
    INVDEL = "INV_del"
    INV = "INVERSION"
    INVDUP = "INV_dup"
    TRANS = "TRANSLOCATION" # translocation will be disabled at first

def generate_sv(ref, sv_type, sv_nbr=10, conda_bin="$HOME/miniconda3/envs/ngmlr/bin"):
    survivor_template = f"""
PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!
DUPLICATION_minimum_length: 100
DUPLICATION_maximum_length: 10000
DUPLICATION_number: 0
INDEL_minimum_length: 20
INDEL_maximum_length: 500
INDEL_number: 0
TRANSLOCATION_minimum_length: 0
TRANSLOCATION_maximum_length: 0
TRANSLOCATION_number: 0
INVERSION_minimum_length: 600
INVERSION_maximum_length: 800
INVERSION_number: 0
INV_del_minimum_length: 600
INV_del_maximum_length: 800
INV_del_number: 0
INV_dup_minimum_length: 600
INV_dup_maximum_length: 800
INV_dup_number: 0
"""

    survivor_template = survivor_template.replace(f"{sv_type.value}_number: 0", f"{sv_type.value}_number: {sv_nbr}")
    # print(f"{sv_type.value}_number: 0")
    # print(survivor_template)

    ref_dir = os.path.abspath(os.path.dirname(ref))
    ref_file = os.path.basename(ref)

    output_prefix = f"{ref_file}_{sv_type.value}_{sv_nbr}"
    survivor_param = f"{ref_dir}/{output_prefix}.svr"
    f = open(f"{survivor_param}", "w")
    f.write(survivor_template)
    f.close()
    snp_rate = 0

    # SURVIVOR simSV PATH/TO/REFERENCE.fasta PATH/TO/param.svr FREQUENCY_OF_SNP(0..1) 0(always using simulated reads) OUTPUT_PREFIX
    print(f"""cd {ref_dir} && {conda_bin}/SURVIVOR simSV {ref} {survivor_param} {snp_rate} 0 {output_prefix}""")
    os.system(f"""cd {ref_dir} && {conda_bin}/SURVIVOR simSV {ref} {survivor_param} {snp_rate} 0 {output_prefix}""")
    
    return output_prefix


def generate_dataset(ref, sv_type, sv_nbr=10, conda_bin="$HOME/miniconda3/envs/ngmlr/bin",
                        error_model="P4C2", 
                        length_mean=9000, length_sd=7000, accuracy_mean=0.85, depth=20, 
                        chromosome_nbr=4):
    diff_ratios = {
        "pacbio": "6:50:54",
        "ont": "23:31:46"
    }

    diff_ratio_name = "pacbio"
    diff_ratio = diff_ratios[diff_ratio_name]
    error_model_file = f"{git_dir()}/scripts/pbsim2_hmm_model/{error_model}.model"
    prefix = "sd"
    ref_dir = os.path.abspath(os.path.dirname(ref))
    sim_dir = f"{ref_dir}/pbsim_{depth}/"

    # if the directory already exists, skip read simulation.
    if not(os.path.exists(sim_dir) and os.path.isdir(sim_dir)):
        os.makedirs(sim_dir)
        os.system(f""" cd {sim_dir} && {conda_bin}/pbsim \
                        --hmm_model {error_model_file} \
                        --depth {depth} \
                        --difference-ratio {diff_ratio} \
                        --length-mean {length_mean} \
                        --length-sd {length_sd} \
                        --accuracy-mean {accuracy_mean} \
                        --prefix {prefix} \
                        {ref}
                """)
    else:
        print(f"Folder {sim_dir} already exists, skipping reads simulation (delete folder to regenerate)")

    dataset = dict()
    # map args to the dict. Easier to use
    dataset['read'] = f"{sim_dir}/sd_{chromosome_nbr:04}.fastq"
    dataset['ref'] =  f"{sim_dir}/sd_{chromosome_nbr:04}.ref"
    dataset['simulatedepth'] = depth
    dataset['threads'] = os.cpu_count()
    dataset['snifflescoverage'] = math.ceil(depth/4)
    dataset['subsegment'] = 256
    dataset['optimize'] = "pacbio"
    dataset['nick'] = ""
    dataset['min_sv_size'] = 100
    dataset['conda_bin'] = "$HOME/miniconda3/envs/ngmlr/bin"

    ## then generate the SVs you want
    if sv_type == SvType.NONE:
        output_prefix = f"sd_{chromosome_nbr:04}.ref"
    else:
        output_prefix = generate_sv(dataset['ref'], sv_type, sv_nbr)
        dataset['ref'] = f"{sim_dir}/{output_prefix}.fasta"


    # add some pre-computed fields to the dict
    dataset['nickname'] = f"{dataset['nick']}-{os.path.basename(dataset['read'])}-{os.path.basename(dataset['ref'])}"
    dataset['gprof_filename'] = f"gprof.{dataset['nickname']}-{dataset['subsegment']}.txt"
    dataset['output_filename'] = f"{dataset['nickname']}-{dataset['subsegment']}.sam"
    dataset['output_dir']= f"out/{dataset['nickname']}/threads-{dataset['threads']}/"
    
    return dataset


def create_sv(dataset):
    print("not implemented yet")


def align(git_dir, dataset):
    if not(os.path.exists(dataset['output_dir']) and os.path.isdir(dataset['output_dir'])):
        print(f"make folder: {dataset['output_dir']}")
        os.makedirs(dataset['output_dir'])

    start_time = time.time()
    os.system(f"""{git_dir}/bin/ngmlr-0.2.8/ngmlr \
        --bam-fix -x {dataset['optimize']} -t {dataset['threads']} \
        --subread-length {dataset['subsegment']} \
        -q {dataset['read']} -r {dataset['ref']} \
        -o {dataset['output_dir']}/{dataset['output_filename']}""")
    end_time= time.time()
    runtime = end_time - start_time
    print(f"+++++ Clock time: {runtime}")
    os.system(f"echo '+++++ Clock time: $runtime' > {dataset['output_dir']}/{dataset['gprof_filename']}")
    os.system(f"gprof {git_dir}/bin/ngmlr-0.2.8/ngmlr >> {dataset['output_dir']}/{dataset['gprof_filename']}")

def bam_convert(dataset):
    conda_bin = dataset['conda_bin']
    output_file_path = f"{dataset['output_dir']}/{dataset['output_filename']}"
    os.system(f"{conda_bin}/samtools view -o {output_file_path}.bam {output_file_path}")
    os.system(f"{conda_bin}/samtools sort -o  {output_file_path}.sort.bam  {output_file_path}.bam")
    os.system(f"{conda_bin}/samtools index  {output_file_path}.sort.bam")

def sv_call(dataset):
    conda_bin = dataset['conda_bin']
    os.system(f"{conda_bin}/sniffles -s {dataset['snifflescoverage']} -m {dataset['output_dir']}/{dataset['output_filename']}.sort.bam -v {dataset['output_dir']}/{dataset['output_filename']}.vcf")

def eval(dataset):
    conda_bin = dataset['conda_bin']
    # TODO: debug why this doesn't run
    res = subprocess.run(f"{conda_bin}/SURVIVOR", 
        "eval",  
        f"{dataset['output_dir']}/{dataset['output_filename']}.vcf",
        f"{dataset['ref'].replace('.fasta', '.bed')}",
        f"{dataset['min_sv_size']}",
        f"{dataset['output_dir']}/eval_{dataset['output_filename']}", stdout=subprocess.PIPE).stdout.decode(utf8)
    print(res)

def bash_run(command):
    return subprocess.run(command.split(" "), stdout=subprocess.PIPE).stdout.decode(utf8)

def git_dir():
    return bash_run("git rev-parse --show-toplevel").replace("\n", "")

# just remapping output_dir
# dataset['output_dir']= f"{git_dir}/docs/{args.nick}/threads-{args.threads}/"


print("align boilerplate code - version alpha 0.1.0")


# prepare_folder(dataset)
# align(git_dir, dataset)
# bam_convert(conda_bin, dataset)
# sv_call(conda_bin, dataset)
# eval(conda_bin, dataset)