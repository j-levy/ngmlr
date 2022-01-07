#/usr/bin/env python3

from data import *

for align_sv_type in sv_types:
    dataset = align.reuse_dataset(default_ref, sv_type=align_sv_type, sv_nbr=sv_nbr, length_mean=length_mean, length_sd=length_sd, depth=depth, chromosome_nbr=chromosome_nbr, force=force) 
    dataset["snifflescoverage"] = 5

    print(json.dumps(dataset, indent=4))

    for subsegment in subsegments:
        print(f"Running with subsegment length {subsegment}\n")
        dataset["subsegment"] = subsegment
        dataset['max_sv_distance'] = 1000
        align.update_dataset(dataset)
        # align.align(align.git_dir(), dataset)
        # align.bam_convert(dataset)
        align.sv_call(dataset)
        align.eval(dataset)

