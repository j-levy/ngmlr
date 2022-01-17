#/usr/bin/env python3

from data import *
import json

for align_sv_type in sv_types:
    for length in lengths_mean:
        dataset = align.reuse_dataset(default_ref, sv_type=align_sv_type, sv_nbr=sv_nbr, length_mean=length, length_sd=length_sd, depth=depth, chromosome_nbr=chromosome_nbr, force=force)
        dataset["snifflescoverage"] = 5

        for subsegment in subsegments:
            print(f"Running with subsegment length {subsegment}\n")
            dataset["subsegment"] = subsegment
            dataset['max_sv_distance'] = 100
            dataset["nick"] = f"{length}"
            align.update_dataset(dataset)
            # print(json.dumps(dataset, indent=4))
            # align.align(align.git_dir(), dataset, is_profiling=True)
            # align.bam_convert(dataset)
            align.sv_call(dataset)
            align.eval(dataset)

