#/usr/bin/env python3
from data import *

for align_sv_type in sv_types:
    dataset = align.reuse_dataset(default_ref, sv_type=align_sv_type, sv_nbr=5, length_mean=9000, length_sd=7000, depth=20, chromosome_nbr=1, force=True)
    dataset["snifflescoverage"] = 5

    print(json.dumps(dataset, indent=4))

    for subsegment in [512]:
        print(f"Running with subsegment length {subsegment}\n")
        dataset["subsegment"] = subsegment
        dataset['max_sv_distance'] = 600
        align.update_dataset(dataset)
        # align.align(align.git_dir(), dataset)
        # align.bam_convert(dataset)
        align.sv_call(dataset)
        align.eval(dataset)

