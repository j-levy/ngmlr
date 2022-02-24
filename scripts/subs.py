#/usr/bin/env python3

from data import *
import json, os

# recompile

os.chdir(f"{align.git_dir()}/bld")
# os.system("rm -r *")
os.system("cmake -DCMAKE_BUILD_TYPE=Debug .. && cmake --build . -j6")
os.chdir(f"{align.git_dir()}")

for align_sv_type in sv_types:
    dataset = align.generate_dataset(default_ref, sv_type=align_sv_type, sv_nbr=sv_nbr, length_mean=length_mean, length_sd=length_sd, depth=depth, chromosome_nbr=chromosome_nbr, force=force)
    dataset["snifflescoverage"] = 3

    for subsegment in subsegments:
        print(f"Running with subsegment length {subsegment}\n")
        dataset["subsegment"] = subsegment
        dataset['max_sv_distance'] = 4
        dataset["nick"] = f"{length_mean}"
        dataset["threads"] = 5
        align.update_dataset(dataset)
        print(json.dumps(dataset, indent=4))
        align.align(align.git_dir(), dataset, is_profiling, args=args)
        align.bam_convert(dataset)
        align.sv_call(dataset)
        align.eval(dataset)

