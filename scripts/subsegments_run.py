#/usr/bin/env python3

import align
import json

default_ref = "/media/DataXFS/ngmlr_data/Genomes/c_elegans/c_elegans.fasta"

sv_types = [align.SvType.INV, align.SvType.DUP, align.SvType.INDEL, align.SvType.INVDEL]
for align_sv_type in sv_types:
    dataset = align.generate_dataset(default_ref, align_sv_type, sv_nbr=15, length_mean=9000, length_sd=7000, depth=20, force=False)

    # dataset = align.reuse_dataset(default_ref, align.SvType.INV, sv_nbr=15, length_mean=9000, length_sd=7000, depth=20)
    dataset["snifflescoverage"] = 5

    print(json.dumps(dataset, indent=4))

    for subsegment in [256, 512, 1024, 2048]:
        print(f"Running with subsegment length {subsegment}\n")
        dataset["subsegment"] = subsegment
        dataset['max_sv_distance'] = 20
        align.update_dataset(dataset)
        align.align(align.git_dir(), dataset)
        align.bam_convert(dataset)
        align.sv_call(dataset)
        align.eval(dataset)

"""
dataset = align.generate_dataset(default_ref, align.SvType.INDEL, sv_nbr=20, length_mean=9000, length_sd=7000, depth=20, force=False)
print(json.dumps(dataset, indent=4))

dataset["snifflescoverage"] = 5

for subsegment in [256, 512, 1024, 2048]:
    print(f"Running with subsegment length {subsegment}\n")
    dataset["subsegment"] = subsegment
    align.update_dataset(dataset)
    align.align(align.git_dir(), dataset)
    align.bam_convert(dataset)
    align.sv_call(dataset)
    align.eval(dataset)
"""