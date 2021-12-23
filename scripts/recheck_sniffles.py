#/usr/bin/env python3

import align
import json

default_ref = "/media/DataXFS/ngmlr_data/Genomes/c_elegans/c_elegans.fasta"

# dataset = align.generate_dataset(default_ref, align.SvType.INV, sv_nbr=20, length_mean=9000, length_sd=7000, depth=20, force=False)

dataset = align.reuse_dataset(default_ref, align.SvType.INV, sv_nbr=20, length_mean=9000, length_sd=7000, depth=20)

print(json.dumps(dataset, indent=4))

for subsegment in [512]:
    print(f"Running with subsegment length {subsegment}\n")
    dataset["subsegment"] = subsegment
    dataset["snifflescoverage"] = 3
    align.update_dataset(dataset)
    # align.align(align.git_dir(), dataset)
    # align.bam_convert(dataset)
    align.sv_call(dataset)
    align.eval(dataset)