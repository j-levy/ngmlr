#/usr/bin/env python3

import align
import json

default_ref = "/media/DataXFS/ngmlr_data/Genomes/c_elegans/c_elegans.fasta"

dataset = align.generate_dataset(default_ref, align.SvType.INV, depth=1) #just checking the whole pipeline *runs*
print(json.dumps(dataset, indent=4))

for subsegment in [256]:
    print(f"Running with subsegment length {subsegment}\n")
    dataset["subsegment"] = subsegment
    align.align(align.git_dir(), dataset)
    align.bam_convert(dataset)
    align.sv_call(dataset)
    align.eval(dataset)