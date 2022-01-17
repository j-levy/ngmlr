import align

"""

    TRANS = "TRANSLOCATION" # translocation will be disabled (using 1 chromosome to make things faster)
"""

default_ref = "/media/DataXFS/ngmlr_data/Genomes/s_cerevisiae/s_cerevisiae.fasta"

sv_types = [align.SvType.DUP, align.SvType.INV, align.SvType.INDEL] #, align.SvType.INV, align.SvType.INVDEL, align.SvType.INVDUP] # , align.SvType.TRANS 
sv_nbr=8
length_mean=9000
lengths_mean = [25000, 30000, 35000, 40000, 50000, 60000]
length_sd=7000
depth=40
chromosome_nbr=1
force=False
is_profiling = False

subsegments = [256] # 128, 256, 384, 512, 768, 1024, 2048