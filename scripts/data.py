import align
import json

"""
class SvType(Enum):
    NONE = "NONE" --> none ok right fine
    DUP = "DUPLICATION" --> not detected???
    INDEL = "INDEL" --> not detected???
    INVDEL = "INV_del" --> not detected? Also spawning N INV_DEL makes 3*N variants appear.
    INV = "INVERSION" --> detected just fine most of the time \o/
    INVDUP = "INV_dup" --> not detected? Also spawning N INV_DUP makes 2*N variants appear.

    TRANS = "TRANSLOCATION" # translocation will be disabled (using 1 chromosome to make things faster)
"""

default_ref = "/media/DataXFS/ngmlr_data/Genomes/s_cerevisiae/s_cerevisiae.fasta"

sv_types = [align.SvType.INDEL, align.SvType.INV, align.SvType.INVDEL, align.SvType.INVDUP, align.SvType.TRANS] #, align.SvType.INV, align.SvType.DUP, align.SvType.INDEL, align.SvType.INVDEL]
sv_nbr=100
length_mean=9000
length_sd=7000
depth=20
chromosome_nbr=None
force=False

subsegments = [128, 256, 384, 512, 768, 1024, 2048]