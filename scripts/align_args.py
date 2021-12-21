import argparse


parser = argparse.ArgumentParser(description='Generate SV, generate reads, align them, call SV and compare the result.')


parser.add_argument("--align")
parser.add_argument("--bam_convert")
parser.add_argument("--sv_call")
parser.add_argument("--eval")

parser.add_argument("--read", dest='read', type=str)
parser.add_argument('--ref', dest='ref', type=str)

parser.add_argument('--simulatedepth', dest='depth', type=int, default=20)
parser.add_argument('--threads', dest='threads', type=int, default=os.cpu_count())
parser.add_argument('--snifflescoverage', dest='coverage', type=int, default=8)
parser.add_argument('--subsegment', dest='subsegment', type=int, default=256)
parser.add_argument('--optimize', dest='optimize', type=str, default='pacbio')
parser.add_argument('--nick', dest='nick', type=str, default='')
parser.add_argument('--min_sv_size', dest='min_sv_size', type=int, default=100)

args = parser.parse_args()

dataset = dict()

# map args to the dict. Easier to use
dataset['read'] = args.read
dataset['ref'] = args.ref
dataset['simulatedepth'] = args.depth
dataset['threads'] = args.threads
dataset['snifflescoverage'] = args.coverage
dataset['subsegment'] = args.subsegment
dataset['optimize'] = args.optimize
dataset['nick'] = args.nick
dataset['min_sv_size'] = args.min_sv_size

# add some pre-computed fields to the dict
dataset['nickname'] = f"{dataset['nick']}-{os.path.basename(dataset['read'])}-{os.path.basename(dataset['ref'])}"
dataset['gprof_filename'] = f"gprof.{dataset['nickname']}-{dataset['subsegment']}.txt"
dataset['output_filename'] = f"{dataset['nickname']}-{dataset['subsegment']}.sam"
dataset['output_dir']= f"out/{dataset['nickname']}/threads-{dataset['threads']}/"