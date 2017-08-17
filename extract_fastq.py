import h5py
import sys
import os
import warnings

# Creates multi-fastq from fast5-files, fastq-files, or a mix of both.

args = sys.argv[1:]
if len(args) > 2:
    raise ValueError('Provide two arguments: input folder and output file name')
in_path = args[0]
if in_path[-1] != '/':
    in_path += '/'
in_files = os.listdir(in_path)
out_file = args[1]

# Remove file if exists
if os.path.isfile(out_file):
    os.remove(out_file)

for f in in_files:
    if f[-6:] == '.fast5':
        with h5py.File(in_path+f) as f5:
            try:
                fastq = f5['Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]
            except KeyError:
                print('Could not retrieve fasta for read %s' % f)
                continue
            fastq = fastq.decode('utf-8')
    elif f[-3:] == '.fq' or f[-6:] == '.fastq':
        with open(in_path+f, 'r') as fq:
            fastq = fq.read()
    else:
        warnings.warn('unrecognized file type for %s' % f, RuntimeWarning)
        continue
    with open(out_file, 'a+') as of:
        of.write(fastq)
