from __future__ import print_function

import argparse
import os
import sys
import os.path
from glob import glob
import subprocess
import re
import itertools
from pkg_resources import resource_filename
from bactpipeline import fix_fastq
from Bio import SeqIO
import csv
from functools import partial
from operator import itemgetter as get
from contracts import contract, new_contract
from functools import reduce
if sys.version_info[0] == 2:
    from itertools import ifilter as filter, imap as map
else:
    # Python3 normalization
    xrange = range


compose2 = lambda f, g: lambda x: f(g(x))
compose  = lambda *f: reduce(compose2, f)
complement = lambda f: lambda x: not f(x)
def mapcat(func, seq):
    return itertools.chain.from_iterable(map(func, seq))
#mapcat = compose(itertools.chain.from_iterable, map)
new_contract('readable', os.path.isfile)
new_contract('exists', os.path.exists)
new_contract('directory', os.path.isdir)
SHEET_COLUMNS = ['sample_directory', 'sample_id', 'primer_file']
SUMMARY_FILE = 'summary.tsv'
#SUMMARY_FIELDS = ['sample_id', 'length', 'contig_num', 'numreads', '%total_reads', 'total_reads', 'N50']
SUMMARY_FIELDS = ['sample_id', 'length', 'contig_num', 'numreads', '%total_reads', 'N50']
SUMMARY_DELIM = '\t'

def run_command(cmd, stderr=None, print_command=True):
    if print_command:
        print(' '.join(cmd))
    p = subprocess.Popen(cmd, stderr=stderr, stdout=subprocess.PIPE)
    sout, serr = p.communicate()
    if sout and isinstance(sout, bytes):
        sout = sout.decode('utf-8')
    if serr and isinstance(serr, bytes):
        serr = serr.decode('utf-8')
    if p.returncode != 0:
        e = subprocess.CalledProcessError(p.returncode, cmd)
        e.output = sout
        raise e
    return sout

@contract
def run_sample_sheet(f, outdir, truseq):
    ''' run the pipeline on each sample within a csv/tsv file.
    sample sheet must have fields sample_directory,sample_id,primer_file
    :type f:      str,readable
    :type outdir: str,!exists
    :type truseq: str,readable
    :rtype:       int
    '''
    os.mkdir(outdir)
    def _run_sample(_dict):
        indir, sample_id, primer_file = get('sample_directory', 'sample_id', 'primer_file')(_dict)
        return run_sample(indir, os.path.join(outdir, sample_id), truseq, sample_id, primer_file)
    sheet = csv.DictReader(open(f))
    results = map(_run_sample, sheet)
    summary_data = itertools.chain.from_iterable(results)
    write_summary(summary_data, os.path.join(outdir, 'summary.tsv'))
    return 0

def main():
    args = parse_args()
    if args.sample_sheet:
        run_sample_sheet(args.sample_sheet, args.outdir, args.truseq)
    else:
        data = run_sample( args.readdir, args.outdir, args.truseq )
        write_summary(data, os.path.join(args.outdir, 'summary.tsv'))

def run_sample( fqdir, outdir, truseq, sample_id=None, primer_file=None ):
    fixf_o = os.path.join( outdir, 'fix_fasta' )
    miseqfq = get_miseq_reads( fqdir )
    fixfqs = fix_fastq.fix_fastqs( fixf_o, miseqfq )
    flash_o = os.path.join( outdir, 'flash' )
    prefix = 'out'
    flasho = flash( *fixfqs, o=prefix, d=flash_o )
    print("flash output files: {0}".format(flasho))

    btrim_o = os.path.join( outdir, 'btrim' )
    bfiles = btrim_files( [f for f in flasho if f.endswith('.fastq')], btrim_o, p=truseq, b=300, P=True, Q=True, S=True, l=100 )
    print("btrim output files: {0}".format(bfiles))
    projdir = os.path.join( outdir, 'newbler_assembly' )
    total_reads = read_count(bfiles)
    run_assembly( bfiles, o=projdir, primer_file=primer_file )
    sample_id = os.path.basename(os.path.normpath(fqdir)) if (not sample_id) else sample_id
    newbler_dir = os.path.join(projdir, 'assembly')
    contig_file = os.path.join(newbler_dir, '454AllContigs.fna')
    summary_data = make_summary(contig_file, total_reads, sample_id)
    write_summary(summary_data, os.path.join(outdir, SUMMARY_FILE))
    write_top_contigs(contig_file, os.path.join(outdir, 'top_contigs.fasta'), sample_id)
    return summary_data

@contract
def write_summary(data, outfile, delim='\t'): # data is 2d list
    ''' write data as a csv/tsv file to outfile.
    :type data:    Iterable
    :type outfile: str,!exists
    :type delim:   str '''
    header = delim.join(SUMMARY_FIELDS)
    with open(outfile, 'w') as out:
        out.write(header + '\n')
        csv.writer(out, delimiter=delim).writerows(data)

@contract
def write_top_contigs(contig_file, outfile, sample_id, top=100):
    '''write the top `top` contigs (by length) to outfile.
    :type contig_file: str,readable
    :type outfile:     str,!exists
    :type sample_id:   str
    :type top:         int,>0
    '''
    top_contigs = sorted(read_fasta(contig_file), key=seqlen)[:top]
    for c in top_contigs:
        c.id = sample_id + '_' + c.id
    SeqIO.write(top_contigs, outfile, 'fasta')

read_fasta = partial(SeqIO.parse, format='fasta')
read_fastq = partial(SeqIO.parse, format='fastq')

seqlen = lambda x: len(x.seq)
@contract
def make_summary(contig_file, total_reads, sample_id, top=100):
    ''' Get the sample id, length, id #, number of contributing reads, percent contributing over total, of each contig
    within the "AllContigs.fna" file within newbler_dir.
    :type contig_file: str,readable
    :type total_reads: int,>0
    :type sample_id:   str
    :type top:         int,>0
    :rtype             list(tuple)
    '''
    recs = list(read_fasta(contig_file))
    #recs = itertools.chain.from_iterable(imap(read_fasta, contig_files))
    lengths = list(map(seqlen, recs))
    n50 = N50(lengths)
    def get_stats(rec):
        contig_num = int(rec.id.split('contig')[-1])
        length, numreads = map(int, re.compile(r'[^\=^\s]+=([0-9]+)').findall(rec.description))
        return sample_id, length, contig_num, numreads, numreads/float(total_reads)*100, n50
    values = map(get_stats, recs)
    return values

@contract
def N_stat(lengths, N):
    '''
    maximum positive integer L such that the total number of nucleotides
    of all contigs having length >= L is at least N% of the sum of contig lengths.
    :type lengths: list[N],N>0
    :type N: float,>0
    :rtype int'''
    def is_candidate(L):
        return (sum(filter(lambda x: x >= L, lengths)) / float(sum(lengths))) >= N
    candidates = filter(is_candidate, xrange(0, sum(lengths)))
    return max(candidates)
N50 = partial(N_stat, N=0.5)

class InvalidFastqException(Exception): pass
def read_count(fq_list):
    ''' count number of all reads in fastq list '''
    _sum = 0
    for fq in fq_list:
        fq_sum = 0
        for line in open(fq):
            fq_sum += 1
        count, remainder = divmod(fq_sum, 4)
        if remainder != 0:
            raise InvalidFastqException
        _sum += count
    return _sum

def run_assembly( fastqs, **options ):
    projdir = new_assembly( options.get('o',None) )
    replace_newbler_settings( projdir, fastqs )
    # In case o was not in options we set it again
    cmd = ['runProject', projdir]
    if options.get('primer_file'):
        cmd += ['-vt', options.get('primer_file')]
    out = run_command(cmd, subprocess.STDOUT)
    if 'Usage:' in out:
        print(out)
        return 1
    else:
        return 0

def new_assembly( projdir=None ):
    cmd = ['newAssembly']
    # If project was specified then
    if projdir is not None:
        cmd += [projdir]

    out = run_command(cmd)
    p = 'Created assembly project directory (.*)'
    projpath = re.match( p, out ).group(1)
    return projpath

def newbler_fastq_read_files( fastqs ):
    filet = '''
    <File>
      <Path>{0}</Path>
      <Type>Fastq</Type>
      <FastqScoreType>Standard</FastqScoreType>
      <FastqAccnoType>General</FastqAccnoType>
      <IsPairedEnd>false</IsPairedEnd>
      <Rescore>false</Rescore>
      <removeFlag>false</removeFlag>
    </File>'''

    xml = '<ReadFiles>\n'
    xml += '\t<SetReadCount>{0}</SetReadCount>'.format(len(fastqs))
    for f in fastqs:
        xml += filet.format(f)
    xml += '\n\t</ReadFiles>'
    return xml

def replace_newbler_settings( projpath, fastqs ):
    xmlfile = os.path.join( projpath, 'assembly/454AssemblyProject.xml' )
    fh = open( xmlfile )
    xml = fh.read()
    fh.close()
    readxml = newbler_fastq_read_files( fastqs )
    rec = re.compile('<ReadFiles>.*</ReadFiles>', re.DOTALL)
    newxml = re.sub( rec, readxml, xml, 1 )
    with open( xmlfile, 'w' ) as fh:
        fh.write( newxml )

def build_options( **options ):
    '''
        Build options list from options
    '''
    cmd = []
    for o,v in options.items():
        if isinstance(v, bool):
            if v:
                cmd.append( '-'+o )
        else:
            cmd += ['-'+o,str(v)]
    return cmd

def btrim_files( fqlist, outdir, **btrimops ):
    ofiles = []
    if not os.path.exists(outdir):
        os.makedirs( outdir )

    for f in fqlist:
        # Have to set the input file in the options
        of = os.path.join( outdir, os.path.basename(f) ).replace( '.fastq', '.btrim.fastq' )
        btrimops['t'] = f
        btrimops['o'] = of
        btrim( **btrimops )
        ofiles.append( of )
    return ofiles

def btrim( **options ):
    cmd = ['btrim'] + build_options( **options )
    # Ensure outdir exists for btrim
    #outpath = options.get('o')
    ## only if outpath is specified and has at least one directory in it
    #if outpath and os.sep in outpath:
        #outdir = os.path.dirname(outpath)
        #if not os.path.isdir(outdir):
            #os.makedirs(outdir)
    try:
        out = run_command(cmd, subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print("Command run: {0}".format(' '.join(cmd) ))
        print("btrim returned with return code ".format( str(e.returncode) ))
        print("and error: " + e.output)
        return e.returncode
    return 0

def flash( mate1, mate2, **options ):
    '''
        Runs flash pre-assembly and returns the location of the output files
        Assumes flash is in the user's PATH

        @param options - Options for flash

        - out.extendedFrags.fastq      The merged reads.
        - out.notCombined_1.fastq      Read 1 of mate pairs that were not merged.
        - out.notCombined_2.fastq      Read 2 of mate pairs that were not merged.
        - out.hist                     Numeric histogram of merged read lengths.
        - out.histogram                Visual histogram of merged read lengths.
    '''
    cmd = ['flash'] + build_options(**options) + [mate1,mate2]
    try:
        out = run_command(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print("flash returned with return code {0}".format( str(e.returncode) ))
        print("and error: " + e.output)
        return []
    prefix = options.get('o','out')
    outdir = options.get('d','')
    files = (
        prefix+'.extendedFrags.fastq',
        prefix+'.notCombined_1.fastq',
        prefix+'.notCombined_2.fastq',
        prefix+'.hist',
        prefix+'.histogram'
    )
    return [os.path.join(outdir,f) for f in files]

def get_miseq_reads( readsdir ):
    mreads = glob( os.path.join( readsdir, '*_S*_L*_R*_*.fastq' ) )
    return mreads

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Runs a sample through all stages'
    )

    parser.add_argument(
        'readdir',
        nargs='?',
        help='Location of read files. Miseq reads will only be used.'
    )

    parser.add_argument(
        '-o',
        '--outdir',
        default='output',
        help='The directory to put everything in for the sample[Default: %(default)s]'
    )

    parser.add_argument(
        '-s',
        '--sample-sheet',
        default=None,
        help='samplesheet here'
    )

    parser.add_argument(
        '-t',
        '--truseq',
        default=resource_filename(__name__, 'truseq.txt'),
        help='Path to truseq.txt[Default: %(default)s]'
    )

    return parser.parse_args( args )

import luigi
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask, OptionallyCluterTask, TORQUE
from os.path import join

class ExternalFileTask(luigi.ExternalTask):
    file = luigi.Parameter()
    def output(self): return luigi.LocalTarget(self.file)

class RunSampleTask(OptionallyCluterTask):
    outdir = luigi.Parameter() # not exists
    readdir = luigi.Parameter() # directory
    sample_id = luigi.Parameter()
    primer_file = luigi.Parameter(default=None)
    local = luigi.BoolParameter()
    base_path = luigi.Parameter()
    truseq = luigi.Parameter() # isfile
    software = luigi.Parameter(default=TORQUE)

    def requires(self):
        return ExternalFileTask(file=self.readdir)

    def output(self):
        # no way to pass summary_data in-memory
        return luigi.LocalTarget(join(self.outdir, SUMMARY_FILE))

    def work(self):
        sys.path.append(self.base_path)
        run_sample(self.readdir, self.outdir, self.truseq, self.sample_id, self.primer_file)
        # get stale file handle error in run_sample

parse_sheet = compose(partial(map, get(*SHEET_COLUMNS)), csv.DictReader, open)
class RunSampleSheet(LocalSGEJobTask):
    sample_sheet = luigi.Parameter() # readable
    outdir = luigi.Parameter() # not exists
    local = luigi.BoolParameter()

    def requires(self):
        sample_dirs, ids, primer_files = zip(*(parse_sheet(self.sample_sheet)))
        outdirs = map( partial(join, self.outdir), ids)
        sample_dirs, outdirs = map(os.path.abspath, sample_dirs), map(os.path.abspath, outdirs)
        truseq = os.path.abspath(resource_filename(__name__, 'truseq.txt'))
        return [RunSampleTask(base_path=os.getcwd(), truseq=truseq, local=self.local, outdir=o, readdir=rf, primer_file=pf, sample_id=id) \
                for (o, rf, pf, id) in zip(outdirs, sample_dirs, primer_files, ids)]

    def output(self):
        def up(p): return os.path.dirname(os.path.normpath(p))
        updir = up(up(self.outdir))
        if updir: os.listdir(updir)
        return luigi.LocalTarget(join(self.outdir, 'full_summary.tsv'))

    def work(self):
        data = mapcat(lambda x: csv.DictReader(x.open('r'),delimiter=SUMMARY_DELIM), self.input())
        with self.output().open('w') as out:
            w = csv.writer(out)
            out.write(SUMMARY_DELIM.join(SUMMARY_FIELDS)+'\n')
            for d in data: w.writerow([d[k] for k in SUMMARY_FIELDS]) #
            #dw = csv.DictWriter(out, SUMMARY_FIELDS) #fails in py26
            #dw.writeheader()
            #dw.writerows(data)

if __name__ == '__main__':
    luigi.run()


