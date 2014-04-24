#!/usr/bin/env python

import argparse
import os
import sys
import os.path
from glob import glob
import subprocess

import fix_fastq

def main( args ):
    run_sample( args.readdir, args.outdir )

def run_sample( fqdir, outdir ):
    fixf_o = os.path.join( outdir, 'fix_fasta' )
    miseqfq = get_miseq_reads( fqdir )
    fixfqs = fix_fastq.fix_fastqs( fixf_o, miseqfq )
    flash_o = os.path.join( outdir, 'flash' )
    prefix = 'out'
    flasho = flash( *fixfqs, o=prefix, d=flash_o )

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
        out = subprocess.check_output( cmd )
    except subprocess.CalledProcessError as e:
        print "flash returned with return code " + e.returncode
        print "and error: " + e.output
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
        help='Location of read files. Miseq reads will only be used.'
    )

    outdir = 'output'
    parser.add_argument(
        '-o',
        '--outdir',
        default=outdir,
        help='The directory to put everything in for the sample[Default: {}]'.format(outdir)
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
