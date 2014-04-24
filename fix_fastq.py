#!/usr/bin/env python

import argparse
import os
import os.path
import sys

def main( args ):
    fix_fastqs( args.outdir, args.fastq )

def fix_fastqs( outdir, fastq ):
    outfq = []
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )
    for fq in fastq:
        newname = os.path.join( outdir, os.path.basename(fq) )
        with open( newname, 'w' ) as fh:
            for seq, id, qual in parse_fq(fq):
                fh.write( mod_fq_read( seq, id, qual ) )
        outfq.append( newname )
    return outfq

def parse_fq( fastq ):
    with open( fastq ) as fh:
        id = ''
        seq = ''
        qual = ''
        for line in fh:
            line = line.rstrip()
            if id == '':
                id = line
                continue
            elif seq == '':
                seq = line
                continue
            elif line.startswith('+'):
                continue
            elif qual == '':
                qual = line
                yield (id,seq,qual)
                # Reset for new sequence
                id = ''
                seq = ''
                qual = ''
            else:
                raise ValueError( 'Incorrect sequnce' )

def mod_fq_read( id, seq, qual ):
    read_template = '{}\n{}\n+\n{}\n'
    newid = miseq_to_newbler_id( id )
    return read_template.format( newid, seq, qual ) 

def miseq_to_newbler_id( id ):
    l,r = id.split()
    r12 = r[0]
    return '{}#0/{} ({})'.format(l,r12,id)

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Changes the seq.id of miseq reads to contain #0/1 or #0/2 at the end for newbler'''
    )

    outdir_default = 'outdir'
    parser.add_argument(
        '-o',
        '--outdir',
        dest='outdir',
        default=outdir_default,
        help='Where to put the modified fastq files[Default: {}]'.format(outdir_default)
    )

    parser.add_argument(
        'fastq',
        nargs='+',
        help='Fastq file paths'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
