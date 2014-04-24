#!/usr/bin/env python

import argparse
import os
import sys
import os.path
from glob import glob
import subprocess
import re

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
    truseq = os.path.join( os.path.dirname(__file__), 'truseq.txt' )
    btrim_o = os.path.join( outdir, 'btrim' )
    bfiles = btrim_files( [f for f in flasho if f.endswith('.fastq')], btrim_o, p=truseq, b=300, P=True, Q=True, S=True, l=100 )
    projdir = os.path.join( outdir, 'newbler_assembly' )
    run_assembly( bfiles, o=projdir )

def run_assembly( fastqs, **options ):
    projdir = new_assembly( options.get('o',None) )
    replace_newbler_settings( projdir, fastqs )
    # In case o was not in options we set it again
    cmd = ['runProject', projdir]
    out = subprocess.check_output( cmd, stderr=subprocess.STDOUT )
    if 'Usage:' in out:
        print ' '.join(cmd)
        print out
        return 1
    else:   
        return 0

def new_assembly( projdir=None ):
    cmd = ['newAssembly']
    # If project was specified then
    if projdir is not None:
        cmd += [projdir]
    
    out = subprocess.check_output( cmd )
    p = 'Created assembly project directory (.*)'
    return re.match( p, out ).group(1)

def newbler_fastq_read_files( fastqs ):
    filet = '''
    <File>
      <Path>{}</Path>
      <Type>Fastq</Type>
      <FastqScoreType>Standard</FastqScoreType>
      <FastqAccnoType>General</FastqAccnoType>
      <IsPairedEnd>false</IsPairedEnd>
      <Rescore>false</Rescore>
      <removeFlag>false</removeFlag>
    </File>'''

    xml = '<ReadFiles>\n'
    xml += '\t<SetReadCount>{}</SetReadCount>'.format(len(fastqs))
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
    newxml = re.sub( '<ReadFiles>\n.*\n.*</ReadFiles>', readxml, xml, 1, re.DOTALL )
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
    print ' '.join(cmd)
    try:
        out = subprocess.check_output( cmd, stderr=subprocess.STDOUT )
    except subprocess.CalledProcessError as e:
        print "Command run: {}".format(' '.join(cmd) )
        print "btrim returned with return code ".format( str(e.returncode) )
        print "and error: " + e.output
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
        out = subprocess.check_output( cmd, stderr=subprocess.STDOUT )
    except subprocess.CalledProcessError as e:
        print "Command run: {}".format(' '.join(cmd) )
        print "flash returned with return code {}".format( str(e.returncode) )
        print "and error: " + e.output
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
