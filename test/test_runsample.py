from imports import *
import common

class Base( common.Base ):
    fixdir = join( dirname(__file__), 'fixtures', 'fix_fastq' )
    fqs = glob( join( fixdir, '*.fastq' ) )
    bnfqs = [basename(fq) for fq in fqs]
    flash_files = (
        '{}.extendedFrags.fastq',
        '{}.notCombined_1.fastq',
        '{}.notCombined_2.fastq',
        '{}.hist',
        '{}.histogram'
    )

    def flash_output_files( self, prefix, outdir ):
        files = []
        for f in self.flash_files:
            files.append( join( outdir, f ).format(prefix) )
        return files

class TestGetMiSeqReads(Base):
    def _C( self, *args, **kwargs ):
        from runsample import get_miseq_reads
        return get_miseq_reads( *args, **kwargs )

    def test_gets_only_miseq( self ):
        r = self._C( self.fixdir )
        eq_( 2, len(r) )

class TestUnitBuildOptions(Base):
    def _C( self, *args, **kwargs ):
        from runsample import build_options
        return build_options( *args, **kwargs )

    def test_handles_flag_options( self ):
        r = self._C( b=True, c=False )
        eq_( ['-b'], r )

    def test_handles_numbers( self ):
        r = self._C( i=5, f=5.5 )
        eq_( ['-i','5','-f','5.5'], r )

class TestUnitFlash(Base):
    def _C( self, *args, **kwargs ):
        from runsample import flash
        return flash( *args, **kwargs )

    def test_runs_with_options( self ):
        ''' -o outprefix -d outdir '''
        efiles = self.flash_output_files( 'test', 'flashout' )
        flashout = self._C( *self.fqs, o='test', d='flashout' )
        eq_( set([]), set(efiles) - set(flashout) )
        for f in flashout:
            ok_( exists( f ), 'Did not create {}'.format(f) )

class TestFunctional(Base):
    def _C( self, *args, **kwargs ):
        script = join( dirname(dirname(__file__)), 'runsample.py' )
        cmd = [script]
        cmd += ['-o', kwargs.get('o','output')]
        cmd += list(args)
        return subprocess.call( cmd )

    def test_fixed_fastqfiles( self ):
        r = self._C( self.fixdir )
        outdir = join( 'output', 'fix_fasta' )
        ok_( exists(outdir), 'Did not create fix_fasta output dir' )
        rfqs = os.listdir( join('output','fix_fasta') )
        eq_( set([]), set(self.bnfqs) - set(rfqs) )

    def test_flash_files( self ):
        r = self._C( self.fixdir )
        outdir = join( 'output', 'flash' )
        prefix = 'out'
        efiles = self.flash_output_files( prefix, outdir )
        for f in efiles:
            ok_( exists(f), 'Did not create flash output file {}'.format(f) )
