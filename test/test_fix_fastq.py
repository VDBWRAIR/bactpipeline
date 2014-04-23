from imports import *
import common

class Base( common.Base ):
    pass

class TestUnitMiSeqToNewbler( Base ):
    def _C( self, *args, **kwargs ):
        from fix_fastq import miseq_to_newbler_id
        return miseq_to_newbler_id( *args, **kwargs )

    def test_r1_correct( self ):
        r = self._C( 'abcd 1' )
        eq_( 'abcd#0/1 (abcd 1)', r )

    def test_r2_correct( self ):
        r = self._C( 'abcd 2' )
        eq_( 'abcd#0/2 (abcd 2)', r )

class TestUnitModFqRead( Base ):
    def _C( self, *args, **kwargs ):
        from fix_fastq import mod_fq_read
        return mod_fq_read( *args, **kwargs )

    def test_mods_correctly( self ):
        from fix_fastq import miseq_to_newbler_id as mtni
        id = 'abcd 1'
        seq = 'ATGC'
        qual = 'IIII'
        r = self._C( id, seq, qual )
        read = '{}\n{}\n+\n{}\n'.format(mtni(id),seq,qual)
        eq_( read, r )

class TestUnitParseFq( Base ):
    def _C( self, *args, **kwargs ):
        from fix_fastq import parse_fq
        return parse_fq( *args, **kwargs )

    def fake_fq( self ):
        with open( 'fake.fq', 'w' ) as fh:
            for i in range( 1, 101 ):
                fh.write( '@abcd:{} {}\n'.format( i, (i%2)+1) )
                fh.write( 'ACGT\n' )
                fh.write( '+\n' )
                fh.write( 'IIII\n' )
        return 'fake.fq'

    def test_parses( self ):
        fq = self.fake_fq()
        r = self._C( fq )
        for id, seq, qual in r:
            ids = id.split()
            x = ids[0].split(':')
            eq_( '@abcd', x[0] )
            eq_( 'ACGT', seq )
            eq_( 'IIII', qual )

class TestFunctional( Base ):
    def sample_files( self ):
        fixdir = join( dirname(__file__), 'fixtures', 'fix_fastq' )
        return glob( join( fixdir, '*.fastq' ) )

    def _C( self, *args, **kwargs ):
        script = join( dirname( dirname( __file__ ) ), 'fix_fastq.py' )
        cmd = [script]
        if kwargs.get('outdir',False):
            cmd += ['-o', kwargs.get('outdir')]
        cmd += list(*args)
        print cmd
        return subprocess.call( cmd )

    def test_runs_correctly( self ):
        fastqs = self.sample_files()
        r = self._C( fastqs )
        eq_( 0, r )
        ok_( exists( 'outdir' ), 'did not create outdir by default' )
        fqs = os.listdir( 'outdir' )
        eq_( set([]), set([basename(fq) for fq in fastqs]) - set(fqs) )
