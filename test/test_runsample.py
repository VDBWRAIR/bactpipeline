from __future__ import print_function
from imports import *
import common

fixdir = join( dirname(__file__), 'fixtures', 'fix_fastq' )
fqs = glob( join( fixdir, '*.fastq' ) )

class Base( common.Base ):
    fixdir = fixdir
    fqs = fqs
    bnfqs = [basename(fq) for fq in fqs]
    truseq = join(dirname(TEST_DIR),'bactpipeline','truseq.txt')
    assemprojxml = join(fixdir,'454AssemblyProject.xml')
    flash_files = (
        '{0}.extendedFrags.fastq',
        '{0}.notCombined_1.fastq',
        '{0}.notCombined_2.fastq',
        '{0}.hist',
        '{0}.histogram'
    )

    def flash_output_files( self, prefix, outdir ):
        files = []
        for f in self.flash_files:
            files.append( join( outdir, f ).format(prefix) )
        return files

    def flash( self, prefix='out', outdir='flash' ):
        from bactpipeline.runsample import flash
        # Files from flash
        return [f for f in flash( *self.fqs, d=outdir, o=prefix ) if f.endswith('.fastq')]

    def fqs_in_xml( self, fqlist, xml ):
        for fqfile in fqlist:
            ok_( fqfile in xml, '{0} did not make it into 454AssemblyProject.xml'.format(fqfile) )
        ok_( '<SetReadCount>{0}</SetReadCount>'.format(len(self.fqs)), 'Did not set read count correctly' )

class TestGetMiSeqReads(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import get_miseq_reads
        return get_miseq_reads( *args, **kwargs )

    def test_gets_only_miseq( self ):
        r = self._C( self.fixdir )
        eq_( 2, len(r) )

class TestUnitBuildOptions(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import build_options
        return build_options( *args, **kwargs )

    def test_handles_flag_options( self ):
        r = self._C( b=True, c=False )
        eq_( ['-b'], r )

    def test_handles_numbers( self ):
        r = self._C( i=5, f=5.5 )
        eq_( sorted(['-i','5','-f','5.5']), sorted(r) )

class TestUnitFlash(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import flash
        return flash( *args, **kwargs )

    def test_runs_with_options( self ):
        ''' -o outprefix -d outdir '''
        efiles = self.flash_output_files( 'test', 'flashout' )
        flashout = self._C( *self.fqs, o='test', d='flashout' )
        eq_( set([]), set(efiles) - set(flashout) )
        for f in flashout:
            ok_( exists( f ), 'Did not create {0}'.format(f) )

    def test_runs_default( self ):
        efiles = self.flash_output_files( 'out', '' )
        flashout = self._C( *self.fqs )
        eq_( set([]), set(efiles) - set(flashout) )
        for f in flashout:
            ok_( exists( f ), 'Did not create {0}'.format(f) )

    def test_fail_run( self ):
        eq_( [], self._C( '1', '2' ), 'Did not return empty list on bad run' )

class TestUnitBtrim(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import btrim
        return btrim( *args, **kwargs )

    def test_runs( self ):
        from bactpipeline.runsample import flash
        flashout = flash( *self.fqs )
        inf = flashout[0]
        outf = inf.replace('.fastq','.btrim.fastq')
        r = self._C( p=self.truseq, t=inf, o=outf, b=300, c=True, Q=True, S=True, l=100 )
        ok_( exists(outf), 'btrim did not create {0}'.format(outf) )
        ok_( os.stat(inf).st_size != os.stat(outf).st_size, 'Did not trim input file. Input and output file are same size still' )

    def test_runs_outdir_noexist( self ):
        # btrim does not create full path to file
        inf = self.flash()[0]
        outf = join('btrim','outfile.fastq')
        r = self._C( p=self.truseq, t=inf, o=outf )
        ok_( 0 != r, 'btrim return code was {0} instead of non-zero'.format(r) )

    def test_returns_retcode( self ):
        from bactpipeline.runsample import flash
        flashout = flash( *self.fqs )
        eq_( 0, self._C( p=self.truseq, t=flashout[0], o='somefile' ), 'Did not return return code' )
        r = self._C()
        ok_( 0 != r, 'btrim return code was {0} instead of non-zero'.format(r) )

class TestUnitBtrimFiles( Base ):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import btrim_files
        return btrim_files( *args, **kwargs )

    def test_creates_outdir( self ):
        outp = join('dir1','dir2')
        r = self._C( self.flash(), outp, p=self.truseq )
        ok_( exists( outp ), 'Did not create {0}'.format(outp) )

    def test_runs_multiple( self ):
        outp = 'btrim'
        inf = self.flash()
        r = self._C( inf, outp, p=self.truseq )
        bfiles = []
        for f in inf:
            ff = join( outp, basename(f).replace('.fastq','.btrim.fastq') )
            bfiles.append(ff)
            ok_( exists( ff ), 'Did not create {0}'.format(ff) )
        eq_( set([]), set(bfiles)-set(r) )

class TestUnitRunAssembly(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import run_assembly
        return run_assembly( *args, **kwargs )

    def test_creates_and_runs_project( self ):
        r = self._C( self.fqs )
        eq_( 0, r, 'runAssembly return non-zero' )
        projdir = glob( 'P_*' )
        eq_( 1, len(projdir), 'Did not create a project' )
        assxml = join(projdir[0],'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {0}'.format(projdir) )

    def test_creates_named_project( self ):
        projdir = 'projdir'
        r = self._C( self.fqs, o=projdir )
        eq_( 0, r, 'runAssembly return non-zero' )
        assxml = join(projdir,'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {0}'.format(projdir) )

class TestUnitNewAssembly(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import new_assembly
        return new_assembly( *args, **kwargs )

    def test_projdir_set( self ):
        projdir = 'testproject'
        r = self._C( projdir )
        assxml = join(projdir,'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {0}'.format(projdir) )

    def test_project_notset( self ):
        r = self._C()
        projdir = glob( 'P_*' )
        eq_( 1, len(projdir), 'Did not create a project' )
        eq_( projdir[0], r, 'Did not return project path correctly. {0} != {1}'.format(projdir,r) )

class TestUnitReplaceNewblerSettings(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import replace_newbler_settings
        return replace_newbler_settings( *args, **kwargs )

    def test_returns_correct_xml( self ):
        from bactpipeline.runsample import new_assembly
        # Create new assembly proj
        projdir = new_assembly( )
        xmlfile = join( projdir, 'assembly', '454AssemblyProject.xml' )
        r = self._C( projdir, self.fqs )
        xml = open(xmlfile).read()
        print(xml)
        eq_( 2, xml.count('ReadFiles'), 'Did not insert ReadFiles correctly' )
        self.fqs_in_xml( self.fqs, xml )

class TestUnitNewblerFastqReadFiles(Base):
    def _C( self, *args, **kwargs ):
        from bactpipeline.runsample import newbler_fastq_read_files
        return newbler_fastq_read_files( *args, **kwargs )

    def test_correct_xml( self ):
        r = self._C( self.fqs )
        self.fqs_in_xml( self.fqs, r )

class TestFunctional(Base):
    def _C( self, *args, **kwargs ):
        script = kwargs.get('script','runsample')
        cmd = [script]
        cmd += ['-o', kwargs.get('o','output')]
        cmd += list(args)
        print(' '.join(cmd))
        r = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        r.wait()
        return r

    def _dump_sout_serr(self, sout, serr):
        print("STDOUT: ")
        print(sout)
        print("STDERR: ")
        print(serr)

    @attr('current')
    def test_writes_summary_and_contig(self):
        r = self._C(self.fixdir)
        self._dump_sout_serr(*r.communicate())
        sumfile = join('output','summary.tsv')
        confile = join('output','top_contigs.fasta')
        ok_(exists(sumfile), 'Did not create summary.tsv')
        ok_(exists(confile), 'Did not create top_contigs.fasta')
        eq_(2, len(open(confile).readlines()))
        stats = open(sumfile).readlines()[1].split()
        eq_(['fix_fastq','10','1','1170','100.0','9'], stats)
        eq_(0, r.returncode)

    def test_fixed_fastqfiles( self ):
        r = self._C( self.fixdir )
        self._dump_sout_serr(*r.communicate())
        outdir = join( 'output', 'fix_fasta' )
        ok_( exists(outdir), 'Did not create fix_fasta output dir' )
        rfqs = os.listdir( join('output','fix_fasta') )
        eq_( set([]), set(self.bnfqs) - set(rfqs) )
        eq_(0, r.returncode)

    def test_flash_files( self ):
        r = self._C( self.fixdir )
        print(r.communicate())
        outdir = join( 'output', 'flash' )
        prefix = 'out'
        efiles = self.flash_output_files( prefix, outdir )
        for f in efiles:
            ok_( exists(f), 'Did not create flash output file {0}'.format(f) )
        eq_(0, r.returncode)

    def test_btrim_files( self ):
        r = self._C( self.fixdir )
        print(r.communicate())
        flasho = filter( lambda x: x.endswith('.fastq'), self.flash_output_files( 'out', 'flash' ) )
        outdir = join( 'output', 'btrim' )
        print([join(root,f) for root,dirs,files in os.walk('.') for f in files])
        for f in flasho:
            bfile = join(outdir,basename(f)).replace('.fastq','.btrim.fastq')
            ok_( exists(bfile), 'Did not create btrim output file {0}'.format(bfile) )
        eq_(0, r.returncode)

    def test_newbler_files( self ):
        r = self._C( self.fixdir, o='output' )
        print(r.communicate())
        projdir = join( 'output', 'newbler_assembly' )
        xmlp = join( projdir, 'assembly', '454AssemblyProject.xml' )
        xml = open(xmlp).read()
        btrimo = join( 'output', 'btrim' )
        btrimfqs = glob( join(btrimo,'*.fastq') )
        self.fqs_in_xml( btrimfqs, xml )
        eq_(0, r.returncode)

from bactpipeline.runsample import N50, N_stat
class StatsTests(unittest.TestCase):
    def setUp(self):
        cs = '''
        GATTACA
        TACTACTAC
        ATTGAT
        GAAGA
        '''
        self.lengths = list(map(len, filter(bool, cs.split())))

    def test_n50(self):
         self.assertEquals(7, N50(self.lengths))

    def test_n75(self):
        self.assertEquals(6,  N_stat(self.lengths, 0.75))

from bactpipeline.runsample import read_count, InvalidFastqException
@attr('current')
class TestReadCount(common.Base):
    def test_counts_correctly(self):
        r = read_count(fqs)
        self.assertEqual(2000, r)

    def test_fq_file_not_4_lines_raises_error(self):
        p = join(self.tdir, 'bad.fastq')
        with open(p, 'w') as fw:
            fw.write('@id1\n')
            fw.write('ATGC\n')
            fw.write('+\n')
        self.assertRaises(InvalidFastqException, read_count, [p])
