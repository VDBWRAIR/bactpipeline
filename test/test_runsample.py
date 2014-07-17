from imports import *
import common

class Base( common.Base ):
    fixdir = join( dirname(__file__), 'fixtures', 'fix_fastq' )
    fqs = glob( join( fixdir, '*.fastq' ) )
    bnfqs = [basename(fq) for fq in fqs]
    truseq = join(dirname(dirname(__file__)),'truseq.txt')
    assemprojxml = join(fixdir,'454AssemblyProject.xml')
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

    def flash( self, prefix='out', outdir='flash' ):
        from runsample import flash
        # Files from flash
        return [f for f in flash( *self.fqs, d=outdir, o=prefix ) if f.endswith('.fastq')]

    def fqs_in_xml( self, fqlist, xml ):
        for fqfile in fqlist:
            ok_( fqfile in xml, '{} did not make it into 454AssemblyProject.xml'.format(fqfile) )
        ok_( '<SetReadCount>{}</SetReadCount>'.format(len(self.fqs)), 'Did not set read count correctly' )

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

    def test_runs_default( self ):
        efiles = self.flash_output_files( 'out', '' )
        flashout = self._C( *self.fqs )
        eq_( set([]), set(efiles) - set(flashout) )
        for f in flashout:
            ok_( exists( f ), 'Did not create {}'.format(f) )

    def test_fail_run( self ):
        eq_( [], self._C( '1', '2' ), 'Did not return empty list on bad run' )

class TestUnitBtrim(Base):
    def _C( self, *args, **kwargs ):
        from runsample import btrim
        return btrim( *args, **kwargs )

    def test_runs( self ):
        from runsample import flash
        flashout = flash( *self.fqs )
        inf = flashout[0]
        outf = inf.replace('.fastq','.btrim.fastq')
        r = self._C( p=self.truseq, t=inf, o=outf, b=300, c=True, Q=True, S=True, l=100 )
        ok_( exists(outf), 'btrim did not create {}'.format(outf) )
        ok_( os.stat(inf).st_size != os.stat(outf).st_size, 'Did not trim input file. Input and output file are same size still' )

    def test_runs_outdir_noexist( self ):
        # btrim does not create full path to file
        inf = self.flash()[0]
        outf = join('btrim','outfile.fastq')
        r = self._C( p=self.truseq, t=inf, o=outf )
        ok_( 0 != r )

    def test_returns_retcode( self ):
        from runsample import flash
        flashout = flash( *self.fqs )
        eq_( 0, self._C( p=self.truseq, t=flashout[0], o='somefile' ), 'Did not return return code' )
        ok_( 0 != self._C() )

class TestUnitBtrimFiles( Base ):
    def _C( self, *args, **kwargs ):
        from runsample import btrim_files
        return btrim_files( *args, **kwargs )

    def test_creates_outdir( self ):
        outp = join('dir1','dir2')
        r = self._C( self.flash(), outp, p=self.truseq )
        ok_( exists( outp ), 'Did not create {}'.format(outp) )

    def test_runs_multiple( self ):
        outp = 'btrim'
        inf = self.flash()
        r = self._C( inf, outp, p=self.truseq )
        bfiles = []
        for f in inf:
            ff = join( outp, basename(f).replace('.fastq','.btrim.fastq') )
            bfiles.append(ff)
            ok_( exists( ff ), 'Did not create {}'.format(ff) )
        eq_( set([]), set(bfiles)-set(r) )

class TestUnitRunAssembly(Base):
    def _C( self, *args, **kwargs ):
        from runsample import run_assembly
        return run_assembly( *args, **kwargs )
    
    @attr('current')
    def test_creates_and_runs_project( self ):
        r = self._C( self.fqs )
        eq_( 0, r, 'runAssembly return non-zero' )
        projdir = glob( 'P_*' )
        eq_( 1, len(projdir), 'Did not create a project' )
        assxml = join(projdir[0],'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {}'.format(projdir) )

    @attr('current')
    def test_creates_named_project( self ):
        projdir = 'projdir'
        r = self._C( self.fqs, o=projdir )
        eq_( 0, r, 'runAssembly return non-zero' )
        assxml = join(projdir,'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {}'.format(projdir) )

class TestUnitNewAssembly(Base):
    def _C( self, *args, **kwargs ):
        from runsample import new_assembly
        return new_assembly( *args, **kwargs )

    def test_projdir_set( self ):
        projdir = 'testproject'
        r = self._C( projdir )
        assxml = join(projdir,'assembly','454AssemblyProject.xml')
        ok_( exists(assxml), 'Did not create project {}'.format(projdir) )

    def test_project_notset( self ):
        r = self._C()
        projdir = glob( 'P_*' )
        eq_( 1, len(projdir), 'Did not create a project' )
        eq_( projdir[0], r, 'Did not return project path correctly. {} != {}'.format(projdir,r) )

class TestUnitReplaceNewblerSettings(Base):
    def _C( self, *args, **kwargs ):
        from runsample import replace_newbler_settings
        return replace_newbler_settings( *args, **kwargs )

    def test_returns_correct_xml( self ):
        from runsample import new_assembly
        # Create new assembly proj
        projdir = new_assembly( )
        xmlfile = join( projdir, 'assembly', '454AssemblyProject.xml' )
        r = self._C( projdir, self.fqs )
        xml = open(xmlfile).read()
        print xml
        eq_( 2, xml.count('ReadFiles'), 'Did not insert ReadFiles correctly' )
        self.fqs_in_xml( self.fqs, xml )

class TestUnitNewblerFastqReadFiles(Base):
    def _C( self, *args, **kwargs ):
        from runsample import newbler_fastq_read_files
        return newbler_fastq_read_files( *args, **kwargs )

    def test_correct_xml( self ):
        r = self._C( self.fqs )
        self.fqs_in_xml( self.fqs, r )

class TestFunctional(Base):
    def _C( self, *args, **kwargs ):
        script = kwargs.get('script',join( dirname(dirname(__file__)), 'bin', 'runsample.py' ))
        cmd = [script]
        cmd += ['-o', kwargs.get('o','output')]
        cmd += list(args)
        return subprocess.call( cmd )

    @attr('current')
    def test_script_path( self ):
        r = self._C( self.fixdir, script=join(dirname(dirname(__file__)), 'runsample.py') )
        eq_( 0, r, 'Did not execute command correctly with actual script path' )

    @attr('current')
    def test_script_symlinkpath( self ):
        r = self._C( self.fixdir, script=join(dirname(dirname(__file__)), 'bin', 'runsample.py') )
        eq_( 0, r, 'Did not execute command correctly with symlink script path' )

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

    def test_btrim_files( self ):
        r = self._C( self.fixdir )
        flasho = filter( lambda x: x.endswith('.fastq'), self.flash_output_files( 'out', 'flash' ) )
        outdir = join( 'output', 'btrim' )
        print [join(root,f) for root,dirs,files in os.walk('.') for f in files]
        for f in flasho:
            bfile = join(outdir,basename(f)).replace('.fastq','.btrim.fastq')
            ok_( exists(bfile), 'Did not create btrim output file {}'.format(bfile) )

    def test_newbler_files( self ):
        r = self._C( self.fixdir, o='output' )
        projdir = join( 'output', 'newbler_assembly' )
        xmlp = join( projdir, 'assembly', '454AssemblyProject.xml' )
        xml = open(xmlp).read()
        btrimo = join( 'output', 'btrim' )
        btrimfqs = glob( join(btrimo,'*.fastq') )
        self.fqs_in_xml( btrimfqs, xml )
