from imports import *

class Base( object ):
    def setUp( self ):
        self.tdir = tempfile.mkdtemp( suffix='test' )
        os.chdir( self.tdir )

    def tearDown( self ):
        os.chdir( '/' )
        shutil.rmtree( self.tdir )
