from __future__ import print_function
from imports import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest

class Base(unittest.TestCase):
    tempdir = join(TEST_DIR, 'testruns')

    def setUp( self ):
        # Ensure test tempdir exists
        try:
            os.mkdir(self.tempdir)
        except OSError as e:
            pass
        self.tdir = tempfile.mkdtemp(suffix='test', dir=self.tempdir)
        print(self.tdir)
        os.chdir( self.tdir )

    def tearDown( self ):
        pass
        #os.chdir( '/' )
        #shutil.rmtree( self.tdir )
