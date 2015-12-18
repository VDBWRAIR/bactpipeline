import pkg_resources
import sys
import os

def flash():
    _exec('dependencies/FLASH-1.2.10/flash')

def btrim():
    _exec('dependencies/btrim/btrim64-static')

def _exec(name):
    _bin = pkg_resources.resource_filename(__name__, name)
    os.execvp(_bin, sys.argv)
