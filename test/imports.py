import os
from os.path import *
import sys
import tempfile
from glob import glob
import subprocess
import shutil
import shlex
from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr
