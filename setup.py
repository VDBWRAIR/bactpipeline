from __future__ import print_function

from setuptools import setup, find_packages
import subprocess
import sys
from os.path import join, exists, basename
import shutil

import bactpipeline

btrim_bin = 'bactpipeline/dependencies/btrim/btrim64-static'
flash_path = 'bactpipeline/dependencies/FLASH-1.2.10'
flash_bin = join(flash_path, 'flash')
if not exists(flash_bin):
    try:
        subprocess.check_call('make', cwd=flash_path)
    except subprocess.CalledProcessError as e:
        print("Failed to build flash prior to install.")
        print("Likely you require some system packages in order to compile it")
        sys.exit(1)

setup(
    name = bactpipeline.__projectname__,
    version = bactpipeline.__release__,
    packages = find_packages(),
    author = bactpipeline.__authors__,
    author_email = bactpipeline.__authoremails__,
    description = bactpipeline.__description__,
    license = "GPLv2",
    keywords = bactpipeline.__keywords__,
    entry_points = {
        'console_scripts': [
            'runsample = bactpipeline.runsample:main',
            'run_bactpipeline = bactpipeline.runsample:luigi_run',
            'fix_fastq = bactpipeline.fix_fastq:main',
            'flash = bactpipeline.util:flash',
            'btrim = bactpipeline.util:btrim',
        ],
    },
    package_data = {
        'bactpipeline': [
            'truseq.txt',
            # Python3 won't copy these if they are in scripts
            # because of encoding issue
            btrim_bin.replace('bactpipeline/',''),
            flash_bin.replace('bactpipeline/','')
        ]
    }
)
