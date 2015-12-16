from setuptools import setup, find_packages

import bactpipeline
import subprocess
import sys
from os.path import join, exists

flash_path = 'dependencies/FLASH-1.2.10'
flash_bin = join(flash_path, 'flash')
if not exists(flash_bin):
    try:
        subprocess.check_call('make', cwd=flash_path)
    except subprocess.CalledProcessError as e:
        print "Failed to build flash prior to install."
        print "Likely you require some system packages in order to compile it"
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
    scripts = [
        'dependencies/FLASH-1.2.10/flash',
        'dependencies/btrim/btrim64-static',
    ],
    entry_points = {
        'console_scripts': [
            'runsample = bactpipeline.runsample:main',
            'fix_fastq = bactpipeline.fix_fastq:main',
        ],
    },
    package_data = {
        'bactpipeline': ['truseq.txt']
    }
)
