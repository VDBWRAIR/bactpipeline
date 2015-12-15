from setuptools import setup, find_packages

import bactpipeline

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
            'runsample = bactpipeline.runsample:main'
        ],
    },
)
