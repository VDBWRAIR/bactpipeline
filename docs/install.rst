=======
Install
=======

Requirements
============

* Roche 454 Analysis Software(gsAssembler, runProject)
* Everything should be bundled with the pipeline

Installation
============

1. Clone repo

.. code-block:: bash

    git clone https://github.com/VDBWRAIR/bactpipeline.git

2. Enter code directory

.. code-block:: bash

    cd bactpipeline

3. Install python virtualenv

.. code-block:: bash

    virtualenv env
    . env/bin/activate

4. Install bactpipeline

.. code-block:: bash

    pip install -r requirements.txt
    python setup.py install
