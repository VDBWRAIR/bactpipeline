# bactpipeline
Bacteria Pipeline modified from MRSN Bacteria Pipeline

## Install

### Requirements

* Roche 454 Analysis Software(gsAssembler, runProject)
* Everything should be bundled with the pipeline

### Installation

1. Clone repo

```
cd ~/ && git clone https://github.com/VDBWRAIR/bactpipeline.git
```

2. Install python virtualenv

```
virtualenv ~/bactpipeline
. ~/bactpipeline/bin/activate
pip install -r ~/bactpipeline/requirements.txt
```

## Running the Pipeline

If you have not already, make sure to check out the [[install]] page first.

You will have to re-source the virtualenv when you want to use the pipeline, 
using `. ~/bactpipeline/bin/activate`

### Running a single sample

Running a single sample is fairly straight forward assuming you have your miseq 
reads parsed out into individual folders for each sample.

#### Usage Example

For this example, assume you have a directory /home/username/reads/sample1 that 
contains MiSeq reads for sample1:

* sample1_S1_L001_R1_001.fastq
* sample1_S1_L001_R2_001.fastq

Simply run the following command to run the pipeline on sample1's data
```
runsample.py -o sample1 /home/username/reads/sample1
```

Multiple samples can be run with the --sample-sheet parameter.
```
./runsample.py --sample-sheet samplesheet.csv -o outdir
./runsample.py -s samplesheet.csv -o outdir
```
See `samplesheet.csv` for the fields that are required in the samplesheet.

## runsample.py

This script takes care of putting all the pieces of the pipeline together

### Pipeline Flow

* fix_fastq.py
* flash
* btrim
* newbler(runAssembly)

### Usage

* Output Directory
   * -o or --output
   * Specifies where to put all the resulting output directories for the various stages
   * Default: output
* readdir
   * Specifies a directory that contains the paired MiSeq reads

```
runsample.py [-o|--output output] readdir
```

### Output Files

* fix_fastq
   * fastq files with same name as were in the readdir argument but with sequence id modified for Newbler
* flash/
   * out.extendedFrags.fastq
     * paired reads combined together
   * out.notCombined_1.fastq
     * R1 reads that did not combine
   * out.notCombined_2.fastq
     * R2 reads that did not combine
   * out.hist
     * Combined read lengths
   * out.histogram
     * Combined read lengths visual
* btrim
    * fastq files with same name as out.*.fastq from flash, but with .btrim.fastq at end
* newbler_assembly

## fix_fastq.py

"Original Perl Version":/attachments/download/25/fixFastq.pl
This script handles renaming sequence identifiers in Illumina reads such that Newbler will use them as paired end correctly.

It addresses [this](http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header)

### Usage

```
fix_fastq.py [-o outdir] fastq [fastq ...]
```

### Example usage

You essentially supply the script with the location of any fastq files you want and it will replace the sequence id in each and copy the modified version into an output directory.

If you have a bunch of fastq files in a directory, lets say /home/username/reads, then you could run it as follows:
```
fix_fastq.py -o newbler_reads /home/username/reads/*.fastq
```

All modified reads would then be placed in a directory called newbler_reads in the current directory you are in
** Newbler assembly project that can be opened by gsAssembler application. See Newbler 2.8 documentation for more information.

