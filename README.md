Sbarcode
==============
### Version: 1.0.0
Split, process and count pacbio data

Third-party
-----------
Sbarcode package includes some third-party software:
* [python](https://www.python.org/)
* [pysam](https://pypi.org/project/pysam/)

Building Requirements
-----------
* Python 2.7 or 3.5+ (with setuptools package installed)

Local building (without installation)
-----------

You may use the package locally without system installation.
To get and compile the latest git version, run:
<pre><code>
git clone https://github.com/zxgsy520/sbarcode
cd sbarcode
</code></pre>

## Instructions
Modify prefix names of files in batch
<pre><code>
rename_lima_file.py -h
usage: rename_lima_file.py [-h] -f FILE [FILE ...] -b STR [-p FILE]

name:
    rename_lima_file.py  Rename the file name after lima split

barcode_name.txt:
sample  R  F
or
sample  barcode_id

attention:
    rename_lima_file.py --files * --bname barcode_name.txt
version: 1.2.0
contact:  Xingguo Zhang <113178210@qq.com>        

optional arguments:
  -h, --help            show this help message and exit
  -f FILE [FILE ...], --files FILE [FILE ...]
                        Input files.
  -b STR, --bname STR   Input sample and barcode number alignment file.
  -p FILE, --path FILE  The path where the input file needs to be moved.
</code></pre>
Convert bam files to fastq files in batch
<pre><code>
bams2fqs.py -h
usage: bams2fqs.py [-h] FILE [FILE ...]

URL: https://github.com/zxgsy520/sbarcode
name:
    bams2fqs.py Convert bam file to fastq file.
     
attention:
    bams2fqs.py *.bam

version: 1.0.0
contact:  Xingguo Zhang <113178210@qq.com>        

positional arguments:
  FILE        Input reads file, format(bam and sam).

optional arguments:
  -h, --help  show this help message and exit
</code></pre>
Run related scripts
<pre><code>
python lib/rename_lima_file.py --files * --bname barcode_name.txt
python lib/bams2fqs.py *.bam
</code></pre>


  
