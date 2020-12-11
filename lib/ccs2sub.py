#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import pysam
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def get_sample(file):

    name = file.split('/')[-1]

    if '--' in name:
        name = name.split('--')[1].split('.bam')[0]
    else:
        name = name.split('.')[0]

    return name


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []

    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq
    fp.close()


def read_bam(file):

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        yield [line.qname, line.seq]

    fh.close()


def read_ccs(file):

    reads = set()

    if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fa.gz") or file.endswith(".fasta.gz"):
        fh = read_fasta(file)
    elif file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fq.gz") or file.endswith(".fastq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".bam") or file.endswith(".sam"):
        fh = read_bam(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        seqid = line[0]
        ccsid = line[0].split('/ccs')[0]

        reads.add(ccsid)

    return reads


def ccs2sub(file, ccs, output):

    reads = read_ccs(ccs)
    name = get_sample(ccs)

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    fname = os.path.abspath(output)

    if os.path.exists(fname):
        os.remove(fname)

    fo = pysam.AlignmentFile(fname, "wb", template=fh)

    for line in fh:
        seqid = line.qname.split('/')
        ccsid = '/'.join(seqid[0:2])

        if ccsid not in reads:
            continue
        fo.write(line)

    fh.close()
    fo.close()


def add_hlep_args(parser):

    parser.add_argument('subreads',
        help='Input subreads file, format(bam and sam).')
    parser.add_argument('-c', '--ccs', metavar='FILE', type=str, required=True,
        help='Input the CCS reads file, format(bam, sam, fastq, fasta).')
    parser.add_argument('-o', '--out', metavar='FILE', type=str, default='out.subreads.bam',
        help='Output file.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    ccs2sub.py Obtain subreads according to CCS reads.

attention:
    ccs2sub.py subreads.bam --ccs ccs.fa --out out.subreads.bam
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    ccs2sub(args.subreads, args.ccs, args.out)


if __name__ == "__main__":

    main()
