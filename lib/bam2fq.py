#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import pysam
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_bam(file):

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        #yield [line.qname, line.seq, pysam.array_to_qualitystring(line.query_qualities), line.get_tag('rq')]
        yield line.qname, line.seq, pysam.array_to_qualitystring(line.query_qualities)

    fh.close()


def bam2fq(file):

    for seqid, seq, quality  in read_bam(file):
        print('@%s\n%s\n+\n%s' % (seqid, seq, quality))


def add_hlep_args(parser):

    parser.add_argument('bam',
        help='Input reads file, format(bam and sam).')

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
    bam2fq.py Convert bam file to fastq file.

attention:
    bam2fq.py input.bam >ouput.fq

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    bam2fq(args.bam)


if __name__ == "__main__":

    main()
