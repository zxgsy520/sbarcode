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
        yield [line.qname, line.seq, pysam.array_to_qualitystring(line.query_qualities), line.get_tag('rq')]

    fh.close()


def filter_bam2fq(file, qvalue=0.9, minlen=500, maxlen=10000):

    for line in read_bam(file):
        if line[3] <qvalue:
             continue
        if len(line[1])<minlen or len(line[1])>maxlen:
             continue
        print('@%s rq=%s\n%s\n+\n%s' % (line[0], line[3], line[1], line[2]))


def add_hlep_args(parser):

    parser.add_argument('bam', 
        help='Input reads file, format(bam and sam).')
    parser.add_argument('-q', '--qvalue', metavar='FLOAT', type=float, default=0.9,
        help='Input filtered quality value, default=0.9.')
    parser.add_argument('--minlen', metavar='INT', type=int, default=500,
        help='Input the minimum length of the filter, default=500.')
    parser.add_argument('--maxlen', metavar='INT', type=int, default=10000,
        help='Input the maximum length of the filter, default=10000.')

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
    filter_bam2fq.py Filter ccs bam and output fastq file

attention:
    filter_bam2fq.py input.bam
    filter_bam2fq.py input.bam --qvalue 0.9 --minlen 500 --maxlen 10000
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    filter_bam2fq(args.bam, args.qvalue, args.minlen, args.maxlen)


if __name__ == "__main__":

    main()
