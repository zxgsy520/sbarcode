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


def get_sample(file):

    name = file.split('/')[-1]

    if '--' in name:
        name = name.split('--')[1].split('.bam')[0]
    else:
        name = name.split('.')[0]

    return name


def filter_bam(file, qvalue=0.9):

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    name = get_sample(file)
    fo = pysam.AlignmentFile("%s.clean.bam" % name, "wb", template=fh)

    for line in fh:
        if line.get_tag('rq')<qvalue:
            continue
        fo.write(line)

    fh.close()
    fo.close()


def filter_bams(files, qvalue):

    for file in files:
        filter_bam(file, qvalue)


def add_hlep_args(parser):

    parser.add_argument('-i', '--input', nargs='+', metavar='FILE', type=str, required=True,
        help='Input reads file, format(bam and sam).')
    parser.add_argument('-q', '--qvalue', metavar='FLOAT', type=float, default=0.9,
        help='input filtered quality value, default=0.9.')

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
    filter_ccs_bam.py Filter ccs bam files based on quality values.

attention:
    filter_ccs_bam.py -i *.bam
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    filter_bams(args.input, args.qvalue)


if __name__ == "__main__":

    main()
