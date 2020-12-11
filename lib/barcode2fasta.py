#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_barc_txt(file):

    data = {}

    for line in read_tsv(file, '\t'):
        if line[0] in data:
            LOG.info('Sample %s has duplicates, please check the input' % line[0])
            continue
        if '-' in line[1]:
            prefix, number = line[1].upper().split('-')
            name = "%s-%02d" % (prefix, int(number))
        else:
            name = line[1]

        data[name] = line[0]

    return data


def barcode2fasta(file, database):

    data = read_barc_txt(file)

    for file in database:
        for line in read_tsv(file, '\t'):
            name = line[1].strip().strip('F')

            if name in data:
                print(">%s\n%s" % (data[name],  line[2].strip()))
            elif line[0] in data:
                print(">%s\n%s" % (data[line[0]],  line[2].strip()))
            else:
                continue

    return 0
              


def add_hlep_args(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Correspondence table of input sample and barcode id.')
    parser.add_argument('-d','--database', nargs='+', metavar='FILE', type=str,
        default="/export/personal/software/Pipeline/sbarcode/v1.0.0/database/*.barcode.txt",
        help='Input barcode database.')

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
    barcode2fasta.py Generate barcode sequence

attention:
    barcode2fasta.py -i barcode.txt >barcode.fa
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    barcode2fasta(args.input, args.database)


if __name__ == "__main__":

    main()
