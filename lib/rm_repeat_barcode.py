#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line)
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def complement(seq):

    cdict = {"A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    seq = list(seq.upper())
    nseq = ""
    for i in seq:
        nseq += cdict[i]

    return nseq


def reverse_complement(seq):

    seq = seq[::-1]

    return complement(seq)


def rm_repeat_barcode(file):

    seqs = {}

    for seqid, seq in read_fasta(file):
        rcseq = reverse_complement(seq)
        if seq not in seqs and rcseq not in seqs:
            print(">%s\n%s" % (seqid, seq))
        else:
            if seq in seqs:
                repeatid = seqs[seq]
            else:
                repeatid = seqs[rcseq]
            LOG.info("%s = %s" % (seqid, repeatid ))
        seqs[rcseq] = seqid
        seqs[seq] = seqid

    return 0


def add_hlep_args(parser):

    parser.add_argument('fasta', metavar='FILE', type=str,
        help='Input barcode sequence')

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
    rm_repeat_barcode.py: Remove repeat barcode sequences
attention:
    rm_repeat_barcode.py barcode.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    rm_repeat_barcode(args.fasta)


if __name__ == "__main__":

    main()
