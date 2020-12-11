#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import json
import gzip
import pysam
import logging
import argparse

import collections
import numpy as np

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


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
        yield [line.qname, line.seq, line.get_tag('rq')]

    fh.close()


def n50(lengths):

    sum_length = sum(lengths)
    accu_length = 0

    for i in sorted(lengths, reverse=True):
        accu_length += i

        if accu_length >= sum_length*0.5:
            return i


def read_file(file):

    if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fa.gz") or file.endswith(".fasta.gz"):
        fh = read_fasta(file)
    elif file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fq.gz") or file.endswith(".fastq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".bam") or file.endswith(".sam"):
        fh = read_bam(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        rq = 1
        if len(line)<3:
            line.append(rq)
        yield line


def stat_quality(file, model="model"):

    data = {}
    total = 0
    length = []

    for line in read_file(file):
        total += 1
        length.append(len(line[1]))
        seqid = "%s [rq=%s]" % (line[0], line[2])
        data[seqid] = [line[1], line[2]]

    if model=="model":
        counts = np.bincount(length)
        mean_len = np.argmax(counts)
    elif model=="median":
        mean_len = np.median(length)
    else:
        mean_len = np.mean(length)

    return data, total, mean_len


def get_sample(file):

    name = file.split('/')[-1]

    if '--' in name:
        name = name.split('--')[1].split('.bam')[0]
    else:
        name = name.split('.')[0]

    return name


def stat_reads(files, out, about, qvalue, model):

    title = ["Sample", "Total_CCS", "Q20_CCS", "Q30_CCS", "Q40_CCS", "Average_length(bp)", "Effective_Rate(%)"]
    sdata = collections.OrderedDict()
    sample = collections.OrderedDict()
    fo = open(out, 'w')
    fo.write('Sample\tTotal_CCS\tQ20_CCS\tQ30_CCS\tQ40_CCS\tAverage_length(bp)\tEffective_Rate(%)\n')
    for file in files:
        name = get_sample(file)
        sdata[name] = collections.OrderedDict()
        sample[name] = [name, "", ""]
        data, total, mean_len = stat_quality(file, model)
        fa = open("%s.clean.fasta" % name, "w")
        q20 = 0
        q30 = 0
        q40 = 0
        length = []

        for seqid in data:
            line = data[seqid]
            if mean_len+about<len(line[0]) or mean_len-about>len(line[0]):
                continue
            if line[1]<qvalue:
                continue
            length.append(len(line[0]))
            fa.write(">%s\n%s\n" % (seqid, line[0]))

            q20 += 1
            if line[1]>=0.999:
                q30 += 1
            if line[1]>=0.9999:
                q40 += 1

        fa.close()

        fo.write('{0}\t{1:,}\t{2:,}\t{3:,}\t{4:,}\t{5:,.2f}\t{6:.2f}\n'.format(name, total, q20, q30, q40, sum(length)*1.0/len(length), q20*100.0/total))
        line = [name, total, q20, q30, q40, round(sum(length)*1.0/len(length), 2), round(q20*100.0/total, 2)]
        for i in range(len(title)-1):
            sdata[name][title[i+1]] = line[i+1]
    fo.close()
    print("field = %s" % json.dumps(title))
    print("summary = %s" % json.dumps(sdata))
    print("sample_lib_lane = %s" % json.dumps(sample))


def add_hlep_args(parser):

    parser.add_argument('-i', '--input', nargs='+', metavar='FILE', type=str, required=True,
        help='Input reads file, format(fasta,fastq,fa.gz,bam and sam).')
    parser.add_argument('-a', '--about', metavar='INT', type=int, default=3,
        help='Allowable read length error range, default=3bp.')
    parser.add_argument('-q', '--qvalue', metavar='FLOAT', type=float, default=0.99,
        help='input filtered quality value, default=0.99.')
    parser.add_argument('-m', '--model',  choices=["model", "median", "mean"], default="model",
        help='Select the reference length for filtering, default=model.')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default="out.tsv",
        help='Out name.')

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
    filter_barcode.py Statistics filtering pacbio ccs reads

attention:
    filter_barcode.py -i *.bam
    filter_barcode.py -i *.fa
    filter_barcode.py -i *.fq
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_reads(args.input, args.out, args.about, args.qvalue, args.model)


if __name__ == "__main__":

    main()
