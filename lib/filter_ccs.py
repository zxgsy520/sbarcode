#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import json
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


def read_bam(file):

    if file.endswith(".bam"):
        fh = pysam.AlignmentFile(file, "rb", check_sq=False)
    elif file.endswith(".sam"):
        fh = pysam.AlignmentFile(file, 'r')
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        #line.get_tags('DB')
        yield [line.qname, line.seq, pysam.array_to_qualitystring(line.query_qualities), line.get_tag('rq')]

    fh.close()


def n50(lengths):

    sum_length = sum(lengths)
    accu_length = 0

    for i in sorted(lengths, reverse=True):
        accu_length += i

        if accu_length >= sum_length*0.5:
            return i


def stat_quality(file, model="model"):

    data = {}
    total = 0
    length = []

    for line in read_bam(file):
        total += 1
        length.append(len(line[1]))
        seqid = "%s [rq=%s]" % (line[0], line[3])
        data[seqid] = [line[1], line[2], line[3]]

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
        fa = open("%s.clean.fastq" % name, "w")
        q20 = 0
        q30 = 0
        q40 = 0
        length = []

        for seqid in data:
            line = data[seqid]
            if mean_len+about<len(line[0]) or mean_len-about>len(line[0]):
                continue
            if line[2]<qvalue:
                continue
            length.append(len(line[0]))
            fa.write("@%s\n%s\n+\n%s\n" % (seqid, line[0], line[1]))

            q20 += 1
            if line[2]>=0.999:
                q30 += 1
            if line[2]>=0.9999:
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
        help='Input reads file, format(bam and sam).')
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
    filter_ccs.py Statistics filtering pacbio ccs reads

attention:
    filter_ccs.py -i *.bam
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_reads(args.input, args.out, args.about, args.qvalue, args.model)


if __name__ == "__main__":

    main()
