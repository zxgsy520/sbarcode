#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
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


def move_file(file, nfile):

    LOG.info("mv %s %s" % (file, nfile))

    shutil.move(file, nfile)

    return 0


def copy_file(file, nfile):

    LOG.info("mv %s %s" % (file, nfile))

    shutil.copy(file, nfile)

    return 0


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def find_file_name(file):

    name = file.split('/')[-1]
    rfp = name.split('.')[1]
    rp, fp = rfp.split('--')

    return {rp, fp}


def read_barcode_name(file):

    data = {}

    for line in read_tsv(file, '\t'):
        if len(line) >= 3:
            data[line[0]] = {line[1], line[2]}
        else:
            data[line[0]] = {line[1], line[1]}

    return data


def judge_primer(rfp, data):
    
    educt = 0

    for i in data:
        if rfp==data[i]:
           educt = i
           break

    return educt


def rename_lima_file(files, bname, path):

    data = read_barcode_name(bname)

    for file in files:
        rfp = find_file_name(file)
        name = judge_primer(rfp, data)

        if not name:
            continue

        forms = file.split('/')[-1]
        forms = forms.split('.')[2::]
        forms = '.'.join(forms)
        path = check_path(path)
        nfile = os.path.join(path, "%s.%s" % (name, forms))
        copy_file(file, nfile)
        #move_file(file, nfile)

    return 0


def add_hlep_args(parser):

    parser.add_argument("-f", "--files", nargs='+', metavar='FILE', type=str, required=True,
        help="Input files.")
    parser.add_argument("-b", "--bname", metavar='STR', type=str, required=True,
        help="Input sample and barcode number alignment file.")
    parser.add_argument("-p", "--path", metavar='FILE', type=str, default='./',
        help="The path where the input file needs to be moved.")

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
    rename_lima_file.py  Rename the file name after lima split

barcode_name.txt:
sample  R  F


attention:
    rename_lima_file.py --files * --bname barcode_name.txt
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    rename_lima_file(args.files, args.bname, args.path)


if __name__ == "__main__":

    main()
