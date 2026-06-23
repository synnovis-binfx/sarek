#!/usr/bin/env python3

__author__ = "alex smith"
__version__ = "1.0.0"
__maintainer__ = "alex smith"
__email__ = "alexander.smith18@nhs.net"
__status__ = "Dev"
__doc__ = '''QC aggregation for nextflow pipeline'''

from uuid import uuid4
import re
import gzip
import os
import inspect
import shutil
import sys
import json
import copy
import time
import glob
from datetime import datetime
from collections import Counter
from pipes import quote

from zipfile import ZipFile
from crimson import fastqc, flagstat, picard, utils

"""Custom metric file parser"""
def parseCustomMetrics(in_data, max_size=1024*1024*1):
    """
    Parses an input Custom Picard-style metrics file into a dictionary.
    Only parses first metrics section!

    :param in_data: Input metrics file.
    :type in_data: str or file handle
    :param max_size: Maximum allowed size of the metrics file (default: 1 MiB).
    :type max_size: int
    :returns: Parsed metrics values.
    :rtype: dict
    """
    # do custom parsing - curenly not using picard parser natively (native picard parser output currently doesnt work well in SQVD)
    with utils.get_handle(in_data) as fh:
        contents = fh.read(max_size)
    sections = contents.strip(os.linesep).split(os.linesep * 2)
    header = [re.compile(r"^#+\s+").sub("", x) for x in sections[0].split(os.linesep)]
    return {
        "header": { "flags": header[1], "time": header[-1] },
        "metrics": picard.parse_metrics(sections[1], "\n"),
        "histogram": None
    }
    
def parseSamtoolsMetrics(in_data, fragments="both"):
    """
    Parses in the input Samtools Cram Stat file to extract % bases over Q30
    
    :param in_data: Input metrics file
    :returns: Parsed metrics values.
    :rtype: dict
    """
    if fragments == "both":
        prefixes = ("FFQ", "LFQ")
    elif fragments in ("FFQ", "LFQ"):
        prefixes = (fragments,)
    else:
        raise ValueError("fragments must be 'both', 'FFQ', or 'LFQ'")

    total = 0
    over_q30 = 0
    with open(in_data) as f:
        for line in f:
            if not line.startswith(prefixes):
                continue
            fields = line.strip().split("\t")
            counts = [int(x) for x in fields[1:]]  # skip cycle number
            for i, count in enumerate(counts):
                total += count
                if i >= 30:  # index 30 = Q30
                    over_q30 += count

    if total == 0:
        raise ValueError("No FFQ/LFQ data found in file.")
    
    pct_q30 = (over_q30 / total) * 100
    
    json_dict = { "type": "samtools cram.stat", "source": in_data, 'data': {"Total Bases": total, "Bases >= Q30": over_q30, "% >= Q30": pct_q30}}

    return json_dict

def parseumihistogram(in_data):
    with open(in_data, 'r') as f:
        lines = f.readlines()
    
    rows = []
    for line in lines[1:]:  # skip header
        values = line.strip().split('\t')
        row = {
            'family_size':  int(values[0]),
            'count':        int(values[1]),
            'fraction':     values[2],
            'fraction_gte': values[3]
        }
        rows.append(row)

    total_count     = sum(r['count'] for r in rows)
    weighted_sum    = sum(r['family_size'] * r['count'] for r in rows)
    avg_family_size = weighted_sum / total_count
    
    json_dict = { "type": "fgio.groupreadsbyumi", "source": in_data, 'data': {"Mean Average": avg_family_size}}
    return json_dict

def aggregateMetrics(input_dir, output_file):
    # parse/get files

    uploadMetrics = []
    # extension for explicit json metrics 'type'
    picard_inserts = ('insert_size_metrics')
    picard_summary = ('alignment_summary_metrics','rna_metrics.txt')
    UMI_depth = ('DUPLICATION','metrics.DUP')
    quality_scores = ('base_quality_metrics')
    custom_extension = ('UMIEXTRACT', 'READS', 'PRIMER', 'ihist.metrics', 'metrics.counts')

    for f in glob.glob(input_dir+'/*'):
        if f.endswith('fastqc.zip'):
            '''FASTQC'''
            if (f):
                with ZipFile(f, 'r') as myzip:
                    for zf in myzip.namelist():
                        if zf.endswith('fastqc_data.txt'):
                            tmpPath = str(uuid4())
                            try:
                                tmpFile = myzip.extract(zf,tmpPath)
                                uploadMetrics.append({
                                    'type': 'fastqc',
                                    'source': f,
                                    'data': dict(fastqc.parse(tmpFile))
                                })
                            except:
                                print ("skip for fastq processing")

        elif f.endswith('flagstat'):
            '''FLAGSTAT'''
            try:
                uploadMetrics.append({
                    'type': 'flagstat',
                    'source': f,
                    'data': flagstat.parse(f)
                })
            except:
                print ("roh ro")

        elif f.endswith('selfSM'):
            '''verifyBamID'''
            try:
                m = MetricsFile(f)
                uploadMetrics.append(m.json(type="verifyBamID"))
            except:
                print ("roh ro")

        elif f.endswith(quality_scores):
            '''Base quality scores'''
            try:
                uploadMetrics.append({
                    'type': 'base_quality_metrics',
                    'source': f,
                    'data': dict(parseCustomMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print ("not formated right for standard picard parsing")

        elif f.endswith(UMI_depth):
            '''UMI Depth'''
            try:
                uploadMetrics.append({
                    'type': 'DUPLICATION',
                    'source': f,
                    'data': dict(parseCustomMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print ("not formated right for standard picard parsing")

        elif f.endswith(picard_summary):
            '''PICARD SUMMARY'''
            try:
                uploadMetrics.append({
                    'type': 'alignment_summary_metrics',
                    'source': f,
                    'data': dict(parseCustomMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print (f, "picard summary formatting incorrect")

        elif f.endswith(picard_inserts):
            '''PICARD INSERTS'''
            try:
                uploadMetrics.append({
                    'type': 'insert_size_metrics',
                    'source': f,
                    'data': dict(parseCustomMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print (f, "not formated right for standard picard parsing")

        elif f.endswith(custom_extension):
            '''CUSTOM'''
            try:
                uploadMetrics.append({
                    'type': 'custom_metrics',
                    'source': f,
                    'data': dict(parseCustomMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print (f, "not formatted right for parsing")
        
        elif f.endswith("cram.stats"):
            '''Samtools Cram Stats'''
            try:
                uploadMetrics.append({
                    'type': 'cram.stats',
                    'source': f,
                    'data': dict(parseSamtoolsMetrics(f))## custom non-standard picard-style metrics parsing
                })
            except:
                print(f, "not formatted right for parsing")
                
        elif f.endswith("umi-grouped_histogram.txt"):
            '''FGBIO Group Reads by UMI'''
            try:
                uploadMetrics.append({
                    'type': 'fgio.groupreadsbyumi',
                    'source': f,
                    'data': dict(parseumihistogram(f))
                })
            except:
                print(f, "not formatted right for parsing")
            
        else:
            print(f+" not parsed")
            

    # write metrics file
    with open(output_file,'w') as jsonMetrics:
        json.dump(uploadMetrics,jsonMetrics,indent=4)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', metavar='FILE', help='output_file_name', type=str, required=True)
    parser.add_argument('--input', metavar='PATH', help='directory of files to parse', type=str, required=True)

    args = parser.parse_args()

    #input_dir = sys.argv [1]
    #output =sys.argv [2]
    aggregateMetrics(args.input, args.output)


if __name__ == "__main__":
    main()
