#############################################################################################
# Module:  CombAlign.py
# Version No.: 1.0
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update:  12 November 2014
#
# Copyright (c) 2015, Lawrence Livermore National Security, LLC. Produced at the Lawrence
# Livermore National Laboratory. Written by Carol L. Ecale Zhou, zhou4@llnl.gov.
# CODE LLNL-CODE-667658 All rights reserved.
# This file is part of CombAlign.
# Please also read README - Our Notice and GNU General Public License.
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License (as published by the Free Software Foundation) version
# 2, dated June 1991.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the terms and conditions of the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA.
#
# Description:  This code parses a set of pairwise structure-based sequence alignments 
# produced by TMalign (see http://zhanglab.ccmb.med.umich.edu/TM-align/), and produces a
# gapped, one-to-many, structure-based sequence alignment.  Although any alignment
# program may be used to produce the pairwise alignments, the format of the input to
# combAlign.py should correspond to that produced by TMalign. A reference
# sequence, corresponding to the structure that has been compared to multiple other
# structures, is used to construct the one-to-many alignment as output. Gaps may be
# introduced into the reference sequence, as needed, to reflect deletions in the
# reference with respect to any other structure or residues in the reference that were 
# not alignmed with any other structure according to the criteria applied in the
# pairwise structure alignments. 
#
# Input requirements:  
# 1) CombAlign.py takes as input a single text file, followed by a positive integer.
# 2) The input file to combAlign.py comprises a single protein
# sequence in fasta format representing the reference structure, followed by a set of
# pairwise alignments. In each pairwise alignment, the first sequence listed is from the
# reference structure. The sequence should be contained on a single line of text, as
# should the data lines for each pairwise alignment. The input file should contain no
# blank lines, but should end in a single carriage return.
# 3) The second input parameter indicates the desired output alignment segment length
# per line of output text. This parameter should be a positive number less than 250. The
# default alignment segment length is 80.   
#
#############################################################################################

import sys
import string
import re, os
import alignment              # alignment.py module by C. Zhou

p_number = re.compile('\d+')  # search pattern for number string

# FILES

inFile   = ""  # user provided
outFile  = "./combAlign.out"
mssaFile = "./combAlign.mssa"

# HELP STRINGS and CONSTANTS

HELP_STRING = """Description:  This code parses a set of pairwise structure-based sequence alignments 
produced by TMalign (see http://zhanglab.ccmb.med.umich.edu/TM-align/), and produces a
gapped, one-to-many, structure-based sequence alignment.  Although any alignment
program may be used to produce the pairwise alignments, the format of the input to
combAlign.py should correspond to that produced by TMalign. A reference
sequence, corresponding to the structure that has been compared to multiple other
structures, is used to construct the one-to-many alignment as output. Gaps may be
introduced into the reference sequence, as needed, to reflect deletions in the
reference with respect to any other structure or residues in the reference that were 
not alignmed with any other structure according to the criteria applied in the
pairwise structure alignments.
""" 

INPUT_STRING = """Input requirements:  
1) CombAlign.py takes as input a single text file, followed by a positive integer.
2) The input file to combAlign.py comprises a single protein
sequence in fasta format representing the reference structure, followed by a set of
pairwise alignments. In each pairwise alignment, the first sequence listed is from the
reference structure. The sequence should be contained on a single line of text, as
should the data lines for each pairwise alignment. The input file should contain no
blank lines, but should end with a single carriage return.
3) The second input parameter indicates the desired output alignment segment length
per line of output text. This parameter should be a positive number less than 250. The
default alignment segment length is 80, and the parameter is optional.
"""   

USAGE_STRING = "python combAlign.py reference_alignment_file (optional)lineLength"

DEFAULT_WIDTH = 80
MAX_WIDTH     = 250

# DATA STRUCTURES
refSeq = {            # Data pertaining to the sequence of the reference structure
    "reference" : "",
    "header"    : "",
    "sequence"  : "",
    }
reference = ""  # Name of the reference structure in the TMalign alignments
                # reference is the header of the fasta provided at top of input alignments file 

# VARIABLES
width = DEFAULT_WIDTH  # width of output alignment segments

# Get input parameter(s) 
# User needs to provide name of file containing reference fasta sequence and list of pairwise alignmentsa
# Optional is the length of the alignment segments

ACCEPTABLE_ARG_COUNT = (2,3) # ie, code + input file or help word + an optional positive integer
argCount = len(sys.argv)
if argCount in ACCEPTABLE_ARG_COUNT:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print HELP_STRING
        print USAGE_STRING
        exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print INPUT_STRING
        exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print USAGE_STRING
        exit(0)
    inFile    = sys.argv[1] 
    if argCount == 3:
        width = sys.argv[2]
        match = re.search(p_number, width)
        if not match or (width < 0 and width > MAX_WIDTH):
            print "Unexpected input parameter", width
            print USAGE_STRING
            exit(0)
else:
    print "Type: python combAlign.py help"

INFILE  = open(inFile,"r")
OUTFILE = open(outFile,"w")
OUTFILE.write("%s%s\n" % ("Name of input file: ", inFile))
OUTFILE.write("%s%s\n" % ("Output mssa is in file: ", mssaFile))
MSSAFILE = open(mssaFile,"w")

# Read reference fasta sequence and create base data structure

fLines = INFILE.read().splitlines()
lineCount = len(fLines)
refSeq["header"]    = fLines[0] 
refSeq["sequence"]  = fLines[1] 
refSeq["reference"] = refSeq["header"].lstrip('>')

OUTFILE.write("%s%s\n" % ("Infile is ", inFile))
OUTFILE.write("%s%s\n" % ("Infile lineCount is: ", lineCount))
OUTFILE.write("%s%s\n" % ("Reference fasta is: ", refSeq["reference"]))
OUTFILE.write("%s\n%s\n" % (refSeq["header"], refSeq["sequence"]))
OUTFILE.write("\n%s\n" % ("Structure-based multiple sequence alignment:"))

alignment = alignment.Alignment() #***
alignment.EnterReference(refSeq)

pairwise = {   # dict holding next pairwise alignment 
    "matchName"          : "",  # name of sequence aligned to reference
    "referenceLine"      : "",  # reference sequence (complete) with possible gaps as '-'
    "correspondenceLine" : "",  # string of characters: blank, '.', or ':'
    "matchLine"          : "",  # match sequence (maybe incomplete) with possible gaps as '-'
    }
failure = 0  # test result from method call

# Insert R-R correspondences for each aligned structure

for i in xrange (2,lineCount-1):
    if i%4 == 0:
        pairwise["correspondenceLine"] = fLines[i]
    elif i%4 == 1:
        pairwise["matchLine"] = fLines[i]
        failure = alignment.AddAlignment(pairwise)  # add next pairwise to alignment object
        if failure == 0:
            print "Pairwise alignment ", pairwise["matchName"], " successfully added."
        else:
            print "Method AddAlignment Failure code ", failure
    elif i%4 == 2:
        pairwise["matchName"] = fLines[i]
    elif i%4 == 3:
        pairwise["referenceLine"] = fLines[i]
        
# Print multiple structure-based sequence alignment

alignment.CreateAlignmentStrings()
alignment.PrintDisplayStrings2file(MSSAFILE,DEFAULT_WIDTH)

# Clean up

INFILE.close()
OUTFILE.close()

