#############################################################################################
# Module:  combAlign.py
# Version No.: 1.1
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update:  22 June 2015
#
# Description:  This code parses a set of pairwise structure- or sequence-based alignments 
# produced by TMalign (see http://zhanglab.ccmb.med.umich.edu/TM-align/), or Dali Lite
# (see http://www.ebi.ac.uk/Tools/struture/dalilite/), and produces a gapped,
# one-to-many, sequence alignment.  Although any alignment program may be used to 
# produce the pairwise alignments, the format of the input to combAlign.py
# should correspond to that produced by TMalign or Dali Lite. A reference fasta sequence
# corresponding to the sequence or structure that has been compared to multiple other
# sequences/structures, is used to construct the one-to-many alignment as output. Gaps may be
# introduced into the reference sequence, as needed, to reflect deletions in the
# reference with respect to any other structure or residues in the reference that were 
# not alignmed with any other sequence/structure according to the criteria applied in the
# pairwise sequence/structure alignments. 
#
# Input requirements:  
# 1) CombAlign.py takes as input a single text file, followed by the format (TM-align
# or DaliLite), plus an optional positive integer.
# 2) The input file to combAlign.py comprises a single protein sequence in fasta format,
# representing the reference sequence/structure, followed by a set of pairwise
# alignments. In each pairwise alignment, the first sequence listed is from the
# reference. 
# 3) The second input parameter indicates the format used to perform the alignments.
# Currently the acceptable formats are TM-align and DaliLite.
# 4) The third (optional) input parameter indicates the desired output alignment segment
# length per line of output text. This parameter should be a positive number less than 250.
# The default alignment segment length will be set to 80, if not provided by the user.   
#
# Programmer's notes:
# This code can be easily modified to accommodate additional input file formats. Any
# input file format should be parsed into the 'pairwise' data structure. If this is done
# properly, then all subsequent calls to class alignment.py will function in exactly
# the same way. To add additional formats, one need only add code under the 'ALIGN'
# block (see line 332, beginning 'if ALIGN'). 
#
# Copyright 2014 by Carol L. Ecale Zhou, Lawrence Livermore National Security. All Rights
# Reserved. Permission to use, copy, modify, and distribute this software and its 
# documentation for educational, research, and not-for-profit purposes, without fee and
# without a signed licensing agreement, is hereby granted, provided that the above
# copyright notice, this paragraph and the following two paragraphs appear in all copies,
# modifications, and distributions. Contact Office of XXXX, Lawrence Livermore National
# Security for commercial licensing opportunities.
#    In no event shall LLNS be liable to any party for direct, indirect, special, incidental,
# or consequential damages, including lost profits, arising out of the use of this
# software and its documentation, even if LLNS has been advised of the possibility of
# such damage.
#    LLNS disclaims any warranties, including, but no limited to, the implied 
# warranties of merchantability and fitness for a particular purpose. The software and
# accompanying documentation, if any, provided hereunder is provided "as is", LLNS has
# no obligation to provide maintenance, support, updates, enhancements, or modifications.
#############################################################################################

import sys
import string
import re, os
import alignment              # alignment.py module by C. Zhou

p_number = re.compile('\d+')  # search pattern for number string

# FILES

inFile   = ""  # user provided
outFile  = "./combAlign.out"
logFile  = "./combAlign.log"
mssaFile = "./combAlign.mssa"

# HELP STRINGS and CONSTANTS

HELP_STRING = """Description:  This code parses a set of pairwise structure-based sequence alignments 
produced by TMalign (see http://zhanglab.ccmb.med.umich.edu/TM-align/), or by 
Dali Lite (see http://www.ebi.ad.uk/Tools/structure/dalilite), and produces a 
gapped, one-to-many, structure-based sequence alignment.  Although any alignment 
program may be used to produce the pairwise alignments, including sequence- 
or structure-based tools, the format of the input to 
combAlign.py should correspond to that produced by TMalign or Dali Lite. A reference
sequence, corresponding to the sequence or structure that has been compared to multiple
others, is used to construct the one-to-many alignment as output. Gaps may be
introduced into the reference sequence, as needed, to reflect deletions in the
reference with respect to any other sequence/structure or residues in the reference that were 
not aligned with any other sequence/structure according to the criteria applied in the 
pairwise alignments. Type: python combAlign.py input, for details about 
the expected input file formats.
""" 

INPUT_STRING = """Input requirements:  

CombAlign.py takes as input a single text file, followed by the format used 
to align the sequences or structures, with an optional third parameter indicating
the desired output segment length. The name of the text file is prefixed by 'file=',
the format is prefixed by 'format=', and the desired line segment length is
prefixed by 'width='. For example: python combAlign.py file=myFile format=TM-align width=80

The input file to combAlign.py comprises a single protein
sequence in fasta format representing the reference sequence or structure, followed by a set of
pairwise alignments. In each pairwise alignment, the first sequence listed is from the
reference structure. Each portion of input data is introduced by a string that signals
the type of data to follow:  the reference fasta is preceded by 'REFERENCE myReference',
where myReference is an optional name of the reference. If a name is not provided, the
fasta header will be substituted as the reference name. Each alignment is introduced by
'ALIGNMENT myNextProtein', where myNextProtein identifies the sequence or structure to
which the reference was aligned. Following the last alignment is the termination signal,
'END', which comprises the last line of the input data file.  For examples of 
alignment data files comprising TM-align or Dali Lite data formats, see files,
sGP_Reston_TM-align_In.txt and GP_Reston_DaliLite_in.txt

The second input parameter indicates the format: either 'TM-align' or 'DaliLite'.

The third input parameter indicates the desired output alignment segment length
per line of output text. This parameter should be a positive number less than 250. The
default alignment segment length is 80, and the parameter is optional.
"""   

USAGE_STRING = """Here is an example for how to run combAlign.py:
python combAlign.py file=my_align_file format=acceptable_format length=80"
or, for information, type \n\'python combAlign.py\' followed by 'help', 'input', 'format'
"""

REQUIRED_PARAMS     = 2  # At least 2 parameters must be provided by user
DEFAULT_WIDTH       = 80
MAX_WIDTH           = 250
TM_ALIGN            = "TM-align"
DALI_LITE           = "DaliLite"
DEFAULT_FORMAT      = TM_ALIGN 
CHATTY              = True  # When True, informative statements will be printed
DL_START            = 6     # column number at which Dali Lite alignment data begins
DL_END              = 66    # col num at which DL alignments end, except final fragment
REFERENCE_LINE      = True  # Bool; controls capture of data line for TM-align format
CORRESPONDENCE_LINE = False #  "
MATCH_LINE          = False #  "

# DATA STRUCTURES
refSeq = {            # Data pertaining to the sequence of the reference structure
    "reference" : "",
    "header"    : "",
    "sequence"  : "",
    }
reference = ""  # Name of the reference structure in the TMalign alignments
                # reference is the header of the fasta provided at top of input alignments file 

# VARIABLES
width  = DEFAULT_WIDTH   # width of output alignment segments
format = DEFAULT_FORMAT  # program used to perform pairwise alignments, begin w/default 
ACCEPTABLE_FORMATS = (TM_ALIGN, DALI_LITE)

# PATTERNS
p_dataFragment = re.compile('[\w\.\-]*')
p_refName      = re.compile('^REFERENCE\s+(\w+)')

# Get input parameter(s) 
# User needs to provide 1) name of file containing reference fasta sequence and 
# list of pairwise alignments, 2) the format of (ie, program used to generate) the 
# alignment, 3) (Optional) desired line length for output.

argCount = len(sys.argv)
if (argCount < REQUIRED_PARAMS):
    print "Insufficient parameters. Type: \'python combAlign.py help\'" 
    exit(0)

for i in range(1,argCount):
    match = re.search("^help$", sys.argv[i].lower())
    if match:
        print HELP_STRING
        print USAGE_STRING
        exit(0)
    match = re.search('^input$', sys.argv[i].lower())
    if match:
        print INPUT_STRING
        exit(0)
    match = re.search('^usage$', sys.argv[i].lower())
    if match:
        print USAGE_STRING
        exit(0)
    match = re.search('^format$', sys.argv[i].lower())
    if match:
        print "Accepatble formats are", ACCEPTABLE_FORMATS
        exit(0)
    match = re.search('=', sys.argv[i])
    if match:
        (parameter,value) = sys.argv[i].split('=')
        if (parameter.lower() == 'format'):
            format = value
            if format not in ACCEPTABLE_FORMATS:
                print "Please use an acceptable format:", ACCEPTABLE_FORMATS
                exit(0)
        if (parameter.lower() == 'file'): 
            inFile = value
        if (parameter.lower() == 'width' or parameter.lower() == 'length'):
            width = value
            match = re.search('[^\d]', width)
            if match:
                print "Choose a more realistic line width."
                print USAGE_STRING
                exit(0)
            if (int(width) < 0 or int(width) > 255):
                print "Choose a more realistic line width."
                print USAGE_STRING
                exit(0)

# Reflect parameters
if CHATTY:
    print "Your input file name is", inFile 
    print "Your input format is", format 
    print "Your desired line width is", width 

INFILE  = open(inFile,"r")
#OUTFILE = open(outFile,"w")
LOGFILE = open(logFile,"w")
LOGFILE.write("%s%s\n" % ("Name of input file: ", inFile))
LOGFILE.write("%s%s\n" % ("Output mssa is in file: ", mssaFile))
MSSAFILE = open(mssaFile,"w")

# Set up data structure for holding alignments
pairwise = {   # dict holding next pairwise alignment 
    "matchName"          : "",  # name of sequence aligned to reference
    "referenceLine"      : "",  # reference sequence (complete) with possible gaps as '-'
    "correspondenceLine" : "",  # string of characters: blank, '.', or ':'
    "matchLine"          : "",  # match sequence (maybe incomplete) with possible gaps as '-'
    }
failure = 0  # will be test result from method call

# Create an Alignment object
myAlignment = alignment.Alignment(format)

if CHATTY:
    print "Acquiring R-R correspondences from input file..."

# Switches that control handling of input file data lines
FASTA = False
ALIGN = False

fLines = INFILE.read().splitlines()
lineCount = len(fLines)

for i in xrange(0,lineCount):
    nextLine = fLines[i]  # Get next data line
    if (nextLine == '' or nextLine == '\^#'):  # comment line starts with '#' 
        continue  # skip blank and comment lines

    # Check if end of data input is reached
    match = re.search('^END', nextLine)
    if match:  
        if ALIGN:   # (should be true) save away current alignment data
            failure = myAlignment.AddAlignment(pairwise)
            if CHATTY:
                if failure == 0:
                    print "Pairwise alignment", pairwise["matchName"], "successfully added."
                else:
                    print "Method AddAlignment failure code", failure, "at", pairwise["matchName"]
            # Reset
            pairwise["referenceLine"] = ''
            pairwise["correspondenceLine"] = ''
            pairwise["matchLine"] = ''
        else:
            if CHATTY:
                print "WARNING:  Problem with input data file at line", i 
            LOGFILE.write("\n%s%s\n" % ("Problem with input data file at line", i))
        break  # jump out of loop, although this should be last iteration anyway

    # Check if alignment is next
    match = re.search('^ALIGNMENT', nextLine)
    if match:

        # If ALIGN flag is 'on', then last data item was the previous alignment
        if ALIGN:  # First, wrap up previous alignment
            failure = myAlignment.AddAlignment(pairwise)
            if CHATTY:
                if failure == 0:
                    print "Pairwise alignment", pairwise["matchName"], "successfully added."
                else:
                    print "Method AddAlignment failure code", failure, "at", pairwise["matchName"]
            # Reset
            pairwise["referenceLine"] = ''
            pairwise["correspondenceLine"] = ''
            pairwise["matchLine"] = ''

            # Capture name of the next aligned sequence/structure from tag 
            (preamble, pairwise["matchName"]) = nextLine.split(' ')

        # If FASTA flag still 'on', then last data item was the reference fasta
        elif FASTA:  # Register the reference fasta before processing alignment 
            myAlignment.EnterReference(refSeq)
            if CHATTY:
                print "Your reference fasta sequence is:"
                print refSeq["header"]
                print refSeq["sequence"] 
            FASTA = False  # done!, and there should be only one fasta in input data file

            # Next order of business, capture the name of the aligned sequence/structure
            ALIGN = True
            (preamble, pairwise["matchName"]) = nextLine.split(' ')
        else:
            if CHATTY:
                print "WARNING: Problem with input file at line", i 
            LOGFILE.write("\n%s%s\n" % ("WARNING:  Problem with input file at line", i))

        continue

    # Check if fasta header is next
    match = re.search('^REFERENCE', nextLine) # Check for reference fasta 
    if match:  # Reference fasta is next
        FASTA = True
        ALIGN = False
        match = re.search('^REFERENCE\s+(\w+)', nextLine)
        if match:
            refSeq["reference"] = match.group(1)
        else:
            refSeq["reference"] = '' 
        continue
    
    if FASTA:  # Capture fasta data

        # Check if it's a header line
        match = re.search('^\>', nextLine)
        if match:
            refSeq["header"] = nextLine
            if (refSeq["reference"] == ''):  # use header as reference name if user did not provide
                refSeq["reference"] = refSeq["header"].lstrip('>')
            continue
 
        # Check if it's a sequence line  
        match = re.search('[\w\*]', nextLine)  # should be letters or '*'
        if match:
            refSeq["sequence"] = refSeq["sequence"] + nextLine
            continue

    if ALIGN:  # Capture alignment data. Additional formats could be accommodated here.

        if (format == TM_ALIGN):
            if REFERENCE_LINE:
                pairwise["referenceLine"] = nextLine
                REFERENCE_LINE = False
                CORRESPONDENCE_LINE = True
            elif CORRESPONDENCE_LINE:
                pairwise["correspondenceLine"] = nextLine
                CORRESPONDENCE_LINE = False
                MATCH_LINE = True
            elif MATCH_LINE:
                pairwise["matchLine"] = nextLine
                MATCH_LINE = False
                REFERENCE_LINE = True
            else:
                LOGFILE.write("%s%s\n" % ("Problem at line", i)) 
            continue

        if (format == DALI_LITE):
            match = re.match('^DSSP', nextLine)
            if match:  # Calculate begin/end column numbers for extracting data fragment
                dataSet = re.findall(p_dataFragment, nextLine) # dataSet[2] contains DSSP result string
                DL_END = len(dataSet[3]) + 6 # calculate column number at which data fragment ends
            match = re.match('^Query', nextLine)
            if match:
                pairwise["referenceLine"] += nextLine[DL_START:DL_END]  # fragment of nextLine containing sequence
            match = re.search('^ident', nextLine)
            if match:
                pairwise["correspondenceLine"] += nextLine[DL_START:DL_END] 
            match = re.search('^Sbjct', nextLine)
            if match:
                pairwise["matchLine"] += nextLine[DL_START:DL_END] 
            continue

    if FASTA:  # Capture fasta data
        match = re.search('^\>', nextLine)
        if match:
            refSeq["header"] = nextLine
            if (refSeq["reference"] == ''):  # use header as reference name if user did not provide
                refSeq["reference"] = refSeq["header"].lstrip('>')
            continue

        match = re.search('[\w\*]', nextLine)  # should be letters or '*'
        if match:
            refSeq["sequence"] += nextLine
            continue

LOGFILE.write("%s%s\n" % ("Infile is ", inFile))
LOGFILE.write("%s%s\n" % ("The format is ", format))
LOGFILE.write("%s%s\n" % ("Infile lineCount is: ", lineCount))
LOGFILE.write("%s%s\n" % ("Reference fasta is: ", refSeq["reference"]))
LOGFILE.write("%s\n%s\n" % (refSeq["header"], refSeq["sequence"]))
LOGFILE.write("\n%s\n" % ("Calculating combined alignment."))

# Print multiple structure-based sequence alignment
if CHATTY:
    print "Creating alignment strings..."
LOGFILE.write("%s\n" % ("Creating alignment strings."))
myAlignment.CreateAlignmentStrings()
if CHATTY:
    print "Printing display strings to the output mssa file..."
LOGFILE.write("%s\n" % ("Printing display strings to output mssa file."))
myAlignment.PrintDisplayStrings2file(MSSAFILE,width)

if CHATTY:
    print "Done!"
LOGFILE.write("%s\n" % ("Done!"))

# Clean up
INFILE.close()
LOGFILE.close()
