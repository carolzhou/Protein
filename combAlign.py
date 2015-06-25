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

HELP_STRING = """Description:  CombAlign parses a set of pairwise structure-based sequence alignments produced by TM-align (see http://zhanglab.ccmb.med.umich.edu/TM-align/), or by DaliLite (see http://www.ebi.ad.uk/Tools/structure/dalilite), and produces a gapped, one-to-many, structure-based sequence alignment (MSSA). Although any alignment program may be used to produce the input pairwise alignments, including sequence- or structure-based tools, the format of the input to combAlign.py should correspond to that produce by TM-align or DaliLite. A reference sequence, corresponding to the sequence or structure that has been compared to multiple others, is used to construct the one-to-many alignment as output. Gaps may be introduced into the reference sequence, as needed, to reflect deletions in the reference with respect to any other sequdnce/structure or residues in the reference that were not aligned with any other sequence/structure according to the criteria aplied in the pairwise alignments. The resulting MSSA may be output in the CombAlign format or in alignedFASTA format.

Type: python combAlign.py input       -for additional information about the required input format
Type: python combAlign.py usage       -for command-line examples for running CombAlign
Type: python combAlign.py in_formats  -for a listing of the acceptable pairwise input formats
Type: python combAlign.py out_formats -for a listing of the acceptable output formats
""" 

INPUT_STRING = """Input requirements:  CombAlign takes as input a single text file containing alignment data. Use 'python combAlign.py' followed by command-line parameters that specify the input file name, input format, output format, and desired line width for output. The input file comprises a single protein sequence in fasta format representing the reference sequence or structure, followed by a set of pairwise alignments. Each pairwise alignment is a comparison of the reference to another protein (in that order). The input file is structured in the following way, regardless of the program used to generated the pairwise alignments:  The fasta reference sequence is preceded by 'REFERENCE' on a single line of text. The word 'REFERENCE' may be followed by a single space and the name of the reference protein. Following the reference fasta is a set of pairwise alignments, each preceded by 'ALIGNMENT' followed by a space and the name of the aligned protein, on a single line of text. The word 'END' on a single line of text signifies the end of the data. A well-formatted input file generated by TM-align will look something like this:

REFRENCE Reston_Ebolavirus_delta_peptide 
>Reston_GPdelta
ELSKEKLATTHPPTTPSWFQRIPLQWFQCSLQDGQRKCRPKV
ALIGNMENT Bundibugyo
EL----------------SKEK----L-ATTHPP-TTP-SWFQRIPLQWFQCSLQDGQRKCRPKV---
                  .::.    : : .::: ::: ::       : : :::    ::.      
--SLPPASPTTKPPRTTKTWFQRIPLQWF-KCETSRGKTQC-------R-P-HPQ----TQS---PQL
ALIGNMENT Sudan
ELSKEKLATTHPPT---------------TPSWF-QRIPLQWFQCSLQ--DGQRKCRPKV---------
                             ::::: :.::. .:: :::  : :    :::         
--------------ELQREESPTGPPGSIRTWFQRIPLGW-FHC-TYQKGK-Q----HCRLRIRQKVEE
ALIGNMENT TaiForest
----E-LSKEKLATTHP-P---TTPSWFQRIPLQ-WFQCSLQDGQRKCRPKV---
    . ::::::: .:: :   ::::       : :  .. ::  :::. .:   
SLLPSPPTTTQPK-TTKNWFQRIPLQ-------WFR--CK-TS--RERT-QCQPQ   
END
"""

USAGE_STRING = """Here are some examples for how to run combAlign.py:

This command will generate an mssa in the CombAlign format with gapped sequences chunked in lengths of 80 characters:
python combAlign.py file=my_align_file1 in_format=DaliLite out_format=standard length=80

This command will generate an mssa in the alignedFASTA format with gapped sequences each in a continuous string 
python combAlign.py file=my_align_file2 in_format=TM-align out_format=aligned_fasta length=0

The input data files must be formatted according to the specified 'format' parameter.
"""

REQUIRED_PARAMS     = 2  # At least 2 parameters must be provided by user
DEFAULT_WIDTH       = 80
MAX_WIDTH           = 250
TM_ALIGN            = "TM-align"
DALI_LITE           = "DaliLite"
COMB_ALIGN          = "combAlign"
ALIGNED_FASTA       = "aligned_fasta"
DEFAULT_FORMAT      = TM_ALIGN 
DEFAULT_OUTPUT_FORMAT = COMB_ALIGN 
CHATTY              = True  # When True, informative statements will be printed
#CHATTY              = False  # When True, informative statements will be printed
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
outFormat = DEFAULT_OUTPUT_FORMAT
ACCEPTABLE_INPUT_FORMATS = (TM_ALIGN, DALI_LITE)
ACCEPTABLE_OUTPUT_FORMATS  = (DEFAULT_OUTPUT_FORMAT, ALIGNED_FASTA)

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
    match = re.search('^in_formats$', sys.argv[i].lower())
    if match:
        print "Acceptable input formats are", ACCEPTABLE_INPUT_FORMATS
        exit(0)
    match = re.search('^out_formats', sys.argv[i].lower())
    if match:
        print "Acceptable output formats are", ACCEPTABLE_OUTPUT_FORMATS
        exit(0)
    match = re.search('=', sys.argv[i])
    if match:
        (parameter,value) = sys.argv[i].split('=')
        if (parameter.lower() == 'in_format' or parameter.lower() == 'input'):
            format = value
            if format not in ACCEPTABLE_INPUT_FORMATS:
                print "Please use an acceptable format:", ACCEPTABLE_INPUT_FORMATS
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
        if (parameter.lower() == 'out_format' or parameter.lower() == 'output'):
            outFormat = value
            if outFormat not in ACCEPTABLE_OUTPUT_FORMATS:
                print "Please request an acceptable output format:", ACCEPTABLE_OUTPUT_FORMATS
                exit(0)

# Reflect parameters
if CHATTY:
    print "Your input file name is", inFile 
    print "Your input format is", format 
    print "Your output format is", outFormat
    print "Your desired line width is", width 

INFILE  = open(inFile,"r")
#OUTFILE = open(outFile,"w")
LOGFILE = open(logFile,"w")
LOGFILE.write("%s%s\n" % ("Name of input file: ", inFile))
LOGFILE.write("%s%s\n" % ("Output mssa is in file: ", mssaFile))
LOGFILE.write("%s%s\n" % ("Input format is: ", format))
LOGFILE.write("%s%s\n" % ("Output format is: ", outFormat))
LOGFILE.write("%s%s\n" % ("Line width is ", width))
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
    print "Setting output format..."
LOGFILE.write("%s\n" % ("Setting output format."))
myAlignment.SetOutputFormat(outFormat)
if CHATTY:
    print "Creating alignment strings..."
LOGFILE.write("%s\n" % ("Creating alignment strings."))
myAlignment.CreateAlignmentStrings()
if CHATTY:
    print "Printing display strings to the output mssa file..."
LOGFILE.write("%s\n" % ("Printing display strings to output mssa file."))
myAlignment.PrintDisplayStrings2file(MSSAFILE,width)

if CHATTY:
    print "Look for your output in file", mssaFile
    print "Done!"
LOGFILE.write("%s\n" % ("Done!"))

# Clean up
INFILE.close()
LOGFILE.close()
