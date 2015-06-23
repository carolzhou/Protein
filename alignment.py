##################################################################################################
# Module:  alignment.py
# Version No.: 1.1
#
# Programmer: Carol L. Ecale Zhou
#
# Most recent update: 22 June 2015
#
# Description: This module contains classes and methods for representing sequence alignment.
# Class: alignment
#
# This code constructs a multiple sequence alignment allowing for gaps in the reference sequence.
# Class Alignment is intended to be used as follows:
#   1) Input a reference sequence in fasta format using method EnterReference
#   2) Input a set of pairwise alignments, one at a time, using method AddAlignment
#   3) Create a multiple alignment using method CreateAlignmentStrings
#   4) Print the multiple alignment using method PrintDisplayStrings
#      The default is to print the multiple alignment horizontally as single strings,
#      however, you may specify the number of desired positions per line of text.
#
# Programmer's notes:
#   a) A terminal '*' is added to the end of the reference sequence and to the ends of
#      the alignment strings, for convenience.
#
# Copyright 2014 by Carol L. Ecale Zhou, Lawrence Livermore National Security.  All Rights Reserved. 
# Permission to use, copy, modify, and distribute this software and its documentation for educational, 
# research, and not-for-profit purposes, without fee and without a signed licensing agreement, is 
# hereby granted, provided that the above copyright notice, this paragraph and the following two 
# paragraphs appear in all copies, modifications, and distributions. Contact Office of XXXX,
# Lawrence Livermore National Security (LLNS) for commercial licensing opportunities.
#    In no event shall LLNS be liable to any party for direct, indirect, special, incidental, or
# consequential damages, including lost profits, arising out of the use of this software and its
# documentation, even if LLNS has been advised of the possibility of such damage.
#    LLNL disclaims any warranties, including, but not limited to, the implied warranties of
# merchantability and fitness for a particular purpose. The software and accompanying documentation
# if any, provided hereunder is provided "as is", LLNS has no obligation to provide maintenance,
# support, updates, enhancements, or modifications.
##################################################################################################

import re, copy
p_comment = re.compile('^#')

TM_ALIGN  = "TM-align"
DALI_LITE = "DaliLite"

class Alignment(object):

    def __init__(self,format):
        self.molecule    = "protein"   # default; 'nucleic acid' or 'protein'
        self.structure   = True        # default; sequence-based if False
        self.reference   = "unknown"   # name of reference structure|sequence
        self.refHeader   = ">undefined"
        self.refSequence = "empty"
        self.alignmentCount = 0        # increments with each added alignment
        self.multiAlignment = []       # list of pairwiseData objects, one per position along reference
        self.pair_i      = 0           # index for self.alignment list (max pair_i is No. of pairwise alignments-1)
        self.matchNameList = []        # captures list of match sequence names
        self.pairwiseData = {          # Each residue of reference is tagged (via array position) with data:
            "refChar"            : "", # Residue in the reference sequence at current position (should mimic self.refSequence)
            "correspondenceList" : [], # list of values for R's that correspond ('.', ':', or ' ')
            "matchList"          : [], # list of R's that correspond to position on reference
            "gapList"            : [], # list of strings, each comrpising a loop in match, gap in reference
            }
        self.loopStrings    = []       # list of strings in match sequences that extend gaps in reference
        self.refDisplayString = ""     # holds final (gapped) reference sequence
        self.displayStrings = []       # list of stringPairs for formatted alignment display
        self.stringPair = {            # a matching sequence and the correspondence values ('.', ':', or ' ')
            "correspondence"  : "",
            "match"           : "",
            }
        if (format == TM_ALIGN):
            self.method = TM_ALIGN
        elif (format == DALI_LITE):
            self.method = DALI_LITE
        else:
            self.method = TM_ALIGN 

    def EnterReference(self,refSeq):   # Enters a reference sequence 
        # Method EnterReference() establishes the data structure for capturing the residue-
        # by-residue correspondences between each reference residue and its counterpart in
        # each aligned sequence or structure. Thus, for each position in the reference
        # sequence, a pairwiseData object is created and assigned to that position,
        # including the terminal '*'. 

        if isinstance(refSeq,dict):
            if "reference" in refSeq:
                self.reference  = refSeq["reference"]
            else:
                return False
            if "header" in refSeq:
                self.refHeader  = refSeq["header"]
            else:
                return False
            if "sequence" in refSeq:
                self.refSequence = refSeq["sequence"]
                if self.refSequence[-1] != '*':
                    self.refSequence += '*'  # add terminal '*'
            else:
                return False

            # Construct self.multiAlignment list of pairwiseData objects, one per reference residue
            # Data details are entered in method AddAlignment
            refLength = len(self.refSequence) # sequence length plus 1 (including the '*')
            for i in xrange(0,refLength):  
                nextPairwise = copy.deepcopy(self.pairwiseData)
                self.multiAlignment.append(nextPairwise)
            return True
        else:
            return False

    def AddAlignment(self,newAlignment):  # Add a new alignment (stringPair) to self.displayStrings
        # Fields reference, correspondence, and match are lines of text
        # containing the alignment strings for a pairwise alignment.
        # Example:  MGPKAKAEA--SKPHQIPQIPVKLPFVTAPDAL  # reference
        #                      ..         .   .  ::::  # correspondence
        #           ---------TDPA---------P---P--PTAL  # match
        # For the multiple alignment, gaps are introduced into the reference sequence string.
        # For example, the 2-residue gap is recorded as "TD" in the gapList field of self.
        # pairwiseData at array position 8 (self.multiAlignment[8]), attached to preceding 'A',
        # in self.pairwiseData["gapList"].
        
        if isinstance(newAlignment,dict):
            if "matchName" in newAlignment:
                matchName = newAlignment["matchName"]
            else:
                return 2
            if "referenceLine" in newAlignment:
                reference = newAlignment["referenceLine"]
                if reference[-1] != '*':
                    reference += '*'
                referenceLen = len(reference)
            else:
                return 2
            if "correspondenceLine" in newAlignment:
                correspondence = newAlignment["correspondenceLine"]
                if correspondence[-1] != '*':
                    correspondence += '*'
                correspondenceLen = len(correspondence)
            else:
                return 2
            if "matchLine" in newAlignment:
                match = newAlignment["matchLine"]
                if match[-1] != '*':
                    match += '*'
                matchLen = len(match)
            else:
                return 2

            # Input alignment should consiste of 3 strings (reference, correspondence, and
            # aligned sequence), which are of equal length
            if referenceLen == correspondenceLen and referenceLen == matchLen and matchLen > 0:
                self.matchNameList.append(matchName)
                gapString = ""
                ref_i = 0  # index of positions along reference sequence (no gaps)
                # i is index of positions along gapped reference from alignment
                for i in xrange(0,matchLen):
                    if reference[i] == '-' or reference[i] == '.':  # gap was introduced into reference sequence
                        gapString += match[i]
                    else:  # next character in reference line is not gap character
                        self.multiAlignment[ref_i]["refChar"] = self.refSequence[ref_i]  #*** redundant: assigns every time
                        self.multiAlignment[ref_i]["matchList"].append(match[i])
                        self.multiAlignment[ref_i]["correspondenceList"].append(correspondence[i])
                        self.multiAlignment[ref_i-1]["gapList"].append(gapString) # associate gapstring w/previous ref R
                        gapString = ""  # reset
                        ref_i += 1
                # Construct empty display strings for correspondence and match\
                self.loopStrings.append("")  # construct/add to list of loopStrings, one for each match sequence
                newStringPair = copy.deepcopy(self.stringPair)
                self.displayStrings.append(newStringPair)
                self.alignmentCount += 1
                return 0
            else:
		print "referenceLen is", referenceLen, "correspondenceLen is", correspondenceLen, "matchLen is", matchLen
                return 3  # error code
        else:
            return 1  # error code

    def IsLoop(self):
        result = False
        for j in xrange(0,self.alignmentCount):
            if self.loopStrings[j]:
                result = True
        return result

    def CreateAlignmentStrings(self):
        for i in xrange(0,len(self.refSequence)): # Add next chars to display strings; iterate through multiAlignment data structure
            for j in xrange(0,self.alignmentCount): # Append gapList at prev refSeq pos to j's current loopString
                self.loopStrings[j] += self.multiAlignment[i-1]["gapList"][j] # recall: gap is associated w/prev pos
            while self.IsLoop():
                self.refDisplayString += '-'  # open/continue gap: reference seq's display string gets gap character
                for j in xrange(0,self.alignmentCount):  
                    if self.loopStrings[j]:  # true if string not empty
                        self.displayStrings[j]["match"] += self.loopStrings[j][0]  # add 1st char from loopString[j]
                        self.loopStrings[j] = self.loopStrings[j][1:]  # remove 1st char from loopString[j]
                        self.displayStrings[j]["correspondence"] += ' '
                    else:
                        self.displayStrings[j]["match"]          += '-'
                        self.displayStrings[j]["correspondence"] += ' '
            self.refDisplayString += self.multiAlignment[i]["refChar"]
            for j in xrange(0,self.alignmentCount):
                self.displayStrings[j]["match"] += self.multiAlignment[i]["matchList"][j]
                self.displayStrings[j]["correspondence"] += self.multiAlignment[i]["correspondenceList"][j]

    def PrintReference(self):
        print "REFERENCE SEQUENCE:"
        print self.refHeader
        print self.refSequence

    def PrintPairwiseData(self):
        position = 1
        print "LIST OF PAIRWISE DATA VALUES:"
        for pairwise in self.multiAlignment:
            print "Position number ", position
            print pairwise["refChar"]
            print pairwise["matchList"]
            print pairwise["correspondenceList"]
            print pairwise["gapList"]
            position += 1

    def PrintDisplayStrings2file(self, OUTFILE, width=0):  # default is to print all sequence lines as single string
        width = int(width) # cast to integer

        if width == 0:  # Simple case:  print as is
            OUTFILE.write("%s%s%s\n" % (self.refDisplayString, " ", self.refHeader))
            for stringPair in self.displayStrings:
                OUTFILE.write("%s\n" % (stringPair["correspondence"]))
                OUTFILE.write("%s\n" % (stringPair["match"]))
                
        else: # Split each reference and correspondence/match
            # We need to determine how many chunks to split the alignment data
            # into, depending on how long each chunk should be. 
            list_i = 0 # need this to access match names in self.matchNameList
            segmentSets = []  # list of stringPair segments
            segmentSetDict = {
                "reference"      : "",
                "correspondence" : "",
                "match"          : "",
                "matchName"      : "",
                }
            segmentCount = 0
            # Split reference, correspondence, and match strings
            # Assuming you have pairwise alignments A, B, C, and D,
            # comprising alignment chunks A1, A2, A3, etc.,
            # This places segmentSet objects as follows in segmentSets list:
            #    A1  A2  A3  B1  B2  B3  C1  C2  C3  D1  D2  D3
            # To be printed as:
            #    A1
            #    B1
            #    C1
            #    D1
            #    A2
            #    B2  etc.
            OUTFILE.write("%s%s\n" % ("Input format was ", self.method))
            OUTFILE.write("%s%s\n" % ("Length of gapped reference: ", len(self.refDisplayString))) 
            if len(self.refDisplayString)%width == 0:
                segmentCount = len(self.refDisplayString)/width
            else:
                segmentCount = len(self.refDisplayString)/width + 1
            OUTFILE.write("%s%s\n" % ("There will be this many segments: ", segmentCount))
            OUTFILE.write("%s%s\n" % ("The reference structure was:  ",self.reference))
            OUTFILE.write("%s%s\n" % ("The reference fasta was:  ",self.refHeader[1:]))
            OUTFILE.write("%s\n" % ("The compared structures were:"))
            for name in self.matchNameList:
                OUTFILE.write("%s%s\n" % ("  ",name))

            # Perform splits: split data into designated chunks
            last = False
            for stringPair in self.displayStrings:  # iterate through alignment strings
                for i in xrange(0,segmentCount):  # len is same for all strings; split segmentCount times
                    segmentSet = copy.deepcopy(segmentSetDict) # create new segmentSet
                    if i == segmentCount-1:  # last segment (may be shorter than width)
                        refSegment   = self.refDisplayString       [i*width:]
                        correSegment = stringPair["correspondence"][i*width:]
                        matchSegment = stringPair["match"]         [i*width:]
                        segmentSet["correspondence"] = correSegment
                        segmentSet["match"]          = matchSegment
                        segmentSet["reference"]      = refSegment
                        segmentSet["matchName"]      = self.matchNameList[list_i] # Name already appended in CreateAlignmentStrings
                    else:
                        refSegment   = self.refDisplayString       [i*width:(i+1)*width]
                        correSegment = stringPair["correspondence"][i*width:(i+1)*width]
                        matchSegment = stringPair["match"]         [i*width:(i+1)*width]
                        segmentSet["correspondence"] = correSegment
                        segmentSet["match"]          = matchSegment
                        segmentSet["reference"]      = refSegment
                        segmentSet["matchName"]      = self.matchNameList[list_i] # need to re-insert name due to line split
                    segmentSets.append(segmentSet)
                list_i += 1
            OUTFILE.write("%s%s\n" % ("There are this many segment sets: ", len(segmentSets)))

            refHeader = self.refHeader[1:]  # trim '>' from header

            for i in xrange(0,segmentCount):
                first = True
                for j in xrange(0,self.alignmentCount):
                    k = (j * segmentCount) + i
                    if first:
                        OUTFILE.write("%s%s%s\n" % (segmentSets[k]["reference"]," ",refHeader)) 
                        first = False
                    OUTFILE.write("%s\n"     % (segmentSets[k]["correspondence"]))
                    OUTFILE.write("%s%s%s\n" % (segmentSets[k]["match"]," ", segmentSets[k]["matchName"]))
                OUTFILE.write("%s\n" % (""))

    def PrintAll(self):
        print "ALL DATA:"
        print "Molecule type: ", self.molecule
        if self.structure:
            print "Structure-based multiple sequence alignment"
        else:
            print "Sequence-based alignment"
        print "Method used: ", self.method
        print "Sequences aligned to reference:"
        for name in self.matchNameList:
            print "   ", name
        self.PrintReference()
        self.PrintPairwiseData()
        self.PrintDisplayStrings()
                
