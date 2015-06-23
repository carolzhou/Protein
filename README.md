# Protein
Code for processing protein data
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
