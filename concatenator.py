#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  untitled.py
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  



    def read_alignments (self, input_list, alignment_format, progress_stat=True):
        """ Function that parses multiple alignment/loci files and returns a
        dictionary with the taxa as keys and sequences as values as well as two
        integers corresponding to the number of taxa and sequence length """

        loci_lengths = [] # Saves the sequence lengths of the
        loci_range = [] # Saves the loci names as keys and their range as values
        main_taxa_order = []

        for infile in input_list:

            # When set to True, this statement produces a progress status on the terminal
            if progress_stat == True:
                print ("\rProcessing file %s out of %s" % (input_list.index(infile)+1,len(input_list)),end="")

            # Parse the current alignment
            current_alignment, taxa_order, current_sequence_len = self.read_alignment(infile,alignment_format)

            # Algorithm that fills absent taxa with missing data
            if loci_lengths == []:
                main_alignment = current_alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
                main_taxa_order = taxa_order
                loci_lengths.append(current_sequence_len)
                loci_range.append((infile.split(".")[0],"1-%s" % (current_sequence_len))) # Saving the range for the first loci

            else:
                for taxa, sequence in current_alignment.items():
                    if taxa in main_alignment:
                        main_alignment[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
                    elif taxa not in main_alignment:
                        main_alignment[taxa] = self.missing*sum(loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment
                        main_taxa_order.append(taxa)

                # Saving the range for the subsequent loci
                loci_range.append((infile.split(".")[0],"%s-%s" % (sum(loci_lengths)+1, sum(loci_lengths)+current_sequence_len)))
                loci_lengths.append(current_sequence_len)

                # Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
                for taxa in main_alignment.keys():
                    if taxa not in current_alignment:
                        main_alignment[taxa] += self.missing*current_sequence_len

        return (main_alignment, main_taxa_order, loci_lengths, loci_range)
