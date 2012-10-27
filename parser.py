#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ElParsito.py v3.0.0-5
#
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

## TODO
## PARSE CONCATENATED FILES (NEXUS, PHYLIP)
## AUTORECOGNITION OF SEQUENCE CODE (DNA/PROTEIN/31|\|4RY)

from pympler.asizeof import asizeof

class SeqUtils ():
    def __init__ (self, missing="X"):
        self.missing = missing

    def rm_illegal (self,name):
    """ Removes any 'illegal' charaters from the taxas' names. """
        warning = ""
        chars = set(" ",":",",",")","(",";","[","]","'")
        newname = []
        for i in name:
            if i in chars:
                newname.append("_")
            else:
                newname.append(i)
        newname = "".join(newname)
        if newname != name:
        #Suggestion - the module should not print any messeges. These should be
        #returned to the main program and let it handle them.
            warning = "WARNING: Replaced illegal characters from the taxa %s" % name
        return newname, warning

    def duplicate_taxa (self, taxa_list):
        """ Function that identifies repeats in taxa names """
        from collections import Counter
        duplicated_taxa = [x for x, y in Counter(taxa_list).items() if y > 1]
        return duplicated_taxa #Returns a list with the names of repeated taxa

    def check_format (self,input_alignment,alignment_format):
        """ This function performs some very basic checks to see if the format
        of the input file is in accordance to the input file format specified
        when the script is executed """
        input_handle = open(input_alignment)
        line = input_handle.readline()
        while line.strip() == "":
            line = next(input_handle)

        if alignment_format == "fasta":
            if line.strip()[0] != ">":
                print ("File not in Fasta format. First non-empty line of the input file %s does not start with '>'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
                raise SystemExit
        elif alignment_format == "nexus":
            if line.strip().lower() != "#nexus":
                print ("File not in Nexus format. First non-empty line of the input file %s does not start with '#NEXUS'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
                raise SystemExit
        elif alignment_format == "phylip":
            try:
                header = line.strip().split()
                int(header[0])
                int(header[1])
            except:
                print ("File not in correct Phylip format. First non-empty line of the input file %s does not start with two intergers separated by whitespace. Please verify the file, or the input format settings\nExiting..." % input_alignment)
                raise SystemExit

    def autofinder (self, infile_name):
        #Autodetects the type of file to be parsed. Based on headers.
        autofind = "unknown"
        infile = open(infile_name,'r')
        header = infile.readline()
        while header.startswith("\n"):
            header = next(infile)
        infile.close()

        if header.upper().startswith("#NEXUS"):
            autofind = "nexus"
            break
        elif header.startswith(">"):
            autofind = "fasta"
            break
        elif len(header.strip().split()) == 2 and phy_header[0].isdigit():
            autofind = "phylip"
            break

        return autofind

    def rm_taxa (self, alignment_dic, taxa_list):
        """ Function that removes specified taxa from the alignment """
        alignment_mod = {}
        taxa_order = []

        for taxa, sequence in alignment_dic.items():
            if taxa not in taxa_list:
                alignment_mod[taxa] = sequence
                taxa_order.append(taxa)

        return alignment_mod, taxa_order

    def pickle_taxa (self, alignment_dic, mode):
        """ Function that exports the list of taxa from an alignment """
        import pickle
        self.taxa_list = []
        if mode == "dump":
            self.taxa_list = [taxa for taxa in alignment_dic.keys()]
            pickle.dump(self.taxa_list, open("taxa_list","wb"))
            print ("Taxa names have been saved in the pickle file 'taxa_list'\nExiting...")
            raise SystemExit
        elif mode == "load":
            self.taxa_list = pickle.load(open("taxa_list","rb"))
        return self.taxa_list

    def import_taxa (self, alignment_dic):
        """ Function that imports new taxa. It mainly exists to complete single
        locus aligments with taxa that are not present in the current alignment
        but occur in other alignments """
        alignment_len = self.loci_lengths[0]
        for taxa in self.taxa_list:
            if taxa not in alignment_dic:
                alignment_dic[taxa] = self.missing*alignment_len
        return alignment_dic, self.taxa_list

    def check_sizes (self, Dict, current_file):
        warning = ""
        length = 0
        for i in Dict.values():
            if length != 0 and len(i) != length:
                print(length)
                print(len(i))
                warning = "Not all of your sequences have the same length.\nYou\
really should look into this as it is a VERY BAD sign that something is wrong \
if you are using these sequences for further analyses."
            length = len(i)
        return warning

    def zorro2rax (self, alignment_file_list, zorro_sufix="_zorro.out"):
        """ Function that converts the floating point numbers contained in the
        original zorro output files into intergers that can be interpreted by
        RAxML. If multiple alignment files are provided, it also concatenates
        them in the same order """
        weigths_storage = []
        for alignment_file in alignment_file_list:
            zorro_file = alignment_file.split(".")[0]+zorro_sufix # This assumes that the prefix of the alignment file is shared with the corresponding zorro file
            zorro_handle = open(zorro_file)
            weigths_storage += [round(float(weigth.strip())) for weigth in zorro_handle]
        return weigths_storage

    def read_alignment (self, input_alignment, alignment_format, size_check=True):
        """ ONLY FOR SINGLE FILE/LOCI INPUT: Function that parses an input file
        alignment and returns a dictionary with the taxa as keys and sequences
        as values """

        self.check_format (input_alignment, alignment_format)

        alignment_storage = {} # Save the taxa and their respective sequences
        taxa_order = [] # Save taxa names to maintain initial order
        file_handle = open(input_alignment)

        # PARSING PHYLIP FORMAT
        if alignment_format == "phylip":
            header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
            self.loci_lengths = int(header[1])
            for line in file_handle:
                if line != "":
                    taxa = line.split()[0].replace(" ","")
                    taxa = self.rm_illegal(taxa)
                    taxa_order.append(taxa)
                    sequence = line.split()[1].strip()
                    alignment_storage[taxa] = sequence

        # PARSING FASTA FORMAT
        elif alignment_format == "fasta":
            for line in file_handle:
                if line.strip().startswith(">"):
                    taxa = line[1:].strip().replace(" ","_")
                    taxa = self.rm_illegal(taxa)
                    taxa_order.append(taxa)
                    alignment_storage[taxa] = ""
                else:
                    alignment_storage[taxa] += line.strip()
            self.loci_lengths = len(list(alignment_storage.values())[0])

        # PARSING NEXUS FORMAT
        elif alignment_format == "nexus":
            counter = 0
            for line in file_handle:
                if line.strip().lower() == "matrix" and counter == 0: # Skips the nexus header
                    counter = 1
                elif line.strip() == ";" and counter == 1: # Stop parser here
                    counter = 0
                elif line.strip() != "" and counter == 1: # Start parsing here
                    taxa = line.strip().split()[0].replace(" ","")
                    taxa = self.rm_illegal(taxa)
                    taxa_order.append(taxa)
                    if taxa in alignment_storage: # This accomodates for the interleave format
                        alignment_storage[taxa] += "".join(line.strip().split()[1:])
                    else:
                        alignment_storage[taxa] = "".join(line.strip().split()[1:])
            self.loci_lengths = len(list(alignment_storage.values())[0])

        # Checks the size consistency of the alignment
        if size_check == True:
            self.check_sizes (alignment_storage, input_alignment)

        # Checks for duplicate taxa
        if len(taxa_order) != len(set(taxa_order)):
            taxa = self.duplicate_taxa(taxa_order)
            print ("WARNING: Duplicated taxa have been found in file %s (%s). Please correct this problem and re-run the program\n" %(input_file,", ".join(duplicated_taxa)))
            raise SystemExit

        return (alignment_storage, taxa_order, self.loci_lengths, None)
