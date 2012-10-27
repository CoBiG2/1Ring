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


## TODO:
## ARLEQUIN
## IMA2/IMA

class writer ():
		
	def __init__ (self, output_file, taxa_order, coding, loci_lengths, loci_range=None, gap = "-", missing = "n", conversion = None):
		self.output_file = output_file
		self.taxa_order = taxa_order
		self.coding = coding
		self.loci_lengths = loci_lengths
		self.loci_range = loci_range
		#print (self.loci_range)
		self.gap = gap
		self.missing = missing
		# The space (in characters) available for the taxon name before the sequence begins
		self.seq_space_nex = 40
		self.seq_space_phy = 30
		self.seq_space_ima2 = 10
		# Cut the taxa names by the following character:
		self.cut_space_nex = 50
		self.cut_space_phy = 50
		self.cut_space_ima2 = 8
			
	def phylip (self, alignment_dic, conversion=None):
		""" Writes a pre-parsed alignment dictionary into a new phylip file """
		out_file = open(self.output_file+".phy","w")
		out_file.write("%s %s\n" % (len(alignment_dic), self.loci_lengths))
		for key in self.taxa_order:
				out_file.write("%s %s\n" % (key[:self.cut_space_phy].ljust(self.seq_space_phy),alignment_dic[key]))
		if conversion == None:
			partition_file = open(self.output_file+"_part.File","a")
			for partition,lrange in self.loci_range:
				partition_file.write("%s, %s = %s\n" % (self.coding,partition,lrange))
		out_file.close()
				
	def fasta (self, alignment_dic, conversion=None):
		""" Writes a pre-parsed alignment dictionary into a new fasta file """
		out_file = open(self.output_file+".fas","w")
		for key in self.taxa_order:
			out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))
		out_file.close()
			
	def nexus (self, alignment_dic, conversion=None):
		""" Writes a pre-parsed alignment dictionary into a new nexus file """
		out_file = open(self.output_file+".nex","w")
		out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment_dic), self.loci_lengths, self.coding, self.gap, self.missing))
		for key in self.taxa_order:
			out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),alignment_dic[key]))
		out_file.write(";\n\tend;")
		if conversion == None:
			out_file.write("\nbegin mrbayes;\n")
			for partition,lrange in self.loci_range:
				out_file.write("\tcharset %s = %s;\n" % (partition,lrange))
			out_file.write("\tpartition part = %s: %s;\n\tset partition=part;\nend;" % (len(self.loci_range),", ".join([part[0] for part in self.loci_range])))
		out_file.close()
			
	def zorro (self, zorro_weigths):
		""" Creates a concatenated file with the zorro weigths for the corresponding alignment files """
		outfile = self.output_file+"_zorro.out"
		outfile_handle = open(outfile,"w")
		for weigth in zorro_weigths:
			outfile_handle.write("%s\n" % weigth)
		outfile_handle.close()
