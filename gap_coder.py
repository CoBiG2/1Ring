#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  gap_coder.py
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



def gap_listing (sequence,gap_symbol):
	""" Function that parses a sequence string and returns the position of indel events. The returned list is composed of tuples with the span of each indel """
	gap = "%s+" % (gap_symbol)
	span_regex = ""
	gap_list,seq_start = [],0
	while span_regex != None:
		span_regex = re.search(gap,sequence)
		if span_regex != None and seq_start == 0:
			gap_list.append(span_regex.span())
			sequence = sequence[span_regex.span()[1]+1:]
			seq_start = span_regex.span()[1]+1
		elif span_regex != None and seq_start != 0:
			gap_list.append((span_regex.span()[0]+seq_start,span_regex.span()[1]+seq_start))
			sequence = sequence[span_regex.span()[1]+1:]
			seq_start += span_regex.span()[1]+1
	return gap_list

def gap_binary_generator (sequence,gap_list):
	""" This function contains the algorithm to construct the binary state block for the indel events """
	for cur_gap in gap_list:
		cur_gap_start,cur_gap_end = cur_gap
		if sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start) and sequence[cur_gap_start-1] != "-" and sequence[cur_gap_end] != "-":
			sequence += "1"
		elif sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start):
			if sequence[cur_gap_start-1] == "-" or sequence[cur_gap_end] == "-":
				sequence += "-"
		elif sequence[cur_gap_start:cur_gap_end] != "-"*(cur_gap_end - cur_gap_start):
			sequence += "0"
	return sequence

def multiSeq_gap_listing (storage,gap_symbol):
	complete_gap_list = []
	for taxa in storage:
		temp_list = gap_listing(storage[taxa],gap_symbol)
		unique_gap = [gap for gap in temp_list if gap not in complete_gap_list]
		complete_gap_list += unique_gap
	return complete_gap_list
