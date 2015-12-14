#! /usr/bin/env python
from __future__ import division
import sys

class SeqObj():
	def __init__(self, name, norm_length, read_cov):
	    self.name = name
	    self.norm_length = norm_length
	    self.read_cov = read_cov
	    self.norm_cov = read_cov/norm_length if (read_cov/norm_length) > 0.0001 else 0.0001

def parse_norm_length(norm_length_f):
	norm_length_d = {}
	with open(norm_length_f) as fh:
		for l in fh:
			name, length = l.rstrip("\n").split("\t") 
			norm_length_d[name] = int(length)
	return norm_length_d

def parse_depth(rnaseq_read_cov_f):
	seq_obj_d = {}
	read_cov_d = {}
	with open(rnaseq_read_cov_f) as fh:
		for l in fh:
			name, pos, depth = l.rstrip("\n").lstrip(" ").split()
			read_cov_d[name] = read_cov_d.get(name, 0) + int(depth)
	for name in norm_length_d:
		if name in read_cov_d:
			seq_obj_d[name] = SeqObj(name, norm_length_d[name], read_cov_d[name])
		else:
			seq_obj_d[name] = SeqObj(name, norm_length_d[name], 0)
	return seq_obj_d

if __name__ == "__main__":
	norm_length_f = sys.argv[1]
	rnaseq_read_cov_f = sys.argv[2]
	norm_length_d = parse_norm_length(norm_length_f)
	seq_obj_d = parse_depth(rnaseq_read_cov_f)
	for name, seqObj in seq_obj_d.items():
		print "%s\t%s\t%s\t%s" % (seqObj.name, seqObj.norm_length, seqObj.read_cov, seqObj.norm_cov)