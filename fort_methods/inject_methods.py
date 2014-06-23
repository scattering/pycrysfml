#!/usr/bin/python
# inject_methods.py
# Syntax: inject_methods.py module.f90 additional_methods.f90
from sys import argv
lines = [line.strip() for line in open(argv[1])]
nlns = [line.strip() for line in open(argv[2])]
out = open(argv[1], "w")
for line in lines:
	if 'end module' in line.lower():
		for c in nlns:
			out.write("\t"+c+"\n")
	out.write(line+"\n")
