#!/usr/bin/python
# fix_line_width.py
# corrects Fortran line length issues
from sys import argv
lines = [line.strip() for line in open (argv[1])]
out = open(argv[1], "w")
for line in lines:
	if len(line) > 72:
		s = line.split()
		nl = ""
		for i in range(len(s)):
			if (len(nl)+len(s[i])) < 70:
				nl = nl+" "+s[i]
			else:
				out.write(nl+" &\n")
				nl = s[i]
		out.write(nl+"\n")
	else:
		out.write(line+"\n")
out.close()
