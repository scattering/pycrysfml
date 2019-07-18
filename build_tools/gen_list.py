#!/usr/bin/python
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
# Generate file list for use with modified FortWrap based on bfiles.txt
from os import getcwd
lines = [line.strip() for line in open("./bfiles.txt")]
for line in lines:
	if '#' not in line:
		s = line.split(" ")
		s1 = s[2].split(".")
		#compilation list
		#print("$COMP", s[2], "$OPT1", s1[0]+".o", sep=" ")
		#fortwrap list
		#print("/home/jel/TestEnv/crysfml/Src/"+s[2])
		print getcwd()+"/Src/"+s[2]
		#print("../"+s1[0]+".o", end=" ")
		#print(s[2])
