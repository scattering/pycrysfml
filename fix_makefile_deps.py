#!/usr/bin/python
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
# hack out EoS and SAn modules from Makefile
lines = [line.strip() for line in open("./Makefile.deps")]
out = open("./Makefile.deps", "w")
for line in lines:
	if 'SOURCES' in line:
		s = line.split()
		for item in s:
			if 'CFML_Optimization_SAn' in item:
				pass
			elif 'CFML_EoS_Mod' in item:
				pass
			else:
				out.write(item + " ")
		out.write("\n")
	else:
		if 'CFML_Optimization_SAn' in line:
			pass
		elif 'CFML_EoS_Mod' in line:
			pass
		else:
			out.write(line+"\n")
out.close()
