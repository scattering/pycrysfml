#!/usr/bin/python3
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
				print(item, end=" ", file=out)
		print(end="\n", file=out)
	else:
		if 'CFML_Optimization_SAn' in line:
			pass
		elif 'CFML_EoS_Mod' in line:
			pass
		else:
			print(line, end="\n", file=out)
out.close()
