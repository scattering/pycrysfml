#!/usr/bin/python
# gen_cpp_methods.py
# Script to facilite the manual generation of string get methods 
# for SWIG
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
lines = [line.strip() for line in open("./ff")]
defs = {}
lname = ""
outfile = open("./cpp_modules/string_wrap.cpp", "w")
hfile = open("./cpp_modules/string_wrap.h", "w")
hlines = []
includes = ["#include <iostream>\n", "#include <string>\n"]
outfile.write("#include \"string_wrap.h\"\n")
hfile.write("#ifndef STRING_WRAP\n#define STRING_WRAP\n")
for line in lines:	
	if '#' in line:
		lname = line.split()[1]
		defs[lname] = []
	else:
		defs[lname].append(line)
for k in defs:
	type_name = k.split('!')[0]
	root_name = type_name.lower().split('_type')[0]
	ff = defs[k]
	#ff = [line.strip() for line in open("./ff")]#+raw_input("Enter name of field file: "))]
	fields = {}
	for line in ff:
	    s = line.split("::")
	    ind = ((s[1].strip()).split('!')[0])
	    names = ind.split(',')
	    if len(names) > 1:
		for name in names:
		    fields[name.split('=')[0].strip()] = s[0].strip()
	    else:
		fields[ind.split('=')[0].strip()] = s[0].strip()
	for key in fields:
		if 'character' in fields[key].lower() and 'dimension' not in fields[key].lower():
			hlines.append("std::string get"+root_name.replace("_", " ").title().replace(" ", "")+"_"+key.lower()+"("+type_name.lower()+"* obj);\n")
			include = "#include \""+type_name.lower()+".h\"\n"
			if include not in includes:
				includes.append(include)
			outfile.write("std::string get"+root_name.replace("_", " ").title().replace(" ", "")+"_"+key.lower()+"("+type_name.lower()+"* obj){\n")
			outfile.write("\tstd::string result = \"\";\n")
			outfile.write("\tobj->get_"+root_name+"_"+key.lower()+"(&result);\n")
			outfile.write("\treturn result;\n")
			outfile.write("}\n")
for line in includes:
	hfile.write(line)
for line in hlines:
	hfile.write(line)
outfile.close()
hfile.write("#endif\n")
hfile.close()
