#!/usr/bin/python
# gen_fort_methods.py
# Script to facilite the manual generation of get and set methods
# as well as constructors in Fortran
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
lines = [line.strip() for line in open("./ff")]
defs = {}
lname = ""
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
	outfile = open("./fort_methods/"+k.split('!')[1]+"/"+type_name+".f90", "w")
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
		if not any(word in fields[key].lower() for word in ['character', 'type', 'dimension']):
		    # write getter
		    outfile.write("function get_"+root_name+"_"+key+"(obj_var)\n")
		    outfile.write("\ttype ("+type_name+") :: obj_var\n")
		    outfile.write("\t"+fields[key]+" :: get_"+root_name+"_"+key+"\n")
		    outfile.write("\tget_"+root_name+"_"+key+" = obj_var%"+key+"\n")
		    outfile.write("end function get_"+root_name+"_"+key+"\n\n")
		else:
		    # write getter
		    # use subroutine and modify arguments for strings, derived types, and arrays
		    outfile.write("subroutine get_"+root_name+"_"+key+"(obj_var, output_value)\n")
		    outfile.write("\ttype ("+type_name+") :: obj_var\n")
		    outfile.write("\t"+fields[key].replace("allocatable,", "")+", intent(out) :: output_value\n")
		    outfile.write("\toutput_value = obj_var%"+key+"\n")
		    outfile.write("end subroutine get_"+root_name+"_"+key+"\n\n")
		# write setter
		outfile.write("subroutine set_"+root_name+"_"+key+"(obj_var, new_value)\n")
		outfile.write("\ttype ("+type_name+") :: obj_var\n")
		outfile.write("\t"+fields[key].replace("allocatable,", "")+", intent(in) :: new_value\n")
		outfile.write("\tobj_var%"+key+" = new_value\n")
		outfile.write("end subroutine set_"+root_name+"_"+key+"\n\n")
	# write constructor
	outfile.write("subroutine "+type_name+"_ctor("+type_name+"_param, ")
	count = 0
	for key in fields:
	    count += 1
	    if count == len(fields):
		outfile.write(key+"_param)\n")
	    else:
		outfile.write(key+"_param, ")
	outfile.write("\ttype ("+type_name+") :: "+type_name+"_param\n")
	for key in fields:
	    outfile.write("\t"+fields[key].replace("allocatable,", "")+", intent(in) :: "+key+"_param\n")
	for key in fields:
	    outfile.write("\t"+type_name+"_param%"+key+" = "+key+"_param\n")
	outfile.write("end subroutine "+type_name+"_ctor\n")
	outfile.close()
