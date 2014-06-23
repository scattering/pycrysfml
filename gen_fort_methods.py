#!/usr/bin/python
# gen_fort_methods.py
# Script to facilite the manual generation of get and set methods
# as well as constructors in Fortran
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
type_name = raw_input("Enter type name: ")
root_name = type_name.lower().split('_type')[0]
ff = [line.strip() for line in open("./ff")]#+raw_input("Enter name of field file: "))]
fields = {}
outfile = open("./fort_methods/"+type_name+".f90", "w")
#field_name = ""
#field_type = ""
#data_types = ['integer', 'logical', 'character', 'real', 'complex']
#def check_type(t):
    #if 'type' in t:
        #return True
    #for c in data_types:
        #if c in t:
            #return True
    #return False
#while field_name != "STOP":
    #field_name = raw_input("Enter field name (STOP when done): ")
    #if field_name == "STOP":
        #break
    #field_type = raw_input("Enter field type: ")
    #if not check_type(field_type.lower()):
        #print "invalid data type: "+field_type
        #continue
    #fields[field_name] = field_type
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
    # write getter
    outfile.write("function get_"+root_name+"_"+key+"(obj_var)\n")
    outfile.write("\ttype ("+type_name+") :: obj_var\n")
    outfile.write("\t"+fields[key]+" :: get"+key+"\n")
    outfile.write("\tget"+key+" = obj_var%"+key+"\n")
    outfile.write("end function get_"+root_name+"_"+key+"\n\n")
    # write setter
    outfile.write("subroutine set_"+root_name+"_"+key+"(obj_var, new_value)\n")
    outfile.write("\ttype ("+type_name+") :: obj_var\n")
    outfile.write("\t"+fields[key]+", intent(in) :: new_value\n")
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
    outfile.write("\t"+fields[key]+", intent(in) :: "+key+"_param\n")
for key in fields:
    outfile.write("\t"+type_name+"_param%"+key+" = "+key+"_param\n")
outfile.write("end subroutine "+type_name+"_ctor\n")
outfile.close()
