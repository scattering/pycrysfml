#!/usr/bin/python
lines = [line.strip() for line in open("./wrap/FortFuncs.h")]
out = open("./wrap/FortFuncs.h", "w")
for line in lines:
    if 'extern "C" {' in line:
        out.write("#define None 0\nstruct complex{\n\tdouble real;\n\tdouble imaginary;\n};\n")
    out.write(line+"\n")
out.close()