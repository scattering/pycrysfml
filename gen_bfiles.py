#!/usr/bin/python
lines = [line.strip() for line in open("./Makefile.deps")]
line = lines[0].replace("SOURCES= ", "").replace("Src/", "").split()
outfile = open("./bfiles.txt", "w")
for l in line:
    if l == "CFML_Sym_Table.f90" or l == "CFML_Chem_Scatt.f90" or l == "CFML_Bonds_Table":
        outfile.write("$COMP $OPT2 "+l+"\n")
    else:
        outfile.write("$COMP $OPT1 "+l+"\n")
outfile.close()
