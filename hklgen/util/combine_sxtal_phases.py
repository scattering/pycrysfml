#!/usr/bin/python
import sys
class Reflection(self):
    def __init__(self, hkl, intensity, error):
        self.hkl = hkl
        self.intensity = intensity
        self.error = error
    def __eq__(self, other):
        return self.hkl == other.hkl
if len(sys.argv) >= 4:
    if sys.argv[1] == "--help":
        print "Usuage is ./combine_sxtal_phases.py nuclear_phase.int magnetic_phase.int combined_out.int"
    else:
        phase1 = [line.strip().split() for line in open(sys.argv[1])]
        phase2 = [line.strip().split() for line in open(sys.argv[2])]
        outfile = open(sys.argv[3], mode='w')
        Nuclear_reflections = []
        mag_reflections = []
        for i in range(len(phase1)):
            if i > 2:
                reflections.append(Reflection([int(phase1[i][0]), int(phase1[i][1]), int(phase1[i][2])], float(phase1[i][3]), float(phase1[i][4])))
                
else:
    print "Usuage is ./combine_sxtal_phases.py nuclear_phase.int magnetic_phase.int combined_out.int"
