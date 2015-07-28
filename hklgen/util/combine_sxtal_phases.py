#!/usr/bin/python
import sys
import numpy as np
class Reflection():
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
        reflections = []
        for i in range(len(phase1)):
            if i > 2:
                reflections.append(Reflection([int(phase1[i][0]), int(phase1[i][1]), int(phase1[i][2])], float(phase1[i][3]), float(phase1[i][4])))
        wavelength = None
        kvec = []
        for j in range(len(phase2)):
            if j == 2:
                wavelength = float(phase2[j][0])
            elif j == 4:
                kvec = np.array([float(phase2[j][1]), float(phase2[j][2]), float(phase2[j][3])])
            elif j > 4:
                reflection = Reflection([int(phase2[j][0]), int(phase2[j][1]), int(phase2[j][2])], float(phase2[j][4]), float(phase2[j][5]))
                if reflection not in reflections:
                    reflections.append(reflection)
                else:
                    reflections[reflections.index(reflection)].intensity += reflection.intensity
                    reflections[reflections.index(reflection)].error = np.sqrt(reflection.error**2+reflections[reflections.index(reflection)].error**2)
        outfile.write(" ".join(phase1[0]) + "\n")
        outfile.write(" ".join([str(n) for n in kvec]) + "\n")
        outfile.write(str(wavelength) + "\n")
        for ref in reflections:
            outfile.write("   "+"   ".join([str(n) for n in ref.hkl])+"   " + str(ref.intensity) + "   " + str(ref.error) + "  1\n")
        outfile.close()
else:
    print "Usuage is ./combine_sxtal_phases.py nuclear_phase.int magnetic_phase.int combined_out.int"
