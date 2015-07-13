import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from fswig_hklgen import *
DATAPATH = os.path.dirname(os.path.abspath(__file__))
infoFile = os.path.join(DATAPATH,r"dy.cfl")
SpG, cell, A = readInfo(infoFile)
funcs.write_crystal_cell(cell)
funcs.write_spacegroup(SpG,int_to_p(0),int_to_p(1))
funcs.write_atom_list(A)
msg, mcell, Am, MGp = readMagInfo(infoFile)
Mh = MagReflection()
ih, ik, il, m = input(" => Enter a magnetic reflection as 4 integers -> (h,k,l,m)=H+sign(m)*k(abs(m)): ")
if m == 0 or abs(ih)+abs(ik)+abs(il) == 0: quit()
j = m
sig = "+"
if j == -1: sig="-"
Mh.set_magh_signp(float(-j))
iv = abs(m)
Mh.set_magh_num_k(iv)
Mh.set_magh_h(FloatVector((np.array([ih,ik,il])-Mh.get_magh_signp()*MGp.kvec(iv-1))))
Mh.set_magh_s(calcS(cell, Mh.hkl))
Mh.set_magh_keqv_minus(funcs.k_equiv_minus_k(FloatVector(MGp.kvec(iv-1)), MGp.lattice))
funcs.calc_magnetic_strf_miv(cell, MGp, Am, Mh)
print "  Reflection: (",ih,ik,il,") ", sig," (", MGp.kvec(iv-1), ")"
print "              (",Mh.hkl,")"
print "  Square of Mag. int.  Vector : ", Mh.sqMiV