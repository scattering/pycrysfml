#!/usr/bin/env python

from ctypes import *
flib = cdll["./fexample.so"]
#clib = ctypes.cdll["./cexample.so"]


class Sym_Oper_Type(Structure):
    _fields_ = [("Rot", c_int*9), ("Tr", c_float*3)]

class Reflection(Structure):
    _fields_ = [("hkl", c_int*3), ("multip", c_int), ("s", c_float)]

class vector(Structure):
    _fields_ = [("x", c_int), ("y", c_int), ("z", c_int)]
    
    def __repr__(self):  return "(%d,%d,%d)"%(self.x,self.y,self.z)

class vectorN(Structure):
    _fields_ = [("x", c_int*5)]
    
    def __repr__(self):  return "(%d)"%(self.x[0])

c_size_t = c_ulong
class DVDim(Structure):
    _fields_ = [("stride_mult", c_size_t), ("lower_bound", c_size_t),
                ("upper_bound", c_size_t)]

class DV(Structure):
    _fields_ = [("base_addr",c_void_p), ("base", c_void_p), ("dtype",c_size_t),
                ("dim", DVDim*7)]
                
    def elemsize(self):
        return self.dv.dtype >> 6
    def rank(self):
        return self.dv.dtype&7
    def __len__(self):
        return self.dim[0].upper_bound - self.dim[0].lower_bound + 1
    def data(self, datatype):
        return (datatype*len(self)).from_address(self.base_addr)

def dv_dtype(size,type,rank): return size*64+type*8+rank
def build_struct_dv(array):
    dv = DV()
    dv.base_addr = ct.addressof(array)
    dv.base = c_void_p()
    dv.dtype = dv_dtype(ct.sizeof(array[0]), 5, 1) # derived type
    dv.dim[0].stride_mult = 1
    dv.dim[0].lower_bound = 1
    dv.dim[0].upper_bound = len(array)
    return dv

def deconstruct_dv(dv, dataType, n):
#    dims = dv.dtype & 7
#    elem_size = dv.dtype >> 6
#    size = dv.dim[0].upper_bound - dv.dim[0].lower_bound + 1
#    array = [c_float() for i in xrange(size)]
#    for i in xrange(size):
#        offset = i * elem_size * dv.dim[0].stride_mult
#        array[i] = dataType.from_address(dv.base_addr + offset)
#    return array
    return (dataType*N).from_address(dv.base_addr)

class Thing(Structure):
    _fields_ = [("size",c_int),("array",DV)]

class RefList(Structure):
    _fields_ = [("numRefs", c_int), ("refs", DV)]

def testDV():
    fn = flib.__struct_MOD_makething
    fn.argtypes = [POINTER(Thing)]
    fn.restype = None
    thing = Thing()
    fn(thing)
    fn(thing)
    print "size attribute = " + str(thing.size)
    print "array = " + str(thing.array)
    dim = thing.array.dim[0]
    print "array size = " + str(dim.upper_bound - dim.lower_bound + 1)
    print "data = " + str(thing.array.data(c_float)[:])
    print
    
    fn = flib.__struct_MOD_makereflist
    fn.argtypes = [POINTER(RefList)]
    fn.restype = None
    refList = RefList()
    fn(refList)
    print "size attribute = " + str(refList.numRefs)
    print "refs = " + str(refList.refs)
    dim = refList.refs.dim[0]
    print "array size = " + str(dim.upper_bound - dim.lower_bound + 1)
    reflections = refList.refs.data(Reflection)
    print "data = " + str(reflections)
    print "hkls: "
    for ref in reflections: print ref.hkl[:], ref.s
    return

def returned_pointer_demo():
    
    make = flib.__struct_MOD_make
    #make = clib.makepcopy
    make.restype = Sym_Oper_Type
    make.argtypes = []
    #made = clib.made
    made = flib.__struct_MOD_made
    made.argtypes = [POINTER(Sym_Oper_Type)]
    made.restype = c_int

    sym = make()
    print "0x%x"%sym.Rot[0]
    print "42 if returning a copy; 43 if returning the original", made(sym)

def passed_structure_demo():
    add = flib.__struct_MOD_add
    #add = clib.add
    add.restype = vector
    add.argtypes = [POINTER(vector), POINTER(vector)]
    #res = vector()
    a,b = vector(1,2,3), vector(4,5,6)
    #add(byref(res), byref(a), byref(b))
    res = add(a, b)
    print a,"+",b,"=",res
    
def length(s,n):
    fn = flib.__struct_MOD_length
    fn.argtypes = [c_char_p, POINTER(c_int), c_int]  
    fn.restype = c_int
    return fn(s, c_int(n), len(s))
    
def length2(s,n=None):
    fn = flib.__struct_MOD_length2
    fn.argtypes = [c_char_p, POINTER(c_int), c_int]  
    #fn.argtypes = [c_char_p, c_int]
    fn.restype = c_int
    return fn(s, c_int(n) if n is not None else None, len(s))

def unit(v):
    fn = flib.__struct_MOD_unitvector
    fn.argtypes = [POINTER(vector)]
    fn.restype = None
    fn(v)
    return v

def structret_demo():
    fn = clib.structret
    fn.restype = POINTER(vectorN);
    fn.argtypes = []
    res = fn()
    print "res",res

def charIndex(c):
    fn = flib.__struct_MOD_charindex
    fn.argtypes = [POINTER(c_char)]
    fn.restype = c_int
    return fn(c)

def arraySum(a):
    fn = flib.__struct_MOD_arraysum
    intArray = c_int*len(a)
    fn.argtypes = [POINTER(intArray), c_int]
    fn.restype = c_int
    return fn(intArray(*a), len(a))

def logicTest(b):
    fn = flib.__struct_MOD_logictest
    fn.argtypes = [POINTER(c_bool)]
    fn.restype = c_bool
    return fn(c_bool(b))

# hklUni: constructs a list of unique reflections in a specified range
#   If code == "r", then d-spacings are used
#   If noOrder is True, then the reflections are not sorted by position
def hklUni(sMin, sMax, code):
    MAX_REFLS = 10
    c_ReflectionArray = Reflection*MAX_REFLS
    crefl = c_ReflectionArray()
    fn = flib.__struct_MOD_hkl_uni_reflect
    fn.argtypes = [POINTER(c_bool),
                   POINTER(c_float), POINTER(c_float), 
                   c_char_p,
                   POINTER(c_int), 
                   POINTER(DV), 
                   POINTER(c_bool),
                   ]
    fn.restype = None
    fn(c_bool(True), c_float(sMin), c_float(sMax),
       code,
       c_int(0),
       build_struct_dv(crefl),
       None,
       )
    print "multip returned",crefl[0].multip

def hklrefs():
    return refs


if __name__ == "__main__":
    testDV()
    #hklUni(0,1,'s')
    #print logicTest(True), logicTest(False)
    #print arraySum(range(5))
    #print charIndex('a')
    #v = vector(0,0,0)
    #print unit(v)
    #print length2("blah")
    #print length2("blah", 3)
    #returned_pointer_demo()
    #passed_structure_demo()
    #structret_demo()