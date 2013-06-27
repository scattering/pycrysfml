#!/usr/bin/python

import re, sys, os

use_line_re = re.compile("^\s*use\s+([a-z0-9_]+)\s*(,|!|$)", re.IGNORECASE)
mod_line_re = re.compile("^\s*module\s+([a-z0-9_]+)\s*(!|$)", re.IGNORECASE)

def get_deps(filename):
    mod = ""
    deps = []
    with open(filename) as fid:
        for line in fid:
            match = use_line_re.search(line)
            if match:
                 deps.append(match.group(1).lower())
                 continue
            match = mod_line_re.search(line)
            if match:
                  mod = match.group(1).lower()
                  continue    
    #print filename, mod, deps
    return mod, deps

def new_ext(filename, ext): return os.path.splitext(filename)[0]+ext

def process_files(files, target_ext=".o"):
    mod_file = {}
    mod_deps = {}
    for f in files:
        mod, deps = get_deps(f)
        if mod in mod_file: 
            print >>sys.stderr, "replacing",mod_file[mod],"with",f
        if mod: mod_file[mod] = f
        if mod and deps: mod_deps[mod] = deps

    print "SOURCES="," ".join(sorted(mod_file.values()))
    print
    #print "\n".join(k+": "+v for k,v in sorted(mod_file.items()))
    for mod, deps in sorted(mod_deps.items()):
        target = new_ext(mod_file[mod], target_ext)
        print target, ":", " ".join(new_ext(mod_file[m],target_ext) 
                                    for m in deps if m in mod_file)

if __name__ == "__main__":
    process_files(sys.argv[1:])
