#!/usr/bin/python
import pycrysfml
import sys
import pydoc

def output_help_to_file(filepath, request):
    f = file(filepath, 'w')
    sys.stdout = f
    pydoc.help(request)
    f.close()
    sys.stdout = sys.__stdout__
    return
output_help_to_file("../help.txt", pycrysfml)
