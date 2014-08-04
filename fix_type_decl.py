#!/usr/bin/python
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
# Fix derived type declarations for use with modified FortWrap
import os
import glob
path = '.'
for filename in glob.glob(os.path.join(path, '*.f90')):
	pub_types = []
	priv_types = []
	lines = [line.strip() for line in open(filename)]
	out = open(filename, 'w')
	for i in range(len(lines)):
		lines[i] = lines[i].replace(",", ", ")
	for line in lines:
		if 'type, public ::' in line:
			pub_types.append(line.split('type, public :: ')[1])
			l1 = line.replace('type, public ::', 'type')
			out.write(l1)
		elif 'type, private ::' in line:
			priv_types.append(line.split('type, private :: ')[1])
			l1 = line.replace('type, private ::', 'type')
			out.write(l1)
		elif 'type,public ::' in line:
			pub_types.append(line.split('type,public :: ')[1])
			l1 = line.replace('type,public ::', 'type')
			out.write(l1)
		elif 'type,private ::' in line:
			priv_types.append(line.split('type,private :: ')[1])
			l1 = line.replace('type,private ::', 'type')
			out.write(l1)
		elif 'type' in line and ('public' in line or 'private' in line) and '::' in line and '(' not in line and '&' not in line:
			name = line.split('::')[1].strip()
			l1 = line.replace(line[line.find('type'):line.find('::')+2], 'type')
			if 'public' in line:
				pub_types.append(name)
			elif 'private' in line:
				priv_types.append(name)	
			out.write(l1)
		else:
			out.write(line)
		out.write('\n')
	out.close()
	lines = [line.strip() for line in open(filename)]
	out = open(filename, 'w')
	pub_types = list(set(pub_types))
	priv_types = list(set(priv_types))
	added = False
	#for line in lines:
		#if 'type' in line and '!' not in line and not added:
			#if len(pub_types) > 0:
				#out.write('public :: ')
				#for i in range(len(pub_types)):
					#if i == len(pub_types)-1:
						#out.write(pub_types[i])
					#else:
						#out.write(pub_types[i]+", ")
				#out.write('\n')
				#added = True
			#if len(priv_types) > 0:
				#out.write('private :: ')
				#for i in range(len(priv_types)):
					#if i == len(priv_types)-1:
						#out.write(priv_types[i])
					#else:
						#out.write(priv_types[i]+", ")
				#out.write('\n')
				#added = True
		#out.write(line)
		#out.write('\n')
	#out.close()
	for line in lines:
		out.write(line+"\n")
		if 'module' in line and 'end module' not in line and '!' not in line:
			if len(pub_types) > 0:
				out.write('public :: ')
				for i in range(len(pub_types)):
					if i == len(pub_types)-1:
						out.write(pub_types[i])
					else:
						out.write(pub_types[i]+", ")
				out.write('\n')
			if len(priv_types) > 0:
				out.write('private :: ')
				for i in range(len(priv_types)):
					if i == len(priv_types)-1:
						out.write(priv_types[i])
					else:
						out.write(priv_types[i]+", ")
				out.write('\n')
	out.close()
