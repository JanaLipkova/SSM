#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 19:17:07 2017

@author: G. Arampatzis (gtarabat@gmail.com)
"""

import sys


if len(sys.argv) < 4:
    print("\n Usage: ./write_xml  sbml.xml  method  epsilon\n")
    exit(1)

old_xml = sys.argv[1]
method  = sys.argv[2]
epsilon = sys.argv[3]


if old_xml.endswith('.xml'):
    old_xml = old_xml[:-4]


new_xml = old_xml + '_' + method + '_' + epsilon


with open(old_xml+'.xml', 'r') as myfile:
    data=myfile.read()
myfile.close()


k=data.find('model');

cnt = k;
for c in data[k+1:]:
    cnt = cnt + 1
    if(c=='\n'):
        break



with open('annotation', 'r') as myfile:
    ann=myfile.read()
myfile.close()

ann = ann.replace( 'Method=XXX',  'Method="'+method+'"' )
ann = ann.replace( 'Epsilon=XXX', 'Epsilon="'+epsilon+'"' )

print(ann)

    
with open( new_xml + '.xml' , 'w+') as f:
    f.write(data[:cnt])
    f.write('\n\n\n')
    f.write(ann)
    f.write('\n\n\n')
    f.write(data[cnt+1:])

    f.seek(0)
    data = f.read()    
