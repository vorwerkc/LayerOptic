#! /usr/bin/python
#generates example input files to compare to isotropic algorithm
import sys
import math
import glob
import numpy as np

hbar=6.58211928*10.**(-16)

w=[370./hbar] #frequency list in Hz
#for j in xrange(0,100):
	#inter=(0.5*j+400.)*10.**(12)
	#w.append(inter)
	

epsilon2=np.array([[1.003,0.,0.],[0.,1.003,0.],[0.,0.,1.003]]) #first layer
epsilon3=np.array([[6.9,0.,0.],[0.,6.9,0.],[0.,0.,6.9]]) #substrat

		
with open('1_11.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon2[0,0], 0.))
		
with open('1_22.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon2[1,1], 0.))

with open('1_33.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon2[2,2], 0.))

with open('2_11.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon3[0,0], 0.))
		
with open('2_22.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon3[1,1], 0.))

with open('2_33.OUT','w') as f:
	for i in xrange(0,len(w)):
		f.write('%g  %1.9e %1.9e\n' % (w[i]*hbar, epsilon3[2,2], 0.))
