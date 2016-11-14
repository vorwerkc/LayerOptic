#! /usr/bin/env python
#
#
# AUTHOR: Christian Vorwerk 2015
# 		  vorwerk(at)physik.hu-berlin.de
# 
# 
# EXPLANATION:
#  applies Euler rotation to the dielectric tensor of a given layer
# 	
#__________________________________________________________________________________________________
import os
import sys
import string
import math
import glob
import copy
import os.path
import numpy as np
from math import sin, cos, pi


if len(sys.argv)>1 and  len(sys.argv) !=5:
	print '\n Please provide FOUR input parameters:'\
		  '\n 3 Euler Angles in Degree'\
		  '\n Layer Number'
	sys.exit()
		  
# Euler Angles
if len(sys.argv)>1:
	alpha=float(sys.argv[1])*math.pi/180.
	beta=float(sys.argv[2])*math.pi/180.
	gamma=float(sys.argv[3])*math.pi/180.
	L=sys.argv[4] # Layer number
else:
	print '\n       _________________________________________________________________ '\
	      '\n      |                                                                 |'\
	      '\n      |                 WELCOME TO LAYEROPTICS (SETUP)                  |'\
	      '\n      |                                                                 |'\
	      '\n      |                 written by Christian Vorwerk                    |'\
	      '\n      |_________________________________________________________________|'\
	
	print '\n                        Please provide the layer index:                  '
	L=raw_input('>>>>>>>>     ')
	print '\n             Please provide the Euler angle ALPHA (deg):             '
	alpha=raw_input('>>>>>>>>     ')
	alpha=float(alpha)*math.pi/180.
	print '\n             Please provide the Euler angle BETA (deg):             '
	beta=raw_input('>>>>>>>>     ')
	beta=float(beta)*math.pi/180.
	print '\n             Please provide the Euler angle GAMMA (deg):             '
	gamma=raw_input('>>>>>>>>     ')
	gamma=float(alpha)*math.pi/180.

	if alpha=='RANGE':
		print '\n          Please type the two values between which the angle '\
			  '\n             should be varied and the number of steps.        '                          
		alpha=raw_input('>>>>>>>>     ')
		if len(alpha.split(" ")) > 3 or len(alpha.split(" ")) < 3 :
			sys.exit('     Please provide two real  numbers and one integer!')
		alpha0=float(alpha.split(" ")[0])*math.pi/180.
		alpha1=float(alpha.split(" ")[1])*math.pi/180.
		N=int(alpha.split(" ")[2])
		alpha=[]
		for i in xrange(0,N):
			alpha.append(alpha0+i*(alpha1-alpha0)/(N-1))
		print '\n           Choose the polarization angle.'\
			  '\n      For an array of beam angles the polarization'\
			  '\n                    has to be fixed'
		beta=raw_input('>>>>>>>>     ')
		beta=float(beta)*math.pi/180.

#Sine and Cosine of Angles
c1=math.cos(alpha)
c2=math.cos(beta)
c3=math.cos(gamma)

s1=math.sin(alpha)
s2=math.sin(beta)
s3=math.sin(gamma)

#define Rotation Matrix A
A=np.zeros((3,3))
A[0,0]=c1*c3-c2*s2*s3
A[0,1]=-c1*s3-c2*c3*s1
A[0,2]=s1*s2

A[1,0]=c3*s1+c1*c2*s3
A[1,1]=c1*c2*c3-s1*s3
A[1,2]=-c1*s2

A[2,0]=s2*s3
A[2,1]=c3*s2
A[2,2]=c2

#calculate A^[-1]
A_inv=np.linalg.inv(A)

# read the full dielectric tensor
inlist=[]
WORKDIR=os.getcwd()
files=os.listdir(WORKDIR)
for i in xrange(0,len(files)):
	if files[i].split('_')[0]=='EPSILON':
		inlist.append(files[i])

# read all frequency lists and elements of the tensors into dictionary
scheme=''
for j in xrange(0,len(inlist[0].split('_'))-1):
	
	scheme=scheme+inlist[0].split('_')[j]+'_'
scheme=scheme+'OC'
#scheme=inlist[0].split('_')[0]+'_'+inlist[0].split('_')[1]+'_'+inlist[0].split('_')[2]+'_OC'
data={}
for i in xrange(0,3):
	for j in xrange(0,3):
	
		name=scheme+str(i+1)+str(j+1)+'.OUT'
		if os.path.isfile(name)== False:
			data['w'+str(i)+str(j)]=[]
			data[str(i)+str(j)]=[]
		else:
			EPS=open(name,'r')
			data['w'+str(i)+str(j)]=[]
			data[str(i)+str(j)]=[]
			while True:
				line=EPS.readline()
				line=line.strip()
				if len(line)==0: break
				data['w'+str(i)+str(j)].append(float(line.split()[0]))
				data[str(i)+str(j)].append(float(line.split()[1])+float(line.split()[2])*1j)
		
# form list of dielectric tensors for each frequency from directory
epsilon=[]
w=[]
for i in xrange(0,len(data['w00'])):
	inter=np.zeros((3,3),'complex')
	for j in xrange(0,3):
		inter[j,j]=1.
	w.append(data['w00'][i])		
	for k in xrange(0,3):
		for l in xrange(0,3):
			if len(data[str(k)+str(l)]) != 0:
				inter[k,l]=data[str(k)+str(l)][i]
	epsilon.append(inter)

#perform rotation
epsilon2=[]
for i in xrange(0,len(epsilon)):
	epsilon[i]=np.dot(A,np.dot(epsilon[i],A_inv))

#write epsilon back to files
for i in xrange(0,3):
	for j in xrange(0,3):
		with open(str(L)+'_'+str(i+1)+str(j+1)+'.OUT','w') as f:
			for k in xrange(0,len(epsilon)):
				if w[k] != 0.0:
					f.write('%g  %1.9e %1.9e\n' % (w[k], epsilon[k][i,j].real, epsilon[k][i,j].imag))

			
	


