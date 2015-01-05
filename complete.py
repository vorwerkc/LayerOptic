#! /usr/bin/python
# test version 1.3.1
# calculates transmission and reflection coefficients for electric fields
# following Puschning & Ambrosch-Draxl(2006), Yeh (1980)
import sys
import os
import math
import cmath
import glob
import time
import numpy as np
import scipy as sp

# DEFINE INPUT PARAMETERS
# General Parameters
c=299792458.0
mu=math.pi*10.**(-7)
hbar=6.58211928*10.**(-16)

#-------------------FUNCTIONS FOR MATRIX--------------
#-----------------------ALGORITHM---------------------
def quartic_solve(alpha, beta, epsilon, omega):
# Solves determinant to obtain kappa
	
	w=omega/c
	A=w**2*epsilon[0,0]-beta**2
	B=w**2*epsilon[0,1]+alpha*beta
	C=w**2*epsilon[0,2]
	D=w**2*epsilon[1,0]+alpha*beta
	E=w**2*epsilon[1,1]-alpha**2
	F=w**2*epsilon[1,2]
	G=w**2*epsilon[2,0]
	H=w**2*epsilon[2,1]
	I=w**2*epsilon[2,2]-alpha**2-beta**2
	

	coeff=[]
	coeff.append(I+alpha**2+beta**2)
	coeff.append(alpha*G+alpha*C+beta*H+beta*F)
	coeff.append(-A*I-E*I+alpha*beta*B+alpha*beta*D+C*G-alpha**2*E-beta**2*A+F*H)
	coeff.append(alpha*B*F+beta*B*G+alpha*D*H+beta*D*C-alpha*E*G-alpha*C*E-beta*A*F-beta*A*H)
	coeff.append(A*E*I+B*F*G+D*H*C-C*E*G-B*D*I-A*F*H)
	
	gamma1=np.roots(coeff)
	gamma=np.zeros((4),'complex')
	for i in xrange(0,4):
		if gamma1[i]< 0.0 and gamma[0]==0.+0.j:
			gamma[0]=gamma1[i]
		elif gamma1[i]< 0.0 and gamma[0]!=0.+0.j:
			gamma[2]=gamma1[i]
		elif gamma1[i]> 0.0 and gamma[1]==0.+0.j:
			gamma[1]=gamma1[i]
		elif gamma1[i]> 0.0 and gamma[1]!=0.+0.j:
			gamma[3]=gamma1[i]

	return gamma
#-----------------------------------------------------------------
def comp_eig(alpha, beta, gamma, epsilon, omega):
# forms the Maxwell matrix to obtain eigenvalues and -vectors
	w=omega/c
	A=np.zeros((3,3), 'complex')
	A[0,0]=w**2*epsilon[0,0]-beta**2-gamma**2
	A[0,1]=w**2*epsilon[0,1]+alpha*beta
	A[0,2]=w**2*epsilon[0,2]+alpha*gamma
	A[1,0]=w**2*epsilon[1,0]+alpha*beta
	A[1,1]=w**2*epsilon[1,1]-alpha**2-gamma**2
	A[1,2]=w**2*epsilon[1,2]+beta*gamma
	A[2,0]=w**2*epsilon[2,0]+alpha*gamma
	A[2,1]=w**2*epsilon[2,1]+beta*gamma
	A[2,2]=w**2*epsilon[2,2]-alpha**2-beta**2
	eigval, eigvector= np.linalg.eig(A)
	 
	return eigval, eigvector
#-----------------------------------------------------------------
def det_multi(eigval):
# determines multiplicity of eigenvalue 0
	m=0
	for n in xrange(0,3):
		if abs(eigval[n])<=10.**(-7):
			m=m+1
	return m
#-----------------------------------------------------------------
def comp_pol(alpha, beta, gamma, epsilon, w):
# computes electric and magnetic polarizations
	P=np.zeros((3,4),'complex')
	Q=np.zeros((3,4),'complex')
	Q1=np.zeros((3,4),'complex')
	multi=np.zeros((4),'int')

	for k in xrange(0,4):

		eigval,eigvector=comp_eig(alpha, beta, gamma[k], epsilon, w)
		multi[k]=det_multi(eigval)
		if multi[k] ==1:
			for i in xrange(0,3):
				if abs(eigval[i])<=10.**(-7):
					P[:,k]=(eigvector[:,i])
		else:
			P[:,0]=[1.,0.,0.]
			P[:,1]=[1.,0.,0.]
			P[:,2]=[0.,1./math.sqrt(abs(beta)**2+abs(gamma[2])**2)*gamma[2],-1./math.sqrt(abs(beta)**2+abs(gamma[2])**2)*beta]
			P[:,3]=[0.,1./math.sqrt(abs(beta)**2+abs(gamma[3])**2)*gamma[3],-1./math.sqrt(abs(beta)**2+abs(gamma[3])**2)*beta]
			if multi[k] != 2:
				print 'Possible Problem with Scaling of Frequency!'
			break
	for l in xrange(0,4):
		k=[alpha, beta, gamma[l]]
		
		Q1[:,l]=np.cross(k,P[:,l])
		Q[:,l]=1./(w*mu)*Q1[:,l]

	return P,Q

#-----------------------------------------------------------------		
def const_matrix(alpha, beta, epsilon, w,t):
#constructs the T-matrix for a given layer
	P=np.zeros((4,4),'complex')
	D=np.zeros((4,4),'complex')
	T=np.zeros((4,4),'complex')
	gamma=quartic_solve(alpha,beta, epsilon,w)
	p,q=comp_pol(alpha, beta, gamma, epsilon, w)
	for i in xrange(0,4):
		for j in xrange(0,4):
			if i==0:
				D[i,j]=p[0,j]
			elif i==2:
				D[i,j]=p[1,j]
			elif i==1:
				D[i,j]=q[0,j]
			elif i==3:
				D[i,j]=q[1,j]
		P[i,i]=cmath.exp(-1.j*gamma[i]*t)
	D_inv=np.linalg.inv(D)
	T=np.dot(np.dot(D,P),D_inv)
	return T

#-----------------------------------------------------------------		

def const_matrix2(alpha, beta, epsilon1, epsilon2, w):
#constructs the D- and D_inv-matrix for vacuum (epsilon1) and substrat (epsilon2)
	D1=np.zeros((4,4),'complex')
	D2=np.zeros((4,4),'complex')
	gamma1=quartic_solve(alpha,beta, epsilon1,w)
	gamma2=quartic_solve(alpha,beta, epsilon2,w)
	p1,q1=comp_pol(alpha, beta, gamma1, epsilon1, w)
	p2,q2=comp_pol(alpha, beta, gamma2, epsilon2, w)
	for i in xrange(0,4):
		for j in xrange(0,4):
			if i==0:
				D1[i,j]=p1[0,j]
				D2[i,j]=p2[0,j]
			elif i==2:
				D1[i,j]=p1[1,j]
				D2[i,j]=p2[1,j]
			elif i==1:
				D1[i,j]=q1[0,j]
				D2[i,j]=q2[0,j]
			elif i==3:
				D1[i,j]=q1[1,j]
				D2[i,j]=q2[1,j]
		
	D_inv=np.linalg.inv(D1)
	
	return D_inv, D2, p2, q2

#-----------------------------------------------------------------		

def T_mult(T,D_inv,D):
#obtain T_ges by multiplication of layer transfer matrices
	T_ges=np.zeros((4,4),'complex')
	for i in xrange(0,4):
		T_ges[i,i]=1.
	for j in xrange(0,len(T)):
		T_ges=np.dot(T_ges,T[j])
	T_ges=np.dot(np.dot(D_inv,T_ges),D)
	return T_ges
#-----------------------------------------------------------------		
def amplitude(A0,T_ges,p,q):
#calculates the Poynting vector in the last layer and obtains R,T
#and the absorbance a=-log(T)
	A0_out=np.zeros((2,1),'complex')
	As_out=np.zeros((2,1),'complex')
	T1=np.array([[T_ges[0,0], T_ges[0,2]],[T_ges[2,0], T_ges[2,2]]])
	T2=np.array([[T_ges[1,0], T_ges[1,2]],[T_ges[3,0], T_ges[3,2]]])
	T1_inv=np.linalg.inv(T1)
	As_out=np.dot(T1_inv,A0)
	A0_out=np.dot(np.dot(T2,T1_inv),A0)
	#Algorithm for anisotropic substrate: Poynting vector can't be divided in two independent contributions
	#E=As_out[0]*p[:,0]+As_out[1]*p[:,2]
	#H=As_out[0]*q[:,0]+As_out[1]*q[:,2]
	#S=np.cross(E,H)
	#T=c*mu*(abs(S[0])+abs(S[1])+abs(S[2]))
	
	#For isotropic substrate: Two independent poynting vectors
	E1=E=As_out[0]*p[:,0]
	E2=As_out[1]*p[:,2]
	H1=As_out[0]*q[:,0]
	H2=As_out[1]*q[:,2]
	S1=np.cross(E1,H1)
	S2=np.cross(E2,H2)
	T1=c*mu*(abs(S1[0])+abs(S1[1])+abs(S1[2]))
	T2=c*mu*(abs(S2[0])+abs(S2[1])+abs(S2[2]))
	R1=abs(A0_out[0])**2
	R2=abs(A0_out[1])**2
	a1=0.
	a2=0.
	return R1,R2,T1,T2,a1,a2
#-----------------------------------------------------------------
def scale(w,t):
# Scales angular frequency and thicknesses in order to secure accurancy for
# different frequencies and angles
	for i in xrange(0,len(w)):
		if w[i] != 0:
			a=w[0]*10.**(-8) 
			for j in xrange(0,len(w)):
				w[j]=w[j]/a
			for k in xrange(0,len(t)):
				t[k]=a*t[k]
			break
	return w,t,a


#-------------------FUNCTIONS TO READ-----------------
#-------------------------INPUT-----------------------
def read(N):
#returns list of frequencies and corresponding dielectric functions
#files have to have the form N_xy.OUT for layer N
	inlist=[]
	WORKDIR=os.getcwd()
	files=os.listdir(WORKDIR)
	for i in xrange(0,len(files)):
		if files[i].split('_')[0]==str(N):
			inlist.append(files[i])

	# read all frequency lists and elements of the tensors into dictionary
	scheme=str(N)+'_'
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
	for i in xrange(0,len(data['w00'])):
		inter=np.zeros((3,3),'complex')
		data['w00'][i]=data['w00'][i]/hbar
		for k in xrange(0,3):
			for l in xrange(0,3):
				if len(data[str(k)+str(l)]) != 0:
					inter[k,l]=data[str(k)+str(l)][i]
		epsilon.append(inter)
	return epsilon, data['w00']
#-----------------------------------------------------------------
def read_vac():
#form the vacuum dielectric tensor for each frequency point of the provided 
#dielectric tensors
	epsilon,w=read(1)
	epsilon=[]
	inter=np.zeros((3,3),'complex')
	for i in xrange(0,3):
		inter[i,i]=1.
	for i in xrange(0,len(w)):
		epsilon.append(inter)
	return epsilon,w
#-----------------------------------------------------------------
		
def read_sub(w):
#form the vacuum dielectric tensor for each frequency point of the provided 
#dielectric tensors
	
	epsilon=[]
	inter=np.zeros((3,3),'complex')
	for i in xrange(0,3):
		inter[i,i]=11.8336
	for i in xrange(0,len(w)):
		epsilon.append(inter)
	return epsilon
					
#-----------------------------------------------------------------
def io():
# interactive function to obtain polarization and beam angle and total number of
# layers and their respective thickness (individual scales can be provided)
	print '\n     There are two options for the angle of the incoming beam:'\
      '\n     1 Choose a fixed angle                               '\
      '\n     2 set a range of numbers by typing "RANGE"       '

	alpha=raw_input('Give the angle of the incoming beam:')

	if alpha=='RANGE':
		print '\n Please type the two values between which the angle '\
			  '\n should be varied and the number of steps.                                  '
		alpha=raw_input('>>>> input the boundaries:')
		if len(alpha.split(" ")) > 3 or len(alpha.split(" ")) < 3 :
			sys.exit('Please provide only two numbers!')
		alpha0=float(alpha.split(" ")[0])*math.pi/180.
		alpha1=float(alpha.split(" ")[1])*math.pi/180.
		N=int(alpha.split(" ")[2])
		alpha=[]
		for i in xrange(0,N):
			alpha.append(alpha0+i*(alpha1-alpha0)/(N-1))
		print '\n Choose the polarization angle.'\
			  '\n For an array of beam angles the polarization'\
			  '\n has to be fixed'
		beta=raw_input('Give the polarization angle:')
		beta=float(beta)*math.pi/180.
	else:
		alpha=float(alpha)*math.pi/180.
		print '\n     There are two options for the polarization angle:'\
			  '\n     1 Choose a fixed angle                           '\
		      '\n     2 set a range of numbers by typing "RANGE"       '
		beta=raw_input('Give the polarization angle:')
		if beta=='RANGE':
			print '\n Please type the two values between which the angle '\
				  '\n should be varied and the number of steps.                                  '
			beta=raw_input('>>>> input the boundaries:')
			if len(beta.split(" ")) > 3 or len(beta.split(" ")) < 3 :
				sys.exit('Please provide only two numbers!')
			beta0=float(beta.split(" ")[0])*math.pi/180.
			beta1=float(beta.split(" ")[1])*math.pi/180.
			N=int(beta.split(" ")[2])
			beta=[]
			for i in xrange(0,N):
				beta.append(beta0+i*(beta1-beta0)/(N-1))
		else:
			beta=float(beta)*math.pi/180.
	print 'Provide the number of layers, not counting vacuum and substrat'
	N=input('>>Number of layers:')
	t=[]
	print 'Please provide in the following the thickness of each layer.'\
		  'Your also asked to provide the unit of length for each layer.'\
		  'Choose between: m, mm, nm'
	for i in xrange(0,N):
		print '*******************************************'
		print 'For Layer ', i+1
		a=raw_input('>>>>Unit:')
		if a =='m':
			b=1.
		elif a=='mm':
			b=10.**(-3)
		elif a=='nm':
			b=10.**(-9) 
		t.append(float(input('>>>>Thickness:'))*b)
	return alpha, beta, t, N

#-------------------------MAIN--PROGRAM----------------------------
		

#time_start=time.clock()

#sigma, beta0,t, N=io() #angle of polarization, angle of incoming beam, thickness of layers, Number of layers
alpha=0.
sigma=0.
beta0=0.

t=[10.*10.**(-9)]
N=1            
EPS=[]
for i in xrange(0,N+2):
	if i==0:
		w=read_vac()[1]
		EPS.append(read_vac()[0])
	elif i==N+1:
		EPS.append(read_sub(w))
	else:
		EPS.append(read(i)[0])

w,t,b=scale(w,t)
#                         1.POLARIZATION DEPENDENCY
if hasattr(sigma,'__len__')==True and hasattr(beta0,'__len__')==False:
	for k in xrange(0,len(sigma)):
		A0=[math.sin(sigma[k]),math.cos(sigma[k])]
		with open('reflection'+str(k)+'.out','w') as f, open('transmission'+str(k)+'.out','w') as g, open('absorbance'+str(k)+'.out','w') as h:
			for i in xrange(0,len(w)):
				beta=w[i]/c*math.sin(beta0)
				R_cof=[]
				T_cof=[]
				a_cof=[]		
				T=[]
				T_ges=[]
				for j in xrange(1,(len(EPS))-1):  #loop over layers
					T.append(const_matrix(alpha, beta, EPS[j][i], w[i],t[j-1]))
				D_inv,D,p,q=const_matrix2(alpha, beta, EPS[0][i], EPS[-1][i], w[i])
				T_ges=T_mult(T,D_inv,D)
				R_cof.append(amplitude(A0,T_ges,p,q)[0])
				R_cof.append(amplitude(A0,T_ges,p,q)[1])
				T_cof.append(amplitude(A0,T_ges,p,q)[2])
				T_cof.append(amplitude(A0,T_ges,p,q)[3])
				a_cof.append(amplitude(A0,T_ges,p,q)[4])
				a_cof.append(amplitude(A0,T_ges,p,q)[5])
				f.write('%g  %1.9e %1.9e\n' % (w[i]*b, R_cof[0], R_cof[1]))
				g.write('%g  %1.9e %1.9e\n' % (w[i]*b, T_cof[0], T_cof[1]))
				h.write('%g  %4.9e %4.9e\n' % (w[i]*b, a_cof[0], a_cof[1]))

#                         1.ANGULAR DEPENDENCY
elif hasattr(sigma,'__len__')==False and hasattr(beta0,'__len__')==True and len(w)!=1:
	A0=[math.sin(sigma),math.cos(sigma)]
	for k in xrange(0,len(beta0)):	
		with open('reflection'+str(k)+'.out','w') as f, open('transmission'+str(k)+'.out','w') as g, open('absorbance'+str(k)+'.out','w') as h:
			for i in xrange(0,len(w)):
				beta=w[i]/c*math.sin(beta0[k])		
				R_cof=[]
				T_cof=[]
				a_cof=[]		
				T=[]
				T_ges=[]
				for j in xrange(1,(len(EPS))-1):  #loop over layers
					T.append(const_matrix(alpha, beta, EPS[j][i], w[i],t[j-1]))
				D_inv,D,p,q=const_matrix2(alpha, beta, EPS[0][i], EPS[-1][i], w[i])
				T_ges=T_mult(T,D_inv,D)
				R_cof.append(amplitude(A0,T_ges,p,q)[0])
				R_cof.append(amplitude(A0,T_ges,p,q)[1])
				T_cof.append(amplitude(A0,T_ges,p,q)[2])
				T_cof.append(amplitude(A0,T_ges,p,q)[3])
				a_cof.append(amplitude(A0,T_ges,p,q)[4])
				a_cof.append(amplitude(A0,T_ges,p,q)[5])
				f.write('%g  %1.9e %1.9e\n' % (w[i]*b*hbar, R_cof[0], R_cof[1]))
				g.write('%g  %1.9e %1.9e\n' % (w[i]*b*hbar, T_cof[0], T_cof[1]))
				h.write('%g  %4.9e %4.9e\n' % (w[i]*b*hbar, a_cof[0], a_cof[1]))
#                         1.NEITHER ANGULAR NOR POLARIZATION DEPENDENCY
elif hasattr(sigma,'__len__')==False and hasattr(beta0,'__len__')==False and len(w)!=1:
	with open('reflection.out','w') as f, open('transmission.out','w') as g, open('absorbance.out','w') as h:
		#A0=[math.sin(sigma),math.cos(sigma)]
		A0=[1.,0.]
		
		for i in xrange(0,len(w)):
			beta=w[i]/c*math.sin(beta0)
			R_cof=[]
			T_cof=[]
			a_cof=[]
			T=[]
			T_ges=[]
			for j in xrange(1,(len(EPS))-1):  #loop over layers
				T.append(const_matrix(alpha, beta, EPS[j][i], w[i],t[j-1]))
			D_inv,D,p,q=const_matrix2(alpha, beta, EPS[0][i], EPS[-1][i], w[i])

			T_ges=T_mult(T,D_inv,D)
			R_cof.append(amplitude(A0,T_ges,p,q)[0])
			R_cof.append(amplitude(A0,T_ges,p,q)[1])
			T_cof.append(amplitude(A0,T_ges,p,q)[2])
			T_cof.append(amplitude(A0,T_ges,p,q)[3])
			a_cof.append(amplitude(A0,T_ges,p,q)[4])
			a_cof.append(amplitude(A0,T_ges,p,q)[5])
			f.write('%g  %1.9e %1.9e\n' % (w[i]*b*hbar, R_cof[0], R_cof[1]))
			g.write('%g  %1.9e %1.9e\n' % (w[i]*b*hbar, T_cof[0], T_cof[1]))
			h.write('%g  %4.9e %4.9e\n' % (w[i]*b*hbar, a_cof[0], a_cof[1]))
#****************************************************************************************
#                         HARD CODED SPECIAL CASES
#                         1.ANGULAR DEPENDENCE AT ONE FREQUENCY
elif hasattr(sigma,'__len__')==False and hasattr(beta0,'__len__')==True and len(w)==1:
	A0=[math.sin(sigma),math.cos(sigma)]
	with open('reflection.out','w') as f, open('transmission.out','w') as g, open('absorbance.out','w') as h:
		for k in xrange(0,len(beta0)):	
			beta=w[0]/c*math.sin(beta0[k])		
			R_cof=[]
			T_cof=[]
			a_cof=[]		
			T=[]
			T_ges=[]
			for j in xrange(1,(len(EPS))-1):  #loop over layers
				T.append(const_matrix(alpha, beta, EPS[j][0], w[0],t[j-1]))
			D_inv,D,p,q=const_matrix2(alpha, beta, EPS[0][0], EPS[-1][0], w[0])
			T_ges=T_mult(T,D_inv,D)
			R_cof.append(amplitude(A0,T_ges,p,q)[0])
			R_cof.append(amplitude(A0,T_ges,p,q)[1])
			T_cof.append(amplitude(A0,T_ges,p,q)[2])
			T_cof.append(amplitude(A0,T_ges,p,q)[3])
			a_cof.append(amplitude(A0,T_ges,p,q)[4])
			a_cof.append(amplitude(A0,T_ges,p,q)[5])
			f.write('%g  %1.9e %1.9e\n' % (beta0[k]*180./math.pi, R_cof[0], R_cof[1]))
			g.write('%g  %1.9e %1.9e\n' % (beta0[k]*180./math.pi, T_cof[0], T_cof[1]))
			h.write('%g  %4.9e %4.9e\n' % (beta0[k]*180./math.pi, a_cof[0], a_cof[1]))
			
#					      2.POLARIZATION AT ONE FREQUENCY
elif hasattr(sigma,'__len__')==True and hasattr(beta0,'__len__')==False and len(w)==1:
	with open('reflection.out','w') as f, open('transmission.out','w') as g, open('absorbance.out','w') as h:
		for k in xrange(0,len(sigma)):
			A0=[math.sin(sigma[k]),math.cos(sigma[k])]
			beta=w[0]/c*math.sin(beta0)
			R_cof=[]
			T_cof=[]
			a_cof=[]		
			T=[]
			T_ges=[]
			for j in xrange(1,(len(EPS))-1):  #loop over layers
				T.append(const_matrix(alpha, beta, EPS[j][0], w[0],t[j-1]))
			D_inv,D,p,q=const_matrix2(alpha, beta, EPS[0][0], EPS[-1][0], w[0])
			T_ges=T_mult(T,D_inv,D)
			R_cof.append(amplitude(A0,T_ges,p,q)[0])
			R_cof.append(amplitude(A0,T_ges,p,q)[1])
			T_cof.append(amplitude(A0,T_ges,p,q)[2])
			T_cof.append(amplitude(A0,T_ges,p,q)[3])
			a_cof.append(amplitude(A0,T_ges,p,q)[4])
			a_cof.append(amplitude(A0,T_ges,p,q)[5])
			f.write('%g  %1.9e %1.9e\n' % (sigma[k]*180./math.pi, R_cof[0], R_cof[1]))
			g.write('%g  %1.9e %1.9e\n' % (sigma[k]*180./math.pi, T_cof[0], T_cof[1]))
			h.write('%g  %4.9e %4.9e\n' % (sigma[k]*180./math.pi, a_cof[0], a_cof[1]))

#***************************************************************************************************
else:
	print 'Incorrect Input: You cannot vary Polarization and Beam Angle at the same time.'
