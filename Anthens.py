# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import newton
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D
import time as time

t0=time.clock()

#Parámetros
s=float(-25)
n=6

print 's=',s,'; n=',n

c=np.arange(0,20,0.01); c=np.delete(c,0)
A=np.arccosh(10**(-s/20))/np.pi
zeros=np.zeros(n+1); u=np.zeros(n+1)

def f(x):
    return sc.special.j1(np.pi*x)

def f0(x):
    return sc.special.j0(x)

i=0;g=0
while zeros[-1]==0:
    try:
        g=np.round(newton(f, i),5)
        if (g in zeros)==False:
            if g>=0:
                zeros[i]=g
    except RuntimeError:
        pass
    i+=1

for i in range(n+1): 
    u[i]=zeros[n]*np.sqrt((A**2+(i-0.5)**2)/(A**2+(n-0.5)**2))    

def G2(m,n):
    if m==0:
        return 1
    else:
        h1=1
        for i in np.arange(1,n):
            h1*=(1-(zeros[m]/u[i])**2)
        h2=1
        for i in np.arange(1,n):
            if m!=i:
                h2*=(1-(zeros[m]/zeros[i])**2)
        return -f0(np.pi*zeros[m])*h1/h2

def g0(x,n):
    aux=0+1j*0
    for i in range(n):
        aux+=((G2(i,n)*f0(x*zeros[i]))/(f0(np.pi*zeros[i])**2))
    return 2*aux/(np.pi**2)

'''
g=np.abs(g0(c,n))
g=g/max(g)

pl.figure(1)
pl.plot(c,g,label=u'Amplitude relativa')
pl.legend(loc='upper left')
pl.xlim(0,np.pi)
pl.ylim(0,1)'''


u=[0+1j*0, 0.5967+1j*0.5225, 1.7837+1j*0.5268, 3.6420+1j*0, 4.3039+1j*0, 5.2129+1j*0]
'''
#Ecuación do raio: rmax=a*(1+cos(beta))/2 beta=[0,np.pi/2]

M=20; a=5; d=0.25 #Distancia ente elementos

rmax=[1*a,0.99*a,0.96*a,0.92*a,0.85*a,0.78*a,0.69*a,0.6*a,0.5*a,0.6*a,0.69*a,0.78*a,0.85*a,0.92*a,0.96*a,0.99*a,a,0.99*a,0.96*a,0.92*a,0.85*a,0.78*a,0.69*a,0.6*a,0.5*a,0.6*a,0.69*a,0.78*a,0.85*a,0.92*a,0.96*a,0.99*a,a,0.99*a,0.96*a,0.92*a,0.85*a,0.78*a,0.69*a,0.6*a,0.5*a,0.6*a,0.69*a,0.78*a,0.85*a,0.92*a,0.96*a,0.99*a,a,0.99*a,0.96*a,0.92*a,0.85*a,0.78*a,0.69*a,0.6*a,0.5*a,0.6*a,0.69*a,0.78*a,0.85*a,0.92*a,0.96*a,0.99*a] 
#rmax=[a,0.96*a,0.69*a,0.5*a,0.69*a,0.96*a,a,0.96*a,0.69*a,0.5*a,0.69*a,0.96*a,a,0.96*a,0.69*a,0.5*a,0.69*a,0.96*a,a,0.96*a,0.69*a,0.5*a,0.69*a,0.96*a] 

lon=len(rmax) #Número de gajos

p=np.arange(1,M+1,1)*np.pi/M

I=np.real(g0(p,n))

beta=np.zeros([M,lon]); elementos=np.zeros(lon); I2=np.zeros(M); N=np.zeros(M)

def pos_sup_max(a):
    aux=np.max(a)
    for i in range(len(a)):
        if aux==a[i]:
            a[i]=np.min(a)-1
            return i

for n in range(lon):
    elementos[n]=np.int(rmax[n]/d)
    beta[:,n]=2*np.pi*pos_sup_max(rmax)/lon

for m in range(M):
    N[m]=np.count_nonzero(elementos) #N é azimutal, elementos radial
    elementos[np.where((m+1)==elementos)]=0
    I2[m]=2*a*p[m]*I[m]/N[m]

def F1(theta,phi,n):
    aux=0+1j*0
    for m in range(M):
        for n in range(int(N[m])):
            for q in np.arange(-20,21,1):
                aux+=(1j)**q*I2[m]*sc.special.jn(q,2*a*np.sin(theta)*p[m])*np.e**(1j*q*(phi-beta[m,n]))
    return aux


pl.figure(2) 

F2complex=np.absolute(F1(c,0,n))
F2complex=20*np.log10(F2complex/np.max(F2complex))
pl.plot(c,F2complex,label=u'0')

F3complex=np.absolute(F1(c,np.pi*4/3,n))
F3complex=20*np.log10(F3complex/np.max(F3complex))
pl.plot(c,F3complex,label=u'4pi/3')

F4complex=np.absolute(F1(c,np.pi*5/3,n))
F4complex=20*np.log10(F4complex/np.max(F4complex))
pl.plot(c,F4complex,label=u'5pi/3')

pl.xlim(0,np.pi/2)
pl.ylim(-50,0)
pl.legend(loc='upper right')
#pl.savefig('phi_cte.png',dpi=300)'''


N=[12,12,11,11,10,8,6,6,5,5,4,2]
M=12

p=np.zeros([M,M]); I=np.zeros([M,M])

for i in range(M):
    for j in range(N[i]):
        p[i,j]=np.sqrt((i+0.5)**2+(j+0.5)**2)*np.pi/M

I=np.real(g0(p,n))

def F2(theta,phi,n):
    aux=0+1j*0
    for i in range(M):
        for j in range(N[i]):
            aux+=I[i,j] * np.cos((i+0.5)*np.pi*np.sin(theta)*np.cos(phi)) * np.cos((j+0.5)*np.pi*np.sin(theta)*np.sin(phi))
    return 4*aux


f=np.absolute(F2(c,0,n))
f=20*np.log10(f/np.max(f))
pl.plot(c,f,label=u'0')

f2=np.absolute(F2(c,np.pi*4/3,n))
f2=20*np.log10(f2/np.max(f2))
pl.plot(c,f2,label=u'4/3')

f3=np.absolute(F2(c,np.pi*5/3,n))
f3=20*np.log10(f3/np.max(f3))
pl.plot(c,f3,label=u'5/3')

pl.xlim(0,np.pi/2)
pl.ylim(-50,0)
pl.legend(loc='upper right')
#pl.savefig('Cartesianas.png',dpi=300)'''



#Paso a 3D ------------------------------------------------------------------
'''
x=np.arange(-np.pi/2,np.pi/2,0.01); x=np.delete(x,np.where(x==0))
X,Y=np.meshgrid(x,x)

phi=np.zeros(np.shape(X))
r=np.sqrt(X**2+Y**2)
nf,nc=np.shape(X)


for i in range(nf):
    for j in range(nc):
        angle=np.abs(np.arctan(Y[i,j]/X[i,j]))
        
        if X[i,j]>0:
            if Y[i,j]>0:
                phi[i,j]=angle
            elif Y[i,j]<0:
                phi[i,j]=-angle
                
        elif X[i,j]<0:
            if Y[i,j]>0:
                phi[i,j]=np.pi-angle
            elif Y[i,j]<0:
                phi[i,j]=angle-np.pi

#for i in range(nf)  #Máis lento
#    for j in range(nc):
#        angle=np.abs(np.arctan(Y[i,j]/X[i,j]))
        
#        if X[i,j]>0:
#            phi[i,j]=np.sign(Y[i,j])*angle
                
#        elif X[i,j]<0:
#            phi[i,j]=np.sign(Y[i,j])*np.abs(np.pi-angle)


Z=F...(r,phi,n)
           
Z=np.absolute(Z)
Z=20*np.log10(Z/np.max(Z))

for i in range(nf):
    for j in range(nc):
        if r[i,j]>np.pi/2 or Z[i,j]<-50:
            Z[i,j]=-50


fig = pl.figure(3)
pl.imshow(Z)
pl.savefig('Corte_horizontal.png',dpi=300)


fig = pl.figure(4)
ax = pl.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.grid(False)
ax.set_visible(True)
ax.set_zlim(-50,0)
pl.savefig('Proba3d_2.png',dpi=300)
#ax.view_init(0, 0)'''




print 'Tempo de execución:',time.clock()-t0,'(s)'

