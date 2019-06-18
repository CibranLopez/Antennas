# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import newton
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D
import time as time

t0=time.clock()

#Parámetros
s=float(-15)
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

def G1(x,n):
    h1=1
    for j in np.arange(1,n):
        h1*=1-(x/u[j])**2
    h2=1
    for j in np.arange(1,n):
        h2*=1-(x/zeros[j])**2
    return (f(x)*h1)/(h2*np.pi*x)

F=np.abs(G1(c,n))
F0=np.abs(G1(c,0))

F=20*np.log10(F/max(F))
F0=20*np.log10(F0/max(F0))
    
'''
pl.figure(1)
pl.plot(c,F0,label=u'Apertura uniforme')
pl.legend(loc='upper right')
pl.xlabel(r'$u=sin(\theta) 2a/\lambda$')
pl.ylabel('dB')
pl.ylim(-50,0)
pl.xlim(0,20)

pl.figure(2)
pl.plot(c,F,label=u'Apertura uniforme')
pl.legend(loc='upper right')
pl.xlabel(r'$u=sin(\theta) 2a/\lambda$')
pl.ylabel('dB')
pl.ylim(-50,0)
pl.xlim(0,10)'''

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

g=np.abs(g0(c,n))
g=g/max(g)

'''
pl.figure(3)
pl.plot(c,g,label=u'Amplitude relativa')
pl.legend(loc='upper left')
pl.xlabel(r'$u=sin(\theta) 2a/\lambda$')
pl.xlabel(r'$p=\rho \pi /a$')
pl.ylabel('Amplitude relativa')
pl.xlim(0,np.pi)
pl.ylim(0,1)'''

#Páxina 254 -----------------------------------------------------------------

N=[4,8,16,20,28,32,40,44,52,56]

d=0.5 #Distancia uniforme entre círculos
M=len(N); ka=int(max(N)/4)
rho=np.zeros(M); beta=np.zeros([M,ka])

for i in range(M):
    for j in range(ka):
        rho[i]=d*(i+0.5)
        beta[i,j]=(2*j+1)*np.pi/N[i]

I=np.absolute(g0(rho*np.pi/5,n))

def F(theta,phi,n): #Páxina 253
    aux=0
    for m in range(M):
        for n in range(int(N[m])/4):                
            aux+=I[m] * np.cos(2*np.pi*rho[m]*np.cos(beta[m,n])*np.sin(theta)*np.cos(phi)) * np.cos(2*np.pi*rho[m]*np.sin(beta[m,n])*np.sin(theta)*np.sin(phi))       
    return 4*aux

F1=np.abs(F(c,0,n))
F1=20*np.log10(F1/max(F1))

'''
pl.figure(4) #Páxina 254
pl.plot(c,F1,label=u'Discretización')
pl.xlabel(r'$\theta$')
pl.ylabel('dB')
pl.legend(loc='upper right')
pl.xlim(0,np.pi/2)
pl.ylim(-50,0)'''


#Caso 1 ---------------------------------------------------------------------

u=[0,1.4839, 1.8933, 2.9268, 3.9622, 5.0416]
caso1=np.abs(G1(c,n)); caso1=caso1/max(caso1)


#Caso 2 ---------------------------------------------------------------------
    
u=[0,0.6322, 1.9308, 3.7674, 4.3929, 5.2633]
caso2=np.abs(G1(c,n)); caso2=caso2/max(caso2)


#Caso 3 ---------------------------------------------------------------------

u=[0, 0.5967, 1.7837, 3.6420, 4.3039, 5.2129]
v=[0, 0.5225, 0.5268, 0, 0, 0]

def Fim(x,n):
    h1=1
    for j in np.arange(1,n):
        h1*=((u[j]**4+v[j]**4+x**4+2*u[j]**2*v[j]**2+2*x**2*(v[j]**2-u[j]**2))/((u[j]**2+v[j]**2)**2))
    h2=1
    for j in np.arange(1,n):
        h2*=1-(x/zeros[j])**2
    return (f(x)/(h2*np.pi*x))**2*h1

caso3=np.abs(Fim(c,n)); caso3=caso3/max(caso3)

for i in range(len(u)):
    u[i]=u[i]+1j*v[i]
    
distrbn=g0(c,n)

phi=np.angle(distrbn)

distrbn=np.absolute(distrbn)
distrbn=distrbn/np.max(distrbn)

caso1=20*np.log10(caso1)
caso2=20*np.log10(caso2)
caso3=10*np.log10(caso3)

'''
pl.figure(5) #Mesmo do principio, modificando diferente a posición dos raíces
pl.plot(c,caso1,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,20)
pl.xlabel('Primeiro caso')

pl.figure(6) #Segundo caso
pl.plot(c,caso2,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,20)
pl.xlabel('Segundo caso')

pl.figure(7) #Terceiro caso
pl.plot(c,caso3,label=u'Raíces complexas')
pl.legend(loc='upper right')
pl.xlabel('u')
pl.ylabel('dB')
pl.ylim(-50,0)
pl.xlim(0,10)

pl.figure(8)
pl.plot(c,distrbn,label=u'Amplitude')
pl.ylabel('Amplitude relativa')
pl.xlabel('u')
pl.legend(loc='upper right')
pl.xlim(0,np.pi)

pl.figure(9)
pl.plot(c,phi,label=u'Fase')
pl.ylabel('Fase relativa')
pl.xlabel('u')
pl.legend(loc='upper right')
pl.xlim(0,np.pi)'''


uc=[0+1j*0, 0.5967+1j*0.5225, 1.7837+1j*0.5268, 3.6420+1j*0, 4.3039+1j*0, 5.2129+1j*0]
l=1; M=20; a=5 #Núemero de puntos e discos

N=np.arange(1,M+1,1)*4*l
p=np.arange(1,M+1,1)*np.pi/M


#Expansión 6.99 -------------------------------------------------------------

'''
I=g0(p,n)
I2=a*np.pi*I/(10*l) 

def Fcomplex(theta,phi,n): #Páxina 259
    aux=0+1j*0
    for m in range(M):
        aux+=N[m]*I2[m]*sc.special.j0(2*a*np.sin(theta)*p[m])
        for s in np.arange(1,20,1):
            aux+=2*(-1)**s*N[m]*I2[m]*sc.special.jn(s*N[m],2*a*np.sin(theta)*p[m])*np.cos(s*N[m]*phi)
    return aux


F1complex=np.absolute(Fcomplex(c,0,n))

F1complex=20*np.log10(F1complex/np.max(F1complex))


pl.figure(10) 
pl.plot(c,F1complex,label=u'Fcomplex')
pl.xlim(0,np.pi/2)
pl.ylim(-50,0)
pl.legend(loc='upper right')
#pl.savefig('---.png',dpi=300)'''


#Desenvolvemento 6.99 -------------------------------------------------------

I=g0(p,n)
I2=a*np.pi*I/(10*l) 

beta=np.zeros([M,np.max(N)])

M2=15; tramos=4  #Discos sen modificar



for i in range(M2):
    aux=0
    for j in range(N[i]):
        beta[i,j]=aux
        aux+=2*np.pi/N[i]
      
for i in np.arange(M2,M,1):
    aux=0
    for t in range(tramos):
        for j in range(N[i]/(2*tramos)):
            beta[i,j]=2*np.pi*t/tramos+aux
            aux+=np.pi/N[i]
'''
for i in range(M):
    aux=0
    if i<(M2-1):
        for j in range(N[i]):
            beta[i,j]=aux
            aux+=2*np.pi/N[i]
    else:
        for t in range(tramos):
            for j in range(N[i]/(2*tramos)):
                beta[i,j]=2*np.pi*t/tramos+aux
                aux+=np.pi/N[i]'''



def Fcomplex_2(theta,phi,n):
    aux=0+1j*0
    for m in range(M):
        for n in range(N[m]):
            for q in np.arange(-20,21,1):
                aux+=(1j)**q*I2[m]*sc.special.jn(q,2*a*np.sin(theta)*p[m])*np.e**(1j*q*(phi-beta[m,n]))
    return aux


pl.figure(10) 

F2complex=np.absolute(Fcomplex_2(c,0,n))
F2complex=20*np.log10(F2complex/np.max(F2complex))
pl.plot(c,F2complex,label=u'Fcomplex')
F3complex=np.absolute(Fcomplex_2(c,np.pi/3,n))
F3complex=20*np.log10(F3complex/np.max(F3complex))
pl.plot(c,F3complex,label=u'Fcomplex')
F4complex=np.absolute(Fcomplex_2(c,2*np.pi/3,n))
F4complex=20*np.log10(F4complex/np.max(F4complex))
pl.plot(c,F4complex,label=u'Fcomplex')

pl.xlim(0,np.pi/2)
pl.ylim(-50,0)
pl.legend(loc='upper right')


#Paso a 3D ------------------------------------------------------------------
'''
x=np.arange(-np.pi/2,np.pi/2,0.01); x=np.delete(x,np.where(x==0))
X,Y=np.meshgrid(x,x)

phi=np.zeros(np.shape(X))
r=np.sqrt(X**2+Y**2)
nf,nc=np.shape(X)

#phi=beta

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


Z=Fcomplex_2(r,phi,n)
           
Z=np.absolute(Z)
Z=20*np.log10(Z/np.max(Z))

for i in range(nf):
    for j in range(nc):
        if r[i,j]>np.pi/2 or Z[i,j]<-50:
            Z[i,j]=-50


fig = pl.figure(11)
pl.imshow(Z)


fig = pl.figure(12)
ax = pl.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.grid(False)
ax.set_visible(True)
ax.set_zlim(-50,0)
#pl.savefig('Proba2.png',dpi=300)
#ax.view_init(0, 0)'''




print 'Tempo de execución:',time.clock()-t0,'(s)'

