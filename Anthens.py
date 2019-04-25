# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
import scipy as sc
import math as math
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
        g=np.round(sc.optimize.newton(f, i),5)
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
F=F/max(F)
F0=F0/max(F0)

for i in range(len(F)):
    F[i]=20*math.log10(F[i])
    F0[i]=20*math.log10(F0[i])
    
'''
pl.figure(1)
pl.plot(c,F0,label=u'n=0')
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,20)
#pl.savefig('Anthens1.png',dpi=300)

pl.figure(2)
pl.plot(c,F,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,10)
#pl.savefig('Anthens2.png',dpi=300)'''

def G2(m,n): #Función F
    if m==0:
        return 1
    else:
        h1=1
        for i in np.arange(1,n):
            h1*=(1-((zeros[m]**2*(A**2+(n-0.5)**2))/(zeros[n]**2*(A**2+(i-0.5)**2))))
        h2=1
        for i in np.arange(1,n):
            if m!=i:
                h2*=(1-(zeros[m]/zeros[i])**2)
        return -f0(np.pi*zeros[m])*h1/h2

def g0(x,n): #Aplicado L'Hôpital
    h=0
    for i in range(n):
        h+=((G2(i,n)*f0(x*zeros[i]))/(f0(np.pi*zeros[i])**2))
    return 4*h/(np.pi**2)

g=np.abs(g0(c,n))
g=g/max(g)

'''
pl.figure(3)
pl.plot(c,g,label=u'Abertura')
pl.legend(loc='upper left')
pl.xlim(0,np.pi)
pl.ylim(0,1)
#pl.savefig('Anthens3.png',dpi=300)'''


matrix=np.loadtxt('datos.txt')
N=matrix[0]; I2=matrix[1]


lamb=2.12 #npi
a=5*lamb; d=a/10 #Distancia uniforme entre círculos
kt=2*np.pi/lamb


M=len(N); ka=int(max(N)) #Para meter de novo simetría facemos ka=int(max(N)/4). Ollo liña 147.
rho=np.zeros([M,ka]); beta=np.zeros([M,ka])
        

for i in range(M): #Páxina 249
    for j in range(ka):
        rho[i,j]=d*np.sqrt((i+0.5)**2+(j+0.5)**2) #Podo eliminar unha dimensión
        beta[i,j]=(2*j+1)*np.pi/N[i]


I=g0(rho*np.pi/a,n)

def F(theta,phi,n): #Páxina 253
    aux=0
    for m in range(M):
        for n in range(int(N[m])):                
            aux+=I[m,0] * np.cos(kt*rho[m,0]*np.cos(beta[m,n])*np.sin(theta)*np.cos(phi)) * np.cos(kt*rho[m,0]*np.sin(beta[m,n])*np.sin(theta)*np.sin(phi))       
    return aux


F1=np.abs(F(c,0,n))
F1=F1/max(F1)

for i in range(len(F1)):
    F1[i]=20*math.log10(F1[i])

#Esta interesa 
'''
pl.figure(5) #Páxina 254
pl.plot(c,F1,label=u'-')
pl.xlabel('$\Theta$ (rad)')
pl.legend(loc='upper right')
pl.xlim(0,np.pi/2)
pl.ylim(-50,0)
#pl.savefig('---.png',dpi=300)'''

#Caso 1 ---------------------------------------------------------------------

u=[0,1.4839, 1.8933, 2.9268, 3.9622, 5.0416]
caso1=np.abs(G1(c,n)); caso1=caso1/max(caso1)


#Caso 2 ---------------------------------------------------------------------
    
u=[0,0.6322, 1.9308, 3.7674, 4.3929, 5.2633]
caso2=np.abs(G1(c,n)); caso2=caso2/max(caso2)


#Caso 3 ---------------------------------------------------------------------
    
u=[0, 0.5967, 1.7837, 3.6420, 4.3039, 5.2129]
v=[0,0.5225, 0.5268, 0, 0, 0]

def Fim(x,n): #Traballando con imaxinarios
    h1=1
    for j in np.arange(1,n):
        h1*=((u[j]**4+v[j]**4+x**4+2*u[j]**2*v[j]**2+2*x**2*(v[j]**2-u[j]**2))/((u[j]**2+v[j]**2)**2))
    h2=1
    for j in np.arange(1,n):
        h2*=1-(x/zeros[j])**2
    return (f(x)/(h2*np.pi*x))**2*h1

caso3=np.abs(Fim(c,n)); caso3=caso3/max(caso3)

for i in range(len(caso3)):
    caso1[i]=20*math.log10(caso1[i])
    caso2[i]=20*math.log10(caso2[i])
    caso3[i]=10*math.log10(caso3[i])

#Estas tres interesan
'''
pl.figure(6) #Mesmo do principio, modificando diferente a posición dos raíces
pl.plot(c,caso1,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,20)
pl.xlabel('Primeiro caso')
pl.savefig('Caso1.png',dpi=300)

pl.figure(7) #Segundo caso
pl.plot(c,caso2,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,20)
pl.xlabel('Segundo caso')
pl.savefig('Caso2.png',dpi=300)

pl.figure(7) #Terceiro caso
pl.plot(c,caso3,label=u'n=%i'%n)
pl.legend(loc='upper right')
pl.ylim(-50,0)
pl.xlim(0,10)
pl.xlabel('Terceiro caso')
#pl.savefig('Caso3.png',dpi=300)'''

'''
x=np.arange(-np.pi/2,np.pi/2,0.01); x=np.delete(x,np.where(x==0))
y=np.arange(-np.pi/2,np.pi/2,0.01); y=np.delete(x,np.where(x==0))
X,Y=np.meshgrid(x,y)
        
phi=np.abs(np.arctan(X/Y))
r=np.sqrt(X**2+Y**2)

z=np.abs(F(r,phi,n)) #Libro

nf,nc=np.shape(z)
Z=np.zeros([nf,nc])

for i in range(nf):
    for j in range(nc):
        Z[i,j]=20*math.log10(z[i,j])

Z=Z-np.max(Z)

#Representación 3D

fig = pl.figure()
ax = pl.axes(projection='3d')
ax.contour3D(X, Y, Z, 100, cmap='viridis')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')'''

#ax.view_init(0, 0)
#pl.savefig('3D.png',dpi=500)




print 'Tempo de execución:',time.clock()-t0,'(s)'
