'''

                            Online Python Compiler.
                Code, Compile, Run and Debug python program online.
Write your code in this editor and press "Run" button to execute it.

'''

import numpy as np
import math as math
import csv

dx=np.pi/10
alpha=0.5
dt=0.0001
T=1.0
L=np.pi
x_size=int(L/dx)
t_size=int(T/dt)
g=list()
g.append(1)
g.append(-alpha)

def f(x,t):
    return ((0.5*math.gamma(3)/math.gamma(3.5))*(t**(1.5))+(math.gamma(2)/math.gamma(2.5))*(t**0.5))*np.sin(x)+(0.5*(t**2)+t)*np.sin(x)
for i in range(x_size):
    print(f(i*dx,1))
def analytic(t,x):
    return (0.5*(t**2)+t)*np.sin(x)

def gamma(j):
    if j>170:
       gam=((j+1)**(-alpha-1))/(math.gamma(-alpha))
    else:
        gam=(math.gamma(j-alpha))/(math.gamma(-alpha)*math.gamma(j+1))
    return gam

s=(t_size+1,x_size+1)
u=np.zeros(s)
u_analytic=np.zeros(x_size+1)

def eval1():
    for j in range(2,t_size+1):
        gam=g[j-1]*((j-alpha-1)/j)
        g.append(gam)
    return g
eval1()

for i in range(x_size+1):
    u[0,i]=0.0
sum = 0.0 ################
eps=0.001
complete=True
coef=dt**alpha
print(coef)

for i in range(t_size):
    u[i,0]=0.0
    u[i,x_size]=0.0
    for j in range(1,x_size):
        sum = 0.0 ############
        u[i+1,j]=coef*(((u[i,j+1]-2.0*u[i,j]+u[i,j-1])/(dx*dx))+f((j*dx),(i+1)*dt))
        #print(i,"asd",j, "fgh",u[i,j+1],"jkl",-2.0*u[i,j],"qwe",u[i,j-1],"rty",u[i+1,j],end=" ")
        for k in range(1,i+2):
            #if k<6: 
                u[i+1,j]-=u[i-k+1,j]*g[k]
        #        sum += u[i-k+1,j]*g[k] ##########
        #print("sum =", sum, "f * dt =", f(j * dx, i * dt) * pow(dt, alpha))
                
        if (abs(u[i+1,j]-u[i,j]))>eps:
            print("df",abs(u[i+1,j]-u[i,j]))
            iter=i+1
            complete=True
    if complete==True:
        continue  
    else:
        break

for i in range((x_size)):
    u_analytic[i]=analytic(iter*dt,i*dx)

#plt.plot(u[iter],label=f'numerical for {alpha}')
print(u[iter], iter, len(u[iter]))
#plt.plot(u_analytic,label='analytic')
print(u_analytic, len(u_analytic))
#plt.legend()
#plt.show()
with open('students1.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(u[iter])
with open('students.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(u_analytic)
    writer.writerow(u[iter])
'''with open('students2.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(u)'''

