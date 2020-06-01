#!/usr/bin/python3

### Create gaussian mixture model with three corrected components and emit proportion of the minor one 
### stdin must be piped into program
### stdin must be tab-separated and contain 3 fields:
### 1st (window name, ignored), 2nd (window length),3d number heterozygous positions in window 

from __future__ import division
import random
import os
from math import (pi,log,exp,sqrt)
import sys

### Requires extra packages: pomegranate & numpy
from pomegranate import *
import numpy as np

random.seed(11)

##### inflexion point of two normals N~(m,s) with mixture proportion prop, not used in final measure
def inflex(m1, s1, m2, s2, prop1, prop2):
        B =(m1/s1**2 - m2/s2**2)
        A = 0.5*(1/s2**2 - 1/s1**2)
        C = 0.5*(m2**2/s2**2 - m1**2/s1**2) - log((s1/s2)*(prop2/prop1))
        return((-B + np.array([1,-1])*sqrt(B**2 - 4*A*C))/(2*A))


#partial derivative of inflexion point with respect to the first mean (interchangeable)
## encoding: c=m1 d=s1 x=m2 b=s2 p=prop1 q=prop2
def deri(c,d,x,b,p,q):
        #return((b**2*c*(c-x)+(2*d**4-2*b**2*d**2)*log((d*p)/(b*q))+d**2*(c**2-3*c*x+2*x**2))/(b**2*d**2*(d**2-b**2)*sqrt(((2*d**2-2*b**2)*log((d*p)/(b*q))+c**2-2*c*x+x**2)/(b**2*d**2))))
        return(((x-c)/(d**2*sqrt(((2*d**2-2*b**2)*log((d*p)/(b*q))+c**2-2*c*x+x**2)/(b**2*d**2)))+1)/(b**2*(1/b**2-1/d**2)))

#partial derivative of inflexion point with respect to the first mixture proportion (interchangeable)
def derip(c,d,x,b,p,q):
        return(1/((p-1)*p*sqrt(((2*d**2-2*b**2)*log((d*p)/(b-b*p))+c**2-2*c*x+x**2)/(b**2*d**2))))

#partial derivative of inflexion point with respect to first standard deviation (interchangeable)
def deris(c,d,x,b,p,q):
        return((-2*c*d**2*b**2*sqrt((c**2-2*c*x+(2*d**2-2*b**2)*log((d*p)/(b-p*b))+x**2)/(d**2*b**2))+2*d**2*x*b**2*sqrt((c**2-2*c*x+(2*d**2-2*b**2)*log((d*p)/(b-p*b))+x**2)/(d**2*b**2))-c**2*d**2-c**2*b**2+2*c*d**2*x+2*c*x*b**2+d**4-d**2*x**2-2*d**2*b**2+(2*d**2*b**2-2*d**4)*log((d*p)/(b-p*b))-x**2*b**2+b**4)/(b*(d**2-b**2)**2*sqrt((c**2-2*c*x+(2*d**2-2*b**2)*log((d*p)/(b-p*b))+x**2)/(d**2*b**2))))

a=[]
c=[]
d=[]

### Read file from standard input (input must be piped):
for i in sys.stdin:    
        try:
            ii=list(map(int,i.strip().split(" ")[1:]))
        except:
            next
        if ii[0]>20000:
            a.append(ii[1]/ii[0])
            c.append(ii[0])
            d.append(ii[1])


### All steps done using Pomegrenate syntax

### Fiting of first component (normal)
#d1=PoissonDistribution(np.mean(a)*50000, frozen=True)
d1=NormalDistribution(np.mean(a),np.std(a), frozen=False)

### Fitting of second component (normal)
#d2=PoissonDistribution(4, frozen=True)
#d2=NormalDistribution(0.00015,0.00015, frozen=True)
d2=NormalDistribution(0.00015,0.00015, frozen=False)

### Fitting of third compnent (outlier collector)
#d3=NormalDistribution(np.max(d),100, frozen=True)
d3=NormalDistribution(np.max(a),0.01, frozen=False)

### Bind distributions into mixture model
model = GeneralMixtureModel([d1, d2, d3])

### EM-fit of mixture model  
model.fit(np.array(a)[np.array(a)<np.percentile(a,95)],weights=np.array(c)/sum(c))


#### correction factor: approximate integral area up to 4 standard deviations
m=d1.parameters[0]
s=d1.parameters[1]
linsp=np.linspace(-m/s-4*s,-m/s,100000)
correct=sum(NormalDistribution(0,1).probability(linsp))*np.diff(linsp)[0]

####Correct by negative distribution offset:
A=model.predict_proba([[i]for i in a])
A[:,0]=A[:,0]+correct
A=A/np.sum(A,1)[0]


### Generate proportions based on posterior classification

#resu=sum(e  for i,e in zip(model.predict([[i] for i in a]),c) if i==1)/sum(c)
resu=sum(e  for i,e in zip([np.where(u==max(u))[0][0] for u in A],c) if i==1)/sum(c)

### Attempt at correction by inflexion point parameters (not used):
try:
        ### Find inflexion point of fitted normals, there are two solutions;
        i1,i2=inflex(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0],model.distributions[1].parameters[1],resu,1-resu)
        
        #Calculate partial derivatives at inflexion point (keep negative value)
        #de=deri(model.distributions[0].parameters[0],model.distributions[0].parameters[1],i1,model.distributions[1].parameters[1], resu,1-resu)/0.1-0.9/0.1
        de=-deri(model.distributions[0].parameters[0],model.distributions[0].parameters[1],i1,model.distributions[1].parameters[1], resu,1-resu)#/0.1-0.9/0.1
        #dep=derip(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0],model.distributions[1].parameters[1], resu,1-resu)/0.00112-0.000717/0.00112
        dep=-derip(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0],model.distributions[1].parameters[1], resu,1-resu)#/0.00112-0.000717/0.00112
        #des=deris(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0],model.distributions[1].parameters[1], resu,1-resu)/1.1879-2.69/1.1879
        des=-deris(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0],model.distributions[1].parameters[1], resu,1-resu)#/1.1879-2.69/1.1879
        j=model.distributions[0].parameters[0]/2+model.distributions[1].parameters[0]/2
        
        ####Use corrected values directly (vestigial):
        #miau=(i1,i2)[abs(i1-j)>abs(i2-j)]+(i1,i2)[abs(i1-j)>abs(i2-j)]*(0.45*de+0.20*dep+0.35*des)#*model.distributions[1].parameters[0]*2
        
        ### Recaclulate inflexion point with corrected values
        i1,i2=inflex(model.distributions[0].parameters[0],model.distributions[0].parameters[1],model.distributions[1].parameters[0]+de*model.distributions[1].parameters[0],model.distributions[1].parameters[1]+des*model.distributions[1].parameters[1],resu-(1-resu)*dep,1-resu+(1-resu)*dep)
        
        ### Keep only the inflexion point between the two distributions:
        cr=(i1,i2)[abs(i1-j)>abs(i2-j)]
        
        ### print result and "scaled corrected" value (not used):
        print(resu,"\t",1/cr/10000)
except:
        ### Any cases where the operations fail will return 0\tNA
        print(0,"\tNA")

