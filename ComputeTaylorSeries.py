# Basic Python 2.7.16 code for computing the latex/tikz figures of the paper 
#
# 'Dispatching to Parallel Servers: Solutions of Poisson's Equation for First-Policy Improvement'
#
# Author: Olivier Bilenne (https://orcid.org/0000-0003-1376-8281)
#
# Generation of intervals for the cost/value functions (Taylor series)
#
import matplotlib.pyplot as plt
import math
from numpy import inf,Inf
#import num2word_base
import series as m
reload(m)


savepath='Figures/'
latexaccesspath=''  
latexcode=True



plt.close("all")

ar=0.5

#npower=5

# exponential st
strate=1.0
expa=0.5*float(strate-ar)
expa_=0.8*float(strate-ar)
expa__=1.0*float(strate-ar)
expasquared=0.2
#expts=m.taylorseries_exponential_service_times(ar,strate,m.oneminusexpmax(a))
expts=m.taylorseries_oneminusexpmax_exponential_service_times(ar,expa,strate)
expts_=m.taylorseries_oneminusexpmax_exponential_service_times(ar,expa_,strate)
expts__=m.taylorseries_oneminusexpmax_exponential_service_times(ar,expa__,strate)
exptssquared=m.taylorseries_oneminusexpmaxsquared_exponential_service_times(ar,expasquared,strate)

# constant st
st=1.0
csta=0.5*1.25643 #  minus the real pole of the LST of the waiting times
cstasquared=0.1
#cstts=m.taylorseries_constant_service_times(ar,st,m.oneminusexpmax(csta))
cstts=m.taylorseries_oneminusexpmax_constant_service_times(ar,csta,st)
csttssquared=m.taylorseries_oneminusexpmaxsquared_constant_service_times(ar,cstasquared,st)

nmax=25
nmin=1
points=1000
xr=[0,50]
yr=[-2,50]
colmin=0.2
colmax=0.95


zoom=2#1/float(5)

expts.display_taylor_series(nmin,nmax,1,latexcode,points,xr,yr,colmin,colmax,savepath,latexaccesspath)
expts_.display_taylor_series(nmin,nmax,2,latexcode,points,xr,yr,colmin,colmax,savepath,latexaccesspath)
expts__.display_taylor_series(nmin,nmax,3,latexcode,points,xr,yr,colmin,colmax,savepath,latexaccesspath)
#exptssquared.display_taylor_series(nmin,25,6,latexcode,points,xr,yr,colmin,colmax,savepath,latexaccesspath)
cstts.display_taylor_series(nmin,nmax,4,latexcode,points,xr,yr,colmin,colmax,savepath,latexaccesspath)
#csttssquared.display_taylor_series(nmin,25,5,latexcode,points,[xr[0]*zoom,xr[1]*zoom],[yr[0]*zoom,yr[1]*zoom],colmin,colmax,savepath,latexaccesspath)



#test=m.taylorseries_constant_service_times(ar,st,m.oneminusexpmax(a))
#test.compute_scenarios(15)
#test.compute_coefs(15)
#test.mz(15)

#ts=m.taylorseries(ar,m.exponential(strate),m.lnaplusx(a))
#ts=m.taylorseries(ar,m.exponential(2.0),m.oneminusexpmax(1.0))
#ts=m.taylorseries_exponential_service_times(ar,1.0,m.lnaplusx(80.1))
#ts=m.taylorseries_exponential_service_times(ar,strate,m.oneminusexpmax(a))
#nmax,ts=npower,m.taylorseries_exponential_service_times(ar,1.0,m.power(npower))
#ts=m.taylorseries_constant_service_times(ar,st,m.oneminusexpmax(a))

#v=[]
#for i in range(len(cstts.moments_z)/2-1):
#    v.append(math.sqrt(i)*(cstts.moments_z[i+1][1])/float(cstts.moments_z[i][1]))
    
#plt.figure(10)
#plt.plot(v)
#plt.show()
