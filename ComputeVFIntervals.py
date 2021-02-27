# Basic Python 2.7.16 code for computing the latex/tikz figures of the paper 
#
# 'Dispatching to Parallel Servers: Solutions of Poisson's Equation for First-Policy Improvement'
#
# Author: Olivier Bilenne (https://orcid.org/0000-0003-1376-8281)
#
# Generation of intervals for the cost/value functions (approximations on finite intervals)
#
import approximate_value_function as avf
import matplotlib.pyplot as plt
approximate_value_function=reload(avf)

savepath='Figures/'



test=avf.korovkin()
test.build(15)


# Bernstein approximation of sine function (test)

ar=1.0
omega=2.0      
      
tau=1
n=14
ab=avf.sin_approximationinterval()
ab.format(n,tau)
ab.compute_bernstein()
eta=ab.bernstein_approximation.errorbound()
points=100
xvec=avf.frange(0,2.2*tau,points)
yvec=[[None for i in xvec] for j in range(3)]
vvec=[None for i in xvec]
solvec=[None for i in xvec]
for i in range(len(xvec)):
    yvec[1][i]=ab.bernstein_approximation.approximation(xvec[i])
    yvec[0][i]=yvec[1][i]-eta
    yvec[2][i]=yvec[1][i]+eta
    solvec[i]=ab.bernstein_approximation.f(xvec[i])    
    dfmrate=ab.bernstein_approximation.vf_approximation_interval_expst_classicqueue(xvec[i],ar,omega)
    vvec[i]=dfmrate
#plt.figure(1)
#plt.plot(xvec,solvec,color='black')  
#plt.plot(xvec,yvec[1],color='blue')  
#plt.plot(xvec,yvec[0],'--',color='blue')  
#plt.plot(xvec,yvec[2],'--',color='blue')  
#plt.figure(2)
#plt.plot(xvec,vvec,color='red')  

print 'varrho:',test.varrho
print 'nu:',test.nu

#############################

# Optimal approximation of x^2/(a^2+x^2)

################################



a=1.0
tau=1.0
n=10
kor=avf.xsquaredoverasquaredplusxsquared_korovkin()
kor.build(n)
kor.set_a(a)
kor.compute_all(tau)

print 'gamma:',kor.gamma




a=1
tau=3.0
n=20#10
ak=avf.xsquaredoverasquaredplusxsquared_korovkin()
ak.build(n)
ak.set_a(a)
ak.compute_all(tau)
points=100
xvec=avf.frange(0,tau,points)
yvec=[[None for i in xvec] for j in range(3)]
vvec=[None for i in xvec]
solvec=[None for i in xvec]
for i in range(len(xvec)):
    yvec[1][i]=ak.fourier_approximation(xvec[i],n)
    yvec[0][i]=ak.approximation(xvec[i],n)#yvec[1][i]-eta
    yvec[2][i]=ak.approximation_bis(xvec[i],n)#yvec[1][i]+eta
    solvec[i]=ak.f(xvec[i])    
    #vvec[i]=ab.bernstein_approximation.vf_approximation_interval_expst_classicqueue(xvec[i],ar,omega)
#plt.figure(3)
#plt.plot(xvec,solvec,color='black')  
##plt.plot(xvec,yvec[1],color='blue')  
#plt.plot(xvec,yvec[0],'--',color='green')  
#plt.plot(xvec,yvec[2],'--',color='magenta')  


points=1000

ar=1.0
omega=2.0      
 
n=15
max_n=50#15
test=avf.xsquaredoverasquaredplusxsquared_approximationinterval(min(n,max_n),1.0)
test.set_a(1)
# test.build_korovkin(n) (done in init)
r=2.2
tau=1.0
test.compute_korovkin(tau,min(n,max_n))
k=n
test.format(min(k,max_n))
test.make_korovkin_approximation(min(k,max_n),k)
xvec=avf.frange(0,r,points)
yvec=[[None for i in xvec] for j in range(3)]
solvec=[None for i in xvec]
for i in range(len(xvec)):
    interval=test.f_interval(xvec[i])
    yvec[0][i]=interval[0]
    yvec[1][i]=(interval[0]+interval[1])*0.5
    yvec[2][i]=interval[1]
    solvec[i]=test.korovkin_approximation.f(xvec[i])
plt.figure(4)
plt.plot(xvec,solvec,color='black')  
#plt.plot(xvec,yvec[1],color='blue')  
plt.plot(xvec,yvec[0],'--',color='green')  
plt.plot(xvec,yvec[2],'--',color='magenta')  
#
yvec=[[None for i in xvec] for j in range(3)]
#solvec=[None for i in xvec]
for i in range(len(xvec)):
    interval=test.df_interval(xvec[i],ar,omega)
    yvec[0][i]=interval[0]
    yvec[1][i]=(interval[0]+interval[1])*0.5
    yvec[2][i]=interval[1]
    #solvec[i]=test.df_interval(xvec[i])
plt.figure(5)
#plt.plot(xvec,solvec,color='black')  
#plt.plot(xvec,yvec[1],color='blue')  
plt.plot(xvec,yvec[0],'--',color='green')  
plt.plot(xvec,yvec[2],'--',color='magenta')  

plt.show()   

#   generation of tables for LaTeX/tikz figures 

test.korovkin_latex(1.0,2.0,1,20,15,300,[0,3],[-1,1],0.2,0.95,savepath,filename='korovkin-graphs',latexaccesspath='figs/')


