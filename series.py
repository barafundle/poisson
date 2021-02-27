# Basic Python 2.7.16 code for computing the latex/tikz figures of the paper 
#
# 'Dispatching to Parallel Servers: Solutions of Poisson's Equation for First-Policy Improvement'
#
# Author: Olivier Bilenne (https://orcid.org/0000-0003-1376-8281)
#
# Library I
#
import math
from numpy import inf,Inf
import matplotlib.pyplot as plt
import pickle
#import num2word_base




class partition():
    def __init__(self,tau_vec,f=None):
        if len(tau_vec):
            tau_vec=sorted(tau_vec)
            self.n=1
            self.parts=[part(0,Inf,f)]
            for tau in tau_vec:
                self.add_part(tau)
    def add_part(self,tau,f=None):
        pos=self.locate(tau)
        self.parts[pos-1].interval[1]=tau
        if pos<self.n:
            self.parts.insert(pos,part(tau,self.parts[pos+1].interval[0],f))
        else:
            self.parts.append(part(tau,Inf,f))
        self.n=self.n+1
        
    def locate(self,tau):
        r=self.n
        for p in range(self.n):
            if self.parts[p].interval[0]>tau:
                r=p
                break
        return r
        
class part():
    def __init__(self,tau,tauu,f=None):
        self.f=f
        self.interval=[tau,tauu]

class mixture():
    np=0
    def __init__(self,attributes):
        self.n=0
        self.a_sc=[]
        if not attributes is None:
            for f in attributes:
                self.add_function(f[0],f[1],f[2])
    def add_function(self,a,n,coef):
        Pfound=False
        for p in self.a_sc:
            if p.a==a:
                Pfound=True
                found=False
                for c in p.components:
                    if c.n==n:
                        c.coef=c.coef+coef
                        found=True
                        break
                if not found:
                    p.components.append(component(n,coef))
                break
        if not Pfound:
            self.a_sc.append(a_specific_components(a,n,coef))
            self.np=self.np+1
    def prnt(self):
        s=''
        dejavu=False
        for p in self.a_sc:
            if dejavu:
                s=s+' + '
            if len(p.components)>1:
                s=s+'['
            dejavu=True
            dv=False
            for c in p.components:
                if dv:
                    s=s+' + '
                dv=True
                s=s+'('+str(c.coef)+')'
                if c.n:
                    if c.n==1:
                        s=s+'*u'
                    else:
                        s=s+'*u^'+str(c.n)
                else:
                    pass
            if len(p.components)>1:
                s=s+']'
            if p.a:
                if p.a==1:
                    s=s+' e^(-u)'
                elif p.a==-1:
                    s=s+' e^u'
                else:                
                    s=s+' e^('+str(-p.a)+'u)'
            else:
                pass
        return (s)
        
                    
class component():
    n=None
    coef=0
    def __init__(self,n,coef):
        self.n=n
        self.coef=coef
class a_specific_components():
    a=None
    def __init__(self,a,n=None,coef=None):
        self.a=a
        if (not n is None) and (not coef is None):
            self.components=[component(n,coef)]
        else:
            self.components=[]
        

##############"""


class taylorseries():
    def __init__(self,ar,distri,cost):
        self.distri=distri
        self.initialize_with_distri(ar,cost)
        self.closedformvf=False
    def initialize_with_distri(self,ar,cost):
        self.moments_z=None
        self.cf_poles=None
        self.ar=float(ar)
        if self.ar*self.distri.moment(1)>1:
            raise IOError('Load > 1.')
        elif self.ar*self.distri.moment(1)==1:
            print('Warning: load == 1.')
        self.cost=cost
        self.poles=self.cost.poles
        self.xak=[[] for i in self.poles]
        self.yak=[[] for i in self.poles]
        self.upperP=None
        self.lowerP=None
        self.P=None
        self.kmax=None
        self.derivtmax=None
        self.coefficient=self.ar/(1-self.ar*self.distri.moment(1))
        #self.Dcofandbounds=None
        self.Dcandbounds=None
        print 'load:', self.ar*self.distri.moment(1)
                          
    def pole_index(self,pole):
        found = False
        for i in range(len(self.poles)):
            if pole==self.poles[i]:
                found= True
                return i
        if not found:
            raise IOError('Pole not found.')
                   
    def compute_xak(self,pole,kmax):
        p=self.pole_index(pole)
        last_i=len(self.xak[p])-1
        if kmax>last_i:
            self.xak[p]=self.distri.xak(-pole,kmax)

    def compute_yak(self,pole,kmax):
        p=self.pole_index(pole)
        last_k=len(self.yak[p])-1
        if kmax>last_k:
            self.compute_xak(pole,kmax+1)
            coef=self.ar/float(1-self.distri.workload(self.ar))
            a=-pole
            if pole==0:
                for k in range(last_k+1,kmax+1):
                    if k==0:
                        self.yak[p]=[1]
                    else:
                        sum=0
                        for t in range(k):
                            sum=sum+ self.xak[p][k-t+1] * self.yak[p][t]
                        self.yak[p].append(coef*sum)                      
            else:
                for k in range(last_k+1,kmax+1):
                    if k==0:
                        self.yak[p]=[self.distri.LStransformWt(self.ar,a)]
                    elif k==1:
                        self.yak[p].append(coef*(self.distri.LStransformWt(self.ar,a)/float(a))**2*(1-self.distri.LStransform(a)-a*self.xak[p][1]))
                    else:
                        sum=(1-self.ar*self.xak[p][1])/float(self.ar)* self.yak[p][k-1]
                        #print sum
                        for t in range(k-1):
                            sum=sum- self.xak[p][k-t] * self.yak[p][t]
                            #print sum
                        self.yak[p].append(coef*sum*self.distri.LStransformWt(self.ar,a)/float(a))
                        #print    self.yak[p]                    
            
    def mz(self,kmax):
        if self.kmax is None:
            self.kmax=-1
            self.moments_z=[[0,0] for i in range(0)]
        return self.part_mz(self.kmax,kmax)

    def part_mz(self,kmin,kmax):
        v=self.moments_z+[[0,0] for i in range(kmax-kmin)]
        #fact=1
        #coef=self.ar/float(1-self.ar*self.distri.moment(1))
        for k in range(kmin+1,kmax+1):
            #fact=fact*(k+1)
            #v[k][0]=moment(k+1)/fact
            v[k][0]=self.distri.momentonfactorial(k+1)
            if k==0:
                v[0][1]=1
            else:
                s=0
                for t in range(k):
                    s=s+v[k-t][0]*v[t][1]
                v[k][1]=self.coefficient*s
        self.moments_z=v
        self.kmax=kmax
        return v
        
    def derivatives(self,t):
        if self.derivtmax is None:
            self.derivtmax=-1
            self.Dcandbounds={'Dc':[],'upperb':[],'lowerb':[]}
        if t>self.derivtmax:
            self.derivatives_part(self.derivtmax,t)
            self.derivtmax=t

#    def derivativesoverfactorials(self,n):
#            self.Dcof[k]=self.cost.Dfof(k,0)
    def derivatives_part(self,tmin,tmax):
        der=self.cost.derivatives(tmin+1,tmax)
        self.Dcandbounds['Dc']=self.Dcandbounds['Dc']+der['Dc']
        self.Dcandbounds['upperb']=self.Dcandbounds['upperb']+der['upperb']
        self.Dcandbounds['lowerb']=self.Dcandbounds['lowerb']+der['lowerb']
        
#    def compute_bounds(self,n):
#        self.mz(n+1)
#        self.derivativesoverfactorials(n)
#        self.upperP=polynom(n+2)
#        self.lowerP=polynom(n+2)
#        upperb=self.Dcofandbounds['upperb']#self.cost.upperDfof(n+1)
#        lowerb=self.Dcofandbounds['lowerb']#self.cost.lowerDfof(n+1)
#        for k in range(n):
#            term=self.Dcofandbounds['Dcof'][k]
#            for t in range(k+1,n+1):
#                term=term-self.Dcofandbounds['Dcof'][t]*self.moments_z[t-k][1]
#            self.upperP.c[k+1]=self.coefficient*(term-upperb*self.moments_z[n+1-k][1])
#            self.lowerP.c[k+1]=self.coefficient*(term-lowerb*self.moments_z[n+1-k][1])
#        self.upperP.c[n+1]=self.coefficient*(self.Dcofandbounds['Dcof'][n]-upperb*self.moments_z[1][1])
#        self.lowerP.c[n+1]=self.coefficient*(self.Dcofandbounds['Dcof'][n]-lowerb*self.moments_z[1][1])
#        self.upperP.c[n+2]=self.coefficient*(upperb)
#        self.lowerP.c[n+2]=self.coefficient*(lowerb)
#        #print self.cost.upperDfof(n+1), self.cost.lowerDfof(n+1)
#        # boundary condition
#        self.upperP.c[1]=0
#        self.lowerP.c[1]=0        
#        return 1

    def compute_bounds(self,n):
        self.mz(n+1)
        self.derivatives(n)
        self.upperP=polynom(n+2)
        self.lowerP=polynom(n+2)
        #print n, len(self.Dcandbounds['upperb']), len(self.Dcandbounds['lowerb'])
        upperb=self.Dcandbounds['upperb'][n]#self.cost.upperDfof(n+1)
        lowerb=self.Dcandbounds['lowerb'][n]#self.cost.lowerDfof(n+1)
        factinv=1.0
        for k in range(n+1):
            factinv=factinv/float(k+1)
            term=self.Dcandbounds['Dc'][k]
            for t in range(k+1,n+1):
                term=term+self.Dcandbounds['Dc'][t]*self.moments_z[t-k][1]
            self.upperP.c[k+1]=self.coefficient*(term+upperb*self.moments_z[n+1-k][1])*factinv
            self.lowerP.c[k+1]=self.coefficient*(term+lowerb*self.moments_z[n+1-k][1])*factinv
        factinv=factinv/float(n+2)
        self.upperP.c[n+2]=self.coefficient*(upperb)*factinv
        self.lowerP.c[n+2]=self.coefficient*(lowerb)*factinv
        #print self.cost.upperDfof(n+1), self.cost.lowerDfof(n+1)
        #
        # boundary condition
        #self.upperP.c[1]=0
        #self.lowerP.c[1]=0        
        #
        return 1

    def showvf(self):
        return self.cost.showvf(self.ar,self.distri)   

    #def valuefunction(self,x):
     #   return 0

    def display_taylor_series(self,nmin=1,nmax=25,fig=1,latexcode=True,points=1000,xr=[0,50],yr=[-2,40],colmin=0.2,colmax=0.95,savepath=None,latexaccesspath='figs/'):
        plt.figure(fig)
        step=(xr[1]-xr[0])/float(points)
        if latexcode:
            filename='ts-'+self.cost.id+'-'+self.distri.filename()+'-'+self.cost.filename()+'-'+str(nmax)
            if not savepath is None:
                f=file(savepath+filename+'.tex','w')
                fcost=file(savepath+filename+'-costfunction.tex','w')
                fUE=file(savepath+filename+'-UE.tex','w')
                fLE=file(savepath+filename+'-LE.tex','w')
            else:
                f=file(filename+'.tex','w')                    
            text=''
        nstep=1.0/float(nmax-nmin+1)
        upperenv=[Inf for i in range(points+1)]
        lowerenv=[-Inf for i in range(points+1)]
        for n in range(nmin,nmax+1): 
            self.compute_bounds(n)
            #
            if latexcode:
                text=text+'\\addplot[id=test,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-n)*nstep*1000)/float(10)))+'!darkcol}] { '
                text=text+self.upperP.show('x','*',False,True)+' };\n%\n'
                text=text+'\\addplot[id=test,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-n)*nstep*1000)/float(10)))+'!darkcol}] { '
                text=text+self.lowerP.show('x','*',False,True)+' };\n%\n'
            #
            k=1
            #Delta=self.upperP.c[k]-self.lowerP.c[k]
            #Deltabis=ar/float(1-ar*st)*(a**(n+1))*self.moments_z[n-k+2][1]/math.factorial(k-1)
            #print Delta, Deltabis
            col=str(colmin+(colmax-colmin)*(nmax-n)*nstep)
            x=[None for i in range(points+1)]
            upperyvec=[None for i in range(points+1)]
            loweryvec=[None for i in range(points+1)]
            for k in range(points+1):
                x[k]=xr[0]+k*step
                upperyvec[k]=self.upperP.P(x[k])    
                loweryvec[k]=self.lowerP.P(x[k])
                upperenv[k]=min(upperenv[k],upperyvec[k])    
                lowerenv[k]=max(lowerenv[k],loweryvec[k]) 
            plt.plot(x,upperyvec,'-',color=col,linewidth=0.5)
            plt.plot(x,loweryvec,'-',color=col,linewidth=0.5)
        plt.plot(x,upperenv,'-',color='black',linewidth=2.0)
        plt.plot(x,lowerenv,'-',color='black',linewidth=2.0)
        #
        if latexcode:
            if self.closedformvf:
                text=text+'\\addplot[id=lcurve,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, {\\solstyle}, smooth,opacity=1.0,line width={\\solwidth}, mark = none, mark options={scale=0.1}, color = black] { '
                text=text+self.showvf()+' };\n%\n'
            #text=text+'\\addplot[id=ucurve,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, smooth,opacity=1.0,line width={\\envelopewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-n)*nstep*1000)/float(10)))+'!darkcol}] { '
            #text=text+self.upperP.show('x','*',False)+' };\n%\n'
            curveUE="{:>15}  {:>15}\n".format('x','y')
            curveLE="{:>15}  {:>15}\n".format('x','y')
            for p in range(len(x)):
                curveUE+="{:>15}  {:>15}\n".format(x[p],upperenv[p])      
                curveLE+="{:>15}  {:>15}\n".format(x[p],lowerenv[p])      
            fUE.write(curveUE)
            fLE.write(curveLE)
            fcost.write('\\addplot[id=costfunction,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, smooth,opacity=1.0,line width={\\envelopewidth}, mark = none, mark options={scale=0.1}, color = black] { '+self.cost.latex_f()+' };\n%\n')
            fUE.close()
            fLE.close()
            fcost.close()
            text=text+'\\addplot[id=uenvelope,domain={\\xlim}:{\\xlimm}, smooth,opacity=1.0,line width={\\envelopewidth}, mark = none, mark options={scale=0.1}, color = black] table{'
            text=text+latexaccesspath+filename+'-UE'+'.tex};\n%\n'
            text=text+'\\addplot[id=lenvelope,domain={\\xlim}:{\\xlimm}, smooth,opacity=1.0,line width={\\envelopewidth}, mark = none, mark options={scale=0.1}, color = black] table{'
            text=text+latexaccesspath+filename+'-LE'+'.tex};\n%\n'
            f.write(text)
            f.close()
        #
        if self.closedformvf:
            sol=[0 for i in range(points+1)]
            for k in range(points+1):
                sol[k]=self.valuefunction(x[k])
            plt.plot(x,sol,':',color='black',linewidth=1.0)
        plt.ylim(yr)
        plt.show()
                
class function():
    def f(self,x):
        return None
    def Df(self,n,x):
        return None
    def Dfzero(self,n):
        return self.Df(n,0)
    def upperDf(self,n):
        return None
    def lowerDf(self,n):
        return None
    def derivatives(self,nmin,nmax): # can be modified if factorials are an issue
        Dc=[None for i in range(nmin,nmax+1)]
        upperb=[None for i in range(nmin,nmax+1)]
        lowerb=[None for i in range(nmin,nmax+1)]
        for k in range(0,nmax-nmin+1):
            Dc[k]=self.Df(nmin+k,0)
            upperb[k]=self.upperDf(nmin+k+1)
            lowerb[k]=self.lowerDf(nmin+k+1)
        return {'Dc':Dc,'upperb':upperb,'lowerb':lowerb}
    def Dfof(self,n,x):
        return None
    def Dfzeroof(self,n):
        return self.Dfof(n,0)
    def upperDfof(self,n):
        return None
    def lowerDfof(self,n):
        return None
    def derivativesoverfactorials(self,n): # can be modified if factorials are an issue
        Dcof=[None for i in range(n+1)]
        for k in range(n+1):
            #Dcof[k]=self.Dfof(k,0)
            Dcof[k]=self.Dfzeroof(k)
        return {'Dcof':Dcof,'upperb':self.upperDfof(n+1),'lowerb':self.lowerDfof(n+1)}

class polynom():
    def __init__(self,n):
        self.n=n
        self.c=[None for i in range(n+1)]
    def coef(self,k,x):
        if k>(self.n):
            self.c=self.c+[None for i in range(k-self.n)]
            self.n=k
        self.c[k]=x
    def P(self,x):
        y=0
        for k in range(self.n+1):
            if not (self.c[k] is None):
                y=y+self.c[k]*x**k
        return y
    def show(self,x,times_op='',header=True, scale=False): 
        if not times_op=='':
            times_op=' '+times_op
        x=str(x)
        if header:
            st='P('+x+') ='
        else:
            st=''
        dejavu=False
        for i in range(len(self.c)):
            if not(self.c[i] is None):
                st=st+' '
                if dejavu:
                    st=st+'+ '
                dejavu=True
                if scale and i>1:
                    if self.c[i]<0:
                        st=st+'(-1)*'
                    st=st+'('+str(abs(self.c[i])**(1/float(i)))
                    st=st+times_op+' '+x
                    st=st+')^'+str(i)
                else:
                    st=st+'('+str(self.c[i])+')'
                    if i:
                        st=st+times_op+' '+x
                        if i>1:
                            st=st+'^'+str(i)
        return(st)
        
    


class distribution():
    def moment(self,n):
        return 0 
    def momentonfactorial(self,n):
        return 0
    
##############################################################################
##############################################################################
################### cost functions ###########################################
##############################################################################
##############################################################################
##############################################################################


class power(function):
    id='power'        
    def __init__(self,n):
        if n>=0:
            self.n=n
            self.poles=[0]
        else:
            raise IOError('Argument of constructor power(n) should be nonnegative.')
    def f(self,x):
        return x**self.n
    def Df(self,n,x):
        if n==self.n:
            return math.factorial(self.n)
        else:
            return 0.0
    def upperDf(self,n):
        if n==self.n:
            return math.factorial(n)
        elif n>self.n:
            return 0
        else:
            return Inf 
    def lowerDf(self,n):
        if n==self.n:
            return math.factorial(n)
        else:
            return 0
    def Dfof(self,n,x):
        if n==self.n:
            return 1.0
        else:
            return 0
    def upperDfof(self,n):
        if n==self.n:
            return 1
        elif n>self.n:
            return 0
        else:
            return Inf 
    def lowerDfof(self,n):
        if n==self.n:
            return 1.0
        else:
            return 0  
      

class oneminusexpmax(function):
    id='oneminusexpmax'        
    def __init__(self,a):
        if a>0:
            self.a=a
            self.poles=[0,-a]
        else:
            raise IOError('Argument of constructor oneminusexpmax(a) should be positive.')
    def filename(self):
        return str(int(self.a*1000))
    def latex_f(self):
        return '1- exp(-'+str(self.a)+'*x)'
    def f(self,x):
        return 1-math.exp(-self.a*x)
    def Df(self,n,x):
        if n==0:
            return self.f(x)
        else:
            return -(-self.a)**n*math.exp(-self.a*x)
    def upperDf(self,n):
        if int(n)%2==0:
            return 0
        else:
            return self.Df(n,0)
    def lowerDf(self,n):
        if int(n)%2==1:
            return 0
        else:
            return self.Df(n,0)
    def Dfof(self,n,x):
        if n==0:
            return self.f(x)
        else:
            return -(-self.a)**n*math.exp(-self.a*x)/math.factorial(n)
    def upperDfof(self,n):
        if int(n)%2==0:
            return 0
        else:
            return self.Dfof(n,0)
    def lowerDfof(self,n):
        if int(n)%2==1:
            return 0
        else:
            return self.Dfof(n,0)
    def derivativesoverfactorials(self,n): # can be modified if factorials are an issue
        Dcof=[None for i in range(n+1)]
        deriv=-1.0
        for k in range(n+1):
            #Dcof[k]=self.Dfof(k,0)
            if k==0:
                Dcof[k]=self.f(0)
            else:
                deriv=deriv*float(-self.a)/float(k)
                #Dcof[k]=-(-self.a)**k/math.factorial(k)
                Dcof[k]=deriv
        deriv=deriv*float(-self.a)/float(n+1)
        upperb=None
        if int(n+1)%2==0:
            upperb=0
        else:
            #upperb=self.Dfof(n+1,0)
            upperb=deriv
        lowerb=None
        if int(n+1)%2==1:
            lowerb=0
        else:
            #lowerb=self.Dfof(n+1,0)
            lowerb=deriv
        return {'Dcof':Dcof,'upperb':upperb,'lowerb':lowerb}
    def showvf(self,ar,distri):
        return('1.0/('+str(self.a)+'/'+str(ar)+'+'+str(distri.LStransform(self.a))+'-1)*(exp(-'+str(self.a)+'* x)-1) + '+str(ar)+'/(1-'+str(ar)+'*'+str(distri.moment(1))+')*x')         

class oneminusexpmaxsquared(function):
    id='oneminusexpmaxsquared'        
    def __init__(self,a):
        if a>0:
            self.a=a
            self.poles=[0,-a]
        else:
            raise IOError('Argument of constructor oneminusexpmaxsquared(a) should be positive.')
    def filename(self):
        return str(int(self.a*1000))
    def latex_f(self):
        return '1- exp(-'+str(self.a)+'*x^2)'
    def f(self,x):
        return 1-math.exp(-self.a*(x**2))
    def Dfzero(self,n):
        if int(n)%2==0:
            if n==0:
                return 0.0
            else:
                return -float(math.factorial(n))/float(math.factorial(n/2))*(-self.a)**(n/2)
        else:
            return 0
    #def upperDf(self,n):
    #    return max(self.Dfzero(n),0)
    #def lowerDf(self,n):
    #    return min(self.Dfzero(n),0)
    def Dfzeroof(self,n):
        if int(n)%2==0:
            if n==0:
                return 0.0
            else:
                return -1.0/float(math.factorial(n/2))*(-self.a)**(n/2)
        else:
            return 0
    #def upperDfof(self,n):
    #    return max(self.Dfzeroof(n),0)
    #def lowerDfof(self,n):
    #    return min(self.Dfzeroof(n),0)
    def derivatives(self,nmin,nmax): # can be modified if factorials are an issue
        Dc=[None for i in range(nmin,nmax+1)]
        upperb=[None for i in range(nmin,nmax+1)]
        lowerb=[None for i in range(nmin,nmax+1)]
        deriv=-1.0
        for k in range(0,nmin):
            if int(k)%2==1:
                deriv=deriv*float(-self.a)*k*2 #deriv*float(-self.a)*(k+1)*k/float((k+1)/2) 
        for k in range(0,nmax-nmin+1):
            if k<5:
                coef=10.0
            else:
                coef=1.0    
            #Dcof[k]=self.Dfzero(k)
            if int(k+nmin)%2==0:
                if (k+nmin)==0:
                    Dc[k]=0
                else:
                    Dc[k]=deriv #-float(math.factorial(k+nmin))/float(math.factorial((k+nmin)/2))*(-self.a)**((k+nmin)/2)
            else:
                deriv=deriv*float(-self.a)*(k+nmin)*2 #deriv*float(-self.a)*(k+nmin+1)*(k+nmin)/float((k+nmin+1)/2)
                Dc[k]=0
            if k>0:
                upperb[k-1]=abs(deriv)*coef
                lowerb[k-1]=-abs(deriv)*coef
        # take as bounds the maximum amplitude of derivative (at 0), given by |deriv| 
        if int(nmax+1)%2==1:
            deriv=deriv*float(-self.a)*2*(nmax+1) #deriv*float(-self.a)*(nmax+2)*(nmax+1)/float((nmax+2)/2)           
        upperb[nmax-nmin]=abs(deriv)*coef
        lowerb[nmax-nmin]=-abs(deriv)*coef
        return {'Dc':Dc,'upperb':upperb,'lowerb':lowerb}
    def showvf(self,ar,distri):
        return None
        
class lnaplusx(function):
    id='lnaplusx'
    def __init__(self,a):
        if a>0:
            self.a=a
        else:
            raise IOError('Argument of constructor lnaplusx(a) should be positive.')
    def f(self,x):
        return math.log(self.a+x)
    def Df(self,n,x):
        if n==0:
            return self.f(x)
        else:
            return (-1)**(n+1)*math.factorial(n-1)/(self.a+x)**n
    def upperDf(self,n):
        if int(n)%2==0:
            return 0
        else:
            return self.Df(n,0)
    def lowerDf(self,n):
        if int(n)%2==1:
            return 0
        else:
            return self.Df(n,0)
    def derivatives(self,n): # can be modified if factorials are an issue
        Dc=[None for i in range(n+1)]
        deriv=1.0/float(self.a)
        for k in range(n+1):
            #Dc[k]=self.Df(k,0)
            if k==0:
                Dc[k]=self.f(0)
            else: 
                 #Dc[k]=(-1)**(k+1)*math.factorial(k-1)/(self.a)**k = -Dc[k-1]*float(k-1)/float(self.a)
                Dc[k]=deriv
                deriv=-deriv*float(k)/float(self.a)
        upperb=None
        if int(n+1)%2==0:
            upperb=0
        else:
            #upperb=self.Df(n+1,0)
            upperb=deriv
        lowerb=None
        if int(n+1)%2==1:
            lowerb=0
        else:
            #lowerb=self.Df(n+1,0)
            lowerb=deriv
        return {'Dc':Dc,'upperb':upperb,'lowerb':lowerb}
    def Dfof(self,n,x):
        if n==0:
            return self.f(x)
        else:
            return (-1)**(n+1)/float(n)/(self.a+x)**n
    def upperDfof(self,n):
        if int(n)%2==0:
            return 0
        else:
            return self.Dfof(n,0)
    def lowerDfof(self,n):
        if int(n)%2==1:
            return 0
        else:
            return self.Dfof(n,0)
            
            
            
#####################################################
#####################################################
################# service time distributions ########
#####################################################
#####################################################


class constant(distribution):  
    id='cst'
    # number of poles
    np=Inf
    def __init__(self,d):
        self.d=float(d) 
    def filename(self):
        return self.id+str(int(self.d*10))
    def moment(self,n):
        return self.d**n
    def workload(self,ar):
        return(self.d*ar)
    def momentonfactorial(self,n):
        return self.d**n/float(math.factorial(n))
    def robust_momentonfactorial(self,n):
        mof=1.0
        for i in range(1,n+1):
            mof=mof*float(self.d)/float(i)
        return mof
    def LStransform(self,s):
        return math.exp(-s*self.d)
    def LStransformWt(self,ar,s):
        if (self.workload(ar)<1):
            if (s==0):
                return 1
            else:
                if (s-ar*(1.0-self.LStransform(s)))==0 :
                    print 'Warning: Calling LST transform of the waiting times at a pole'
                    return Inf
                else:
                    return (1-ar*self.d)*s/float(s-ar*(1.0-self.LStransform(s)))
        else:
            raise IOError('Workload >= 1')
    def xak(self,a,kmax):
        vec=[None for i in range(kmax+1)]
        x=1
        if a==0:
            emad=1
        else:
            emad=math.exp(-a*self.d)
        for i in range(kmax+1):
            vec[i]=x*emad
            x=x*self.d/float(i+1)
        return vec
                
class exponential(distribution): 
    id='exp' 
    def __init__(self,lamb):
        self.lamb=float(lamb)       
    def filename(self):
        return self.id+str(int(self.lamb*10))
    def moment(self,n):
        return math.factorial(n)/float(self.lamb**n) 
    def workload(self,ar):
        return(ar/float(self.lamb))
    def momentonfactorial(self,n):
        return 1/float(self.lamb**n)
    def LStransform(self,s):
        # \int_0^\infty lamb*\exp(-lamb*x) * \exp(-s*x) dx = lamb/(lamb+s) (lamb+as0)
        if (self.lamb+s)>0:
            return self.lamb/float(self.lamb+s)
        else:
            raise IOError('Exponential distribution parameter and a must be positive')
    def LStransformWt(self,ar,s):
        if (self.workload(ar)<1):
            if (s+self.lamb-ar==0):
                print 'Warning: Calling LST transform of the waiting times at a pole'
                return Inf
            else:
                return (self.lamb-ar)*(s+self.lamb)/float(self.lamb*(s+self.lamb-ar))
        else:
            raise IOError('Exponential distribution parameter must be larger than the arrival rate')
    #number of poles
    np=1
    # [[pole_1,degree_1],...,[pole_np,degree_np]]
    def poles(self,ar):
        return([[ar-self.lamb,1]])
    # coefficient in (s-p_np)^(-q) of the Laurent series at p_np of LStransformWt * (-1)^(degree of p_np) (q=0,...,degree of p_np-1)
    def LaurentWt(self,ar,pn,q):
        if pn<self.np:
            if q<self.poles(ar)[pn][1]:
                return (ar/float(self.lamb)-1.0)*ar
            else:
                raise IOError('q out of range.')
        else:
            raise IOError('pn out of range.')
    # W(-s)=E[exp(s*X)]=\int w exp(-(w-s)*x) dx = w/(w-s)
    # d^k/ds^k W(-s) =E[x^k exp(s*X)] = k! w/(w-s)^(k+1)
    # xa:k= E[x^k exp(-a*X)] / k! = w/(w+a)^(k+1)
    def xak(self,a,kmax):
        vec=[None for i in range(kmax+1)]
        x=self.lamb/float(self.lamb+a)
        for i in range(kmax+1):
            vec[i]=x
            x=x/float(self.lamb+a)
        return vec


##################################################################
###################################################################
##################################################################

class taylorseries_constant_service_times(taylorseries):
    def __init__(self,ar,d,cost):
        self.distri=constant(d)
        self.initialize_with_distri(ar,cost)
        self.scenarios=None
        self.scenarios_of=None
        self.coefs=None
    def compute_scenarios(self,n):
        start=0
        if self.scenarios is None:
            self.scenarios=[None for i in range(n+1)]
            self.scenarios_of=[None for i in range(n+1)]
            start=2
        else:
            start=len(self.scenarios)
            self.scenarios=self.scenarios+[None for i in range(n-start+1)]
            self.scenarios_of=self.scenarios_of+[None for i in range(n-start+1)]
        for i in range(start,n+1):
            self.scenarios[i]=[None for j in range(i/2+1)]
            self.scenarios[i][1]=1
            self.scenarios_of[i]=[None for j in range(i/2+1)]
            if i==2:
                self.scenarios_of[i][1]=0.5
            else:
                self.scenarios_of[i][1]=self.scenarios_of[i-1][1]/float(i)
            for j in range(2,i/2+1):
                m=j-1
                sum=0
                sum_of=0
                binomial=float(i)
                factnmp=float(1)
                for p in range(i-2,2*m-1,-1):
                    binomial=binomial*(p+1)/float(i-p)
                    sum=sum+binomial*self.scenarios[p][m]
                    factnmp=factnmp/float(i-p)
                    sum_of=sum_of+factnmp*self.scenarios_of[p][m]                    
                self.scenarios[i][m+1]=sum
                self.scenarios_of[i][m+1]=sum_of
        return 0
    def compute_coefs(self,n):
        coef=self.ar/float(1-self.ar*self.distri.moment(1))
        start=0
        coefproduct=1.0
        if self.coefs is None:
            self.coefs=[None for i in range(n+1)]
        else:
            start=len(self.coefs)
            coefproduct=self.coefs[-1]*coef
            self.coefs=self.coefs+[None for i in range(n-start+1)]
        for i in range(start,n+1):
            self.coefs[i]=coefproduct
            coefproduct=coefproduct*coef
    def part_mz(self,kmin,kmax):
        self.moments_z=self.moments_z+[[None,None] for i in 2*range(kmax-kmin)]
        # compute moments 
        mof=1.0
        if kmin>=0:
            mof=self.moments_z[2*kmin+1][0]*float(self.distri.d)#/float(kmin+1)
        for t in range(2*kmin+2,2*kmax+2):
            self.moments_z[t][0]=mof
            mof=mof*float(self.distri.d)#/float(t+1)
        # moments on factorials ready
        self.compute_scenarios(2*kmax)
        self.compute_coefs(kmax)
        for k in range(kmin+1,kmax+1):
            if k==0:
                self.moments_z[0][1]=1
            else:
                sum=0
                sl=0
                su=0
                for t in range(1,k+1):
                    sum=sum+self.coefs[t]*self.scenarios_of[k+t][t]*self.moments_z[k+t][0]              
                    if t<=math.sqrt(k):
                        sl=sl+self.coefs[t]*self.scenarios_of[k+t][t]*self.moments_z[k+t][0]
                        tt=int(math.floor(t*math.sqrt(k)))
                        #print (self.coefs[tt]*self.scenarios_of[k+tt][tt]*self.moments_z[k+tt][0])/float(self.coefs[t]*self.scenarios_of[k+t][t]*self.moments_z[k+t][0])
                    else:
                        su=su+self.coefs[t]*self.scenarios_of[k+t][t]*self.moments_z[k+t][0]
                #print su/float(sl)
                self.moments_z[k][1]=sum
        self.kmax=kmax
        return self.moments_z
    # the next function is an early , equivalent version of part_mz but with numerical issues due to the factorials
    def part_mz_old(self,kmin,kmax):
        self.moments_z=self.moments_z+[[None,None] for i in 2*range(kmax-kmin)]
        # compute moments on factorials
        mof=1.0
        if kmin>=0:
            mof=self.moments_z[2*kmin+1][0]*float(self.distri.d)/float(kmin+1)
        for t in range(2*kmin+2,2*kmax+2):
            self.moments_z[t][0]=mof
            mof=mof*float(self.distri.d)/float(t+1)
        # moments on factorials ready
        self.compute_scenarios(2*kmax)
        self.compute_coefs(kmax)
        for k in range(kmin+1,kmax+1):
            if k==0:
                self.moments_z[0][1]=1
            else:
                sum=0
                for t in range(1,k+1):
                    sum=sum+self.coefs[t]*self.scenarios[k+t][t]*self.moments_z[k+t][0]              
                self.moments_z[k][1]=sum
        self.kmax=kmax
        return self.moments_z


class taylorseries_exponential_service_times(taylorseries):
    def __init__(self,ar,lamb,cost):
        self.distri=exponential(lamb)
        self.initialize_with_distri(ar,cost)
    def part_mz(self,kmin,kmax):
        v=self.moments_z+[[0,0] for i in range(kmax-kmin)]
        #fact=1
        #coef=self.ar/float(1-self.ar*self.distri.moment(1))
        prod=1/float(self.distri.lamb-self.ar)**(kmin)
        for k in range(kmin+1,kmax+1):
            #fact=fact*(k+1)
            #v[k][0]=moment(k+1)/fact
            v[k][0]=self.distri.momentonfactorial(k+1)
            if k==0:
                v[0][1]=1
                prod=1/float(self.distri.lamb-self.ar)**(kmin+1)
            else:
                ### faster computation
                #s=0
                #for t in range(k):
                #    s=s+v[k-t][0]*v[t][1]
                #v[k][1]=self.coefficient*s
                prod=prod/float(self.distri.lamb-self.ar)
                v[k][1]=self.ar/self.distri.lamb*prod
        self.moments_z=v
        self.kmax=kmax
        return v

class taylorseries_oneminusexpmax(taylorseries):
    closedformvf=True
    def __init__(self,ar,distri,a):
        self.distri=distri
        self.initialize_with_distri(ar,oneminusexpmax(a))
    def denominator(self):
        return self.cost.a/float(self.ar)+self.distri.LStransform(self.cost.a)-1.0
        # for exponential service times: denom = a/ar + strate/(strate+a) - 1 
        #                                      = [a*(strate+a)+strate*ar-ar*(strate+a)]/ar/(strate+a)
        #                                      = a*[strate+a-ar]/ar/(strate+a)
    def valuefunction(self,x):
        return (math.exp(-self.cost.a*x)-1.0)/float(self.denominator())+(float(self.ar)/(1-float(self.ar)*self.distri.moment(1)))*x

class taylorseries_oneminusexpmax_constant_service_times(taylorseries_constant_service_times,taylorseries_oneminusexpmax):
    def __init__(self,ar,a,d):
        self.distri=constant(d)
        self.initialize_with_distri(ar,oneminusexpmax(a))
        self.scenarios=None
        self.scenarios_of=None
        self.coefs=None

class taylorseries_oneminusexpmax_exponential_service_times(taylorseries_exponential_service_times,taylorseries_oneminusexpmax):
    def __init__(self,ar,a,lamb):
        self.distri=exponential(lamb)
        self.initialize_with_distri(ar,oneminusexpmax(a))

class taylorseries_oneminusexpmaxsquared(taylorseries):
    closedformvf=False
    def __init__(self,ar,distri,a):
        self.distri=distri
        self.initialize_with_distri(ar,oneminusexpmaxsquared(a))
    def valuefunction(self,x):
        return None
        
class taylorseries_oneminusexpmaxsquared_constant_service_times(taylorseries_constant_service_times,taylorseries_oneminusexpmaxsquared):
    def __init__(self,ar,a,d):
        self.distri=constant(d)
        self.initialize_with_distri(ar,oneminusexpmaxsquared(a))
        self.scenarios=None
        self.scenarios_of=None
        self.coefs=None

class taylorseries_oneminusexpmaxsquared_exponential_service_times(taylorseries_exponential_service_times,taylorseries_oneminusexpmaxsquared):
    def __init__(self,ar,a,lamb):
        self.distri=exponential(lamb)
        self.initialize_with_distri(ar,oneminusexpmaxsquared(a))

    
