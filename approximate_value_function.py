# Basic Python 2.7.16 code for computing the latex/tikz figures of the paper 
#
# 'Dispatching to Parallel Servers: Solutions of Poisson's Equation for First-Policy Improvement'
#
# Author: Olivier Bilenne (https://orcid.org/0000-0003-1376-8281)
#
# Library II
#
from math import cos,acos,sin,pi,factorial,exp,sqrt,atan
from scipy.special import binom,gamma


def frange(start, stop, n):
    v=[]
    step=(stop-start)/float(n)
    for i in range(n+1):
        v.append(i*step+start)
    return v
    

class korovkin:
    n=0
    built_for=-1
    def __init__(self,n=1):
        self.n=0
        self.varrho=[[1]] # varrho_{n,k}
        self.nu=[[1]]  # nu(k,q)
        self.gamma=None # gamma(n,t)
        self.tildegamma=None # tildegamma(n,t)
        #unnamed=[[1]]  # unnamed coefffcient (q,j)
        self.laurent=None
        self.fourier_coefs_ready=False
        self.format(n)
    def format(self,n=None,tau=None):
        self.fourier_coefs_ready=False
        if n>self.built_for:
            self.build(n)
            self.built_for=n
    def is_even(self,x):
        if (x%2==0):
            return(True)
        else:
            return(False)
    def fourier_beta(self,k):
        if not self.fourier_coefs_ready:
            print "fourier_beta coefficients not yet computed. First use compute_fourier_coefficients(self,tau)."
            return(0)
        else:
            return(self.fourier_beta_coef[k])
    def compute_fourier_coefficients(self,tau):
        print "compute_fourier_coefficients() not yet specified."
    def compute_all(self,tau):
        self.tau=tau
        self.compute_fourier_coefficients(tau)
        self.compute_gamma()
    def build(self,n=1):
        if n>0:
            for i in range(self.n+1,n+1):
                rho=[[] for j in range(i+1)]
                rho[0]=1
                rho[1]=cos(pi/float(i+2))
                denom=0
                for q in range(i+1):
                    denom=denom+sin(pi*(q+1)/float(i+2))**2
                for k in range(2,i+1):
                    num=0
                    for q in range(0,i-k+1):
                        num=num+sin(pi*(q+1)/float(i+2))*sin(pi*(q+k+1)/float(i+2))                     
                    rho[k]=num/float(denom)
                self.varrho.append(rho)   
            #     
            self.n=n
            #
            self.nu=[[1]]
            matrix=[[None for t in range((n/2+1))] for q in range((n/2+1))]
            for k in range(1,n+1):
                nuk=[None for q in range(k/2+1)]
                for q in range(k/2+1):
                    if (q==0):
                        if (self.is_even(k)):
                            matrix[0][0]=(-1)**(k/2)
                        else:
                            matrix[0][0]=(-1)**(k/2)*k
                    else:
                        c=(-1)^(k/2-q)                                
                        for t in range(q):
                            matrix[q][t]=(-1)*matrix[q-1][t]*(k/2-(q-1))/(q-t)
                        if (self.is_even(k)):
                            matrix[q][q]=(-1)*matrix[q-1][q-1]*(k-2*q+2)*(k-2*q+1)/(2*q*(2*q-1))
                        else:
                            matrix[q][q]=(-1)*matrix[q-1][q-1]*(k-2*q+1)*(k-2*q)/((2*q+1)*(2*q))
                    nuk[q]=0
                    for t in range(q+1):
                        nuk[q]+=matrix[q][t]
                self.nu.append(nuk)
            #
        self.fourier_coefs_ready=False
    def compute_gamma(self):
        self.gamma=[None for n in range(self.n+1)]      
        for n in range(self.n+1):
            self.gamma[n]=[None for t in range(n+1)]
            for t in range (n+1):
                self.gamma[n][t]=0
                for k in range(t,n+1,2):
                    #self.gamma[n][t]+=1*self.fourier_beta(k)*self.nu[k][t/2]  # !!!!!!!
                    self.gamma[n][t]+=self.varrho[n][k]*self.fourier_beta(k)*self.nu[k][t/2]
        self.tildegamma=[None for n in range(self.n+1)]      
        for n in range(self.n+1):
            self.tildegamma[n]=[None for k in range(n+1)]
            coef=1;
            for k in range (n+1):
                bino=1
                self.tildegamma[n][k]=0
                for t in range(n-k+1):
                    self.tildegamma[n][k]+=bino*self.gamma[n][t+k]
                    bino=-bino*(t+1+k)/(t+1)
                self.tildegamma[n][k]*=coef
                coef*=2/float(self.tau)
    def f(self,x):
        return 1
    def continuitymodulus(self,delta,tau):
        return 0
    def errorbound(self,n,tau=None):
        if tau is None:
            tau=self.tau
        return 6*self.continuitymodulus(tau/float(2*n),tau)
    def fourier_approximation(self,x,n=None):
        if not self.fourier_coefs_ready:
            self.compute_all()
        if n is None:
            n=self.n
        else:
            n=min(n,self.n)
        s=0
        theta=acos(2*x/float(self.tau)-1)
        for k in range(n+1):
            s+=self.fourier_beta(k)*cos(k*theta)
        return s
    def approximation(self,x,n=None):
        if not self.fourier_coefs_ready:
            self.compute_all()
        if n is None:
            n=self.n
        else:
            n=min(n,self.n)
        s=0
        powx=1.0
        for t in range(n+1):
            s+=self.gamma[n][t]*powx   #
            powx=powx*float(2*x/float(self.tau)-1)
        return s
    def approximation_bis(self,x,n=None):
        if not self.fourier_coefs_ready:
            self.compute_all()
        if n is None:
            n=self.n
        else:
            n=min(n,self.n)
        s=0
        powx=1.0
        for k in range(n+1):
            s+=self.tildegamma[n][k]*powx   #
            powx=powx*x
        return s
    def coefficient(self,n,k):
        if (n<=self.n and k<=n):
            return self.tildegamma[n][k]
        else:
            print "Wrong arguments",n,k
            
        

    
class bernstein:
    def __init__(self,n=None,tau=None):
        self.tau=tau
        self.n=n
        self.beta=[]
        self.format(n,tau)
        self.ready=False
    def format(self,n=None,tau=None):
        if not n is None:
            self.n=n
        if not tau is None:
            self.tau=tau
        if not self.n is None:
            self.beta=[None for i in range(self.n+1)]
            self.smallvalues=[None for i in range(self.n+1)]
            for i in range(self.n+1):
                self.smallvalues[i]=realexponentialcomponent(None,i,0)
        self.ready=False
    def f(self,x):
        return 1
    def continuitymodulus(self,delta,tau):
        return 0
    def errorbound(self):
        return 1.5*self.continuitymodulus(self.tau/float(sqrt(self.n)),self.tau)
    def compute_beta(self,tau=None,n=None):
        if not tau is None:
            self.tau=tau
        if not n is None:
            self.n=n
        self.format()
        if not (self.n is None or self.tau is None):
            for k in range(self.n+1):
                s=0
                denom=1/float((-self.tau)**k)
                for t in range(k+1):            
                    #s=s+binom(self.n,t)*binom(self.n-t,k-t)*(-1)**t*self.f(t*self.tau/float(self.n))
                    s=s+binom(k,t)*(-1)**t*self.f(t*self.tau/float(self.n))
                    #print binom(self.n,t)*binom(k,t)/float(binom(self.n,t)*binom(self.n-t,k-t)),self.n,k,t
                self.beta[k]=binom(self.n,k)*s*denom
                self.smallvalues[k].k=self.beta[k]
            self.ready=True
    def approximation(self,x):
        if not self.ready:
            self.compute_beta()
        s=0
        powx=1.0
        for k in range(self.n+1):
            s=s+self.beta[k]*powx
            powx=powx*float(x)
        return s
    def vf_approximation_interval_expst_classicqueue(self,x,ar,omega):
        if not self.ready:
            self.compute_beta()
        s=0
        for k in range(self.n+1):
            [df,rate]=self.smallvalues[k].interval_vf_expst_classicqueue(x,self.tau,ar,omega)
            s+=df-rate
        return s
        
    
 
    
class realexponentialcomponent:
    def __init__(self,k=None,n=None,a=None):
        self.k=k
        self.n=n
        self.a=a
    def show(self):
        return [self.k,self.n,self.a]
    def f(self,x):
        return(self.k*(x**self.n))*exp(-self.a*x)
    def latex_f(self):
        return str(self.k)+'*(x^'+str(self.n)+')*exp(-'+str(self.a)+'*x)'
    def largevalues_f(self,x,tau):
        if x>tau:
            return(self.f(x))
        else:
            return(0)
    def largevalues_vf_expst_classicqueue(self,u,tau,ar,omega):
        if self.k==0:
            # cost function is 0 for u > tau
            return [0,0]
        else:
            # cost function is exponential (k,a) for u > tau
            dif=float(omega-ar)
            difpa=dif+self.a
            x=min(u,tau)
            coef=ar*ar*exp(-difpa*tau)/float(difpa)
            s=0
            sterm=1.0
            for t in range(self.n+1):
                coef*=t/float(difpa)
                s+=sterm
                sterm*=(tau*difpa/float(t+1)) 
            # computation of the derivative at zero of the drain function
            deriv_at_0=coef*s
            # computation of the drain function
            #   from small values
            df=deriv_at_0*(exp(dif*x)-1)/float(dif) 
            #   from large values
#             table=[None for t in range(self.n+1)]
#             table[self.n]=ar/float(difpa)
#             for q in range(self.n-1,-1,-1):
#                 table[q]=table[q+1]*(q+1)/float(difpa*tau)
#             for q in range(self.n-1,-1,-1):
#                 table[q]+=table[q+1]
#             if (self.a>0):
#                 mul=1.0
#                 for j in range(self.n+1):
#                     table[j]=(table[j]+1)*mul
#                     mul*=(j+1)/float(self.a*tau)
#                 for t in range(self.n-1,-1,-1):    
#                     table[t]+=table[t+1]
#                 s=0
#                 mul=1
#                 y=self.a*max(u,tau)
#                 for t in range(self.n+1):
#                     s=table[t]*mul
#                     mul*=(y-self.a*tau)/float(t+1)
#                 df=df+ar*(tau**self.n)/float(self.a)*(-exp(-y)*s+exp(-self.a*tau)*table[0])  
#             else:   # self.a==0
#                 z=max((u-tau)/float(tau),0)
#                 mul=z
#                 for t in range(self.n+1):
#                     table[t]=(table[t]+1)*mul
#                     mul*=z
#                 for t in range(self.n+1):
#                     table[t]=table[t]/float(t+1) 
#                 df=df+ar*(tau**(self.n+1))*sum(table)
            # 
            if (self.a>0):
                maxx=max(u,tau)
                expmat=exp(-self.a*tau)
                expmmaxx=exp(-self.a*maxx)
                term1=[None for t in range(self.n+1)]
                term2=[None for t in range(self.n+1)]
                term3=[None for t in range(self.n+1)]
                term4=[None for t in range(self.n+1)]
                term1[self.n]=1/float(self.a)
                term2[self.n]=1/float(difpa)
                term3[0]=expmat
                term4[0]=expmmaxx
                s=(term1[self.n]-term2[self.n])*(term3[0]-term4[0])
                for q in range(1,self.n+1):
                    term1[self.n-q]=term1[self.n-q+1]/float(self.a)*(q+1)
                    term2[self.n-q]=term2[self.n-q+1]/float(difpa)*(q+1)
                    term3[q]=term3[q-1]*tau
                    term4[q]=term3[q-1]*maxx
                    s=s+(term1[self.n-q]-term2[self.n-q])*(term3[q]-term4[q])
                df=df+ar*ar/float(dif)*s
            else: # self.a==0
                relx=max(u/float(tau),1) 
                df=df+ar*(tau**(self.n+1))/float(self.n+1)*(relx**(self.n+1)-1)
                table=[None for t in range(self.n+1)]
                table[0]=dif*tau
                s=table[0]*(relx-1)
                for q in range(1,self.n+1):
                    table[q]=table[q-1]*dif*tau/float(q+1)
                    s=(s+table[q]*(relx**(q+1)-1))
                s=ar*ar*s/float(dif*dif)
                for q in range(1,self.n+1):
                    s=s*q/float(dif)
                df=df+s
            return  [self.k*df,self.k*deriv_at_0*u]
    def interval_f(self,x,tau):
        if x<=tau:
            return(self.f(x))
        else:
            return(0)
    def interval_vf_expst_classicqueue(self,u,tau,ar,omega):
        if not self.a==0:
            print 'error'
        else:
            if self.k==0:
                return [0,0]
            else:
                dif=float(omega-ar)
                x=min(u,tau)
                # computation of the derivative at zero of the drain function
                deriv_at_0=(self.n==0)*ar
                s=0
                coef=exp(-dif*tau)
                for t in range(self.n+1):
                    s+=(t==0)-coef*(dif*tau)**t/float(factorial(t))          
                deriv_at_0+=ar*ar*factorial(self.n)/float(dif**(self.n+1))*s
                # computation of the drain function
                df=ar*(x**(self.n+1))/float(self.n+1)
                coef=(exp(dif*x)-1)/float(exp(dif*tau))
                s=0
                for t in range(self.n+1):
                    s+=(dif*x)**(t+1)/float(factorial(t+1)) - coef*(dif*tau)**t/float(factorial(t))
                df=df+ar*ar*factorial(self.n)/float(dif**(self.n+2))*s
                #
                #return self.k*(df-deriv_at_0*u)
                return [self.k*df,self.k*deriv_at_0*u]
                
class realexponentialmixture:
    def __init__(self,v=None):
        if v is None:
            self.components=[]
            self.n=0
        else:
            self.n=len(v)
            self.components=[None for i in range(self.n)]
            for i in range(self.n):
                self.components[i]=realexponentialcomponent(v[i][0],v[i][1],v[i][2])
    def show(self):
        ans=[None for i in range(self.n)]
        for i in range(self.n):
            ans[i]=self.components[i].show()
        return ans
    def f(self,x):
        f=0
        for i in range(self.n):
            f+=self.components[i].f(x)
        return(f)
    def latex_f(self):
        text=''
        for i in range(self.n):
            text+=self.components[i].latex_f()
        return text 
    def largevalues_f(self,x,tau):
        f=0
        for i in range(self.n):
            f+=self.components[i].largevalues_f(x,tau)
        return(f)
    def largevalues_vf_expst_classicqueue(self,u,tau,ar,omega):
        vf=[0,0]
        for i in range(self.n):
            interv=self.components[i].largevalues_vf_expst_classicqueue(u,tau,ar,omega)
            vf[0]+=interv[0]
            vf[1]+=interv[1]
        return(vf)
    def interval_f(self,x,tau):
        f=0
        for i in range(self.n):
            f+=self.components[i].interval_f(x,tau)
        return(f)
    def interval_vf_expst_classicqueue(self,u,tau,ar,omega):
        vf=[0,0]
        for i in range(self.n):
            interv=self.components[i].interval_vf_expst_classicqueue(u,tau,ar,omega)
            vf[0]+=interv[0]
            vf[1]+=interv[1]
        return(vf)

                
class approximationinterval:
    bernstein_approximation=None
    korovkin_approximation=None    
    def __init__(self,n=None,tau=None): # initialize with maximum n so that korovkin is built only once
        if not tau is None:
            self.tau=tau
        if not n is None:
            self.format(n,tau)
        self.largevalues_interval=[realexponentialmixture(),realexponentialmixture()] # interval for large backlog values
        self.smallvalues=None # polynomial approximation for small backlog values
        self.smallvalues_interval=[realexponentialcomponent(),realexponentialcomponent()] # uniform error interval for small backlog values
        self.approx_init(n,tau)
    def adjust_largevalues_bounds(self,tau):
        pass
    def approx_init(self,n,tau):
        if not tau is None:
            self.tau=tau
            self.adjust_largevalues_bounds(tau)
        self.bernstein_approximation=bernstein(n,tau)
        self.korovkin_approximation=korovkin(n)
        self.korovkin_approximation.build(n)
    def format(self,n,tau=None):
        self.n=n
        if not tau is None:
            self.tau=tau
        if not self.bernstein_approximation is None:     
            self.bernstein_approximation.format(n,tau)
        if not self.korovkin_approximation is None:     
            self.korovkin_approximation.format(n)
        self.coefs=[None for i in range(n+1)]
        self.smallvalues=[None for i in range(n+1)]
        for i in range(n+1):
            self.smallvalues[i]=realexponentialcomponent(None,i,0)
    def set_coef(self,i,coef):
        if i<=self.n and i>=0:
            self.coefs[i]=coef
            self.smallvalues[i].k=coef
    def show_bernstein(self):
        print self.bernstein_approximation.beta,'-',self.largevalues.show()
    def compute_bernstein(self):
        self.bernstein_approximation.compute_beta()    
    def build_korovkin(self,n): # use first build with with n_max, computes all coefficients for n<=n_max
        self.korovkin_approximation.format(n) 
    def compute_korovkin(self,tau,n): # use first build with with n_max, computes all coefficients for n<=n_max
        #self.korovkin_approximation.format(n) # then this does not build anymore
        self.adjust_largevalues_bounds(tau)
        self.korovkin_approximation.compute_all(tau) 
    def make_korovkin_approximation(self,n,true_inaccessible_n=None):
        if true_inaccessible_n is None:
            true_inaccessible_n=n
        for k in range(n+1):
            self.set_coef(k,self.korovkin_approximation.coefficient(n,k))
            eta=self.korovkin_approximation.errorbound(true_inaccessible_n)
            self.smallvalues_interval=[realexponentialcomponent(-eta,0,0),realexponentialcomponent(eta,0,0)]
    def f_interval(self,u):
        fi=[0,0]
        # add each exact component from small values
        for k in range(self.n+1):
            f=self.smallvalues[k].interval_f(u,self.korovkin_approximation.tau)
            fi[0]+=f
            fi[1]+=f
        # add uniform interval from small values
        f=self.smallvalues_interval[1].interval_f(u,self.korovkin_approximation.tau)
        fi[0]-=f
        fi[1]+=f
        # add (exponential) interval from large values
        fi[0]+=self.largevalues_interval[0].largevalues_f(u,self.korovkin_approximation.tau)
        fi[1]+=self.largevalues_interval[1].largevalues_f(u,self.korovkin_approximation.tau)
        return fi 
    def latex_f_interval(self):
        f=self.smallvalues_interval[1].latex_f()
        text_interval=['-'+f,f]
        for k in range(self.n+1):
            f=self.smallvalues[k].latex_f()
            text_interval[0]+=' + '+f
            text_interval[1]+=' + '+f
        text_largevalues_interval=[self.largevalues_interval[0].latex_f(),self.largevalues_interval[1].latex_f()]
        return [text_interval,text_largevalues_interval]   
    # clearance function interval       
    def df_interval(self,u,ar,omega):
        dfi=[0,0]
        # add each exact component from small values
        for k in range(self.n+1):
            [df,rate]=self.smallvalues[k].interval_vf_expst_classicqueue(u,self.korovkin_approximation.tau,ar,omega)
            dfi[0]+=(df-rate)
            dfi[1]+=(df-rate)
        # add uniform interval from small values
        [df,rate]=self.smallvalues_interval[1].interval_vf_expst_classicqueue(u,self.korovkin_approximation.tau,ar,omega)
        #print dfi[0],"-",(df+rate)
        #print dfi[1],"+",(df+rate)
        dfi[0]-=(df+rate)
        dfi[1]+=(df+rate)
        if (df<0 or rate<0):
            print 'uniform error negative:',df,rate
        # add (exponential) interval from large values
        [dfl,ratel]=self.largevalues_interval[0].largevalues_vf_expst_classicqueue(u,self.korovkin_approximation.tau,ar,omega)
        [dfu,rateu]=self.largevalues_interval[1].largevalues_vf_expst_classicqueue(u,self.korovkin_approximation.tau,ar,omega)
        #print dfi[0],"-",(df+rate)
        #print dfi[1],"+",(df+rate)
        dfi[0]+=(dfl-rateu)
        dfi[1]+=(dfu-ratel) 
        if (rateu<ratel):
            print 'rateu=',rateu, '<', 'ratel=',ratel      
        if (dfu<dfl):
            print 'dfu=',dfu, '<', 'dfl=',dfl      
        return dfi
    def admission_cost_interval(self,u,ar,omega,st):
        cost=self.f(u)
        df=self.df_interval(u,ar,omega)
        dfx=self.df_interval(u+st,ar,omega)
        return([cost+dfx[0]-df[1],cost+dfx[1]-df[0]])
    def f(self,x):
        return(0)
    def korovkin_latex(self,ar=1.0,omega=2.0,nmin=1,nmax=15,maxn=15,points=1000,xr=[0,3],yr=[-1,1],colmin=0.2,colmax=0.95,savepath=None,filename='korovkin-graphs',latexaccesspath='figs/'):
        interval_points=int(min(points,points*(self.tau-xr[0])/(xr[1]-xr[0]))+1)
        interval_step=(self.tau-xr[0])/float(interval_points-1)
        x=range(interval_points+1)
        for p in x:
            x[p]=x[p]/float(interval_points)*(self.tau-xr[0])
        if xr[1]>self.tau:
            largevalues_points=int(points-(interval_points-1)+1)
            largevalue_step=(xr[1]-self.tau)/float(largevalues_points-1)
            y=range(largevalues_points+1)
            for p in y:
                y[p]=self.tau+y[p]/float(largevalues_points)*(xr[1]-self.tau)
            x=x+y
        filename=filename+'-'+str(self.tau)+'-'+str(nmax)
        if not savepath is None:
            #f#cost=file(savepath+filename+'-costf.tex','w')
            #fcostU=file(savepath+filename+'-Ucostf.tex','w')
            #fcostE=file(savepath+filename+'-Ecostf.tex','w')
            #vfU=file(savepath+filename+'-vfU.tex','w')
            #vfL=file(savepath+filename+'-vfL.tex','w')
            fcost=file(savepath+filename+'-fcost.tex','w')
            vf=file(savepath+filename+'-vf.tex','w')
            vftablel=[file(savepath+filename+'-vftablel-'+str(k)+'.tex','w') for k in range(nmin,nmax+1)]
            vftableu=[file(savepath+filename+'-vftableu'+str(k)+'.tex','w') for k in range(nmin,nmax+1)]
        else:
            f=file(filename+'.tex','w')
        fcost.write('\\addplot[id=costfunction,domain={\\xlim}:{\\xlimm}, samples={\\nbsamples}, {\\solstyle}, smooth,opacity=1.0,line width={\\solwidth}, mark = none, mark options={scale=0.1}, color = black] { '+self.latex_f()+' };\n%\n')
        self.compute_korovkin(self.tau,min(nmax,maxn))
        nstep=1.0/float(nmax-nmin+1)
        i=0
        for k in range(nmin,nmax+1):
            col=str(colmin+(colmax-colmin)*(nmax-k)*nstep)
            self.format(min(k,maxn))
            self.make_korovkin_approximation(min(k,maxn),k)
            #
            lfi=self.latex_f_interval()
            fcost.write('\\addplot[id=costfunctioniE,domain={\\xlim}:{'+str(self.tau)+'}, samples={\\nbsamplesinterval}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] { '+lfi[1][0]+' };\n%\n')
            fcost.write('\\addplot[id=costfunctioniU,domain={\\xlim}:{'+str(self.tau)+'}, samples={\\nbsamplesinterval}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1},  color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] { '+lfi[1][1]+' };\n%\n')
            fcost.write('\\addplot[id=costfunctionlvE,domain={'+str(self.tau)+'}:{\\xlimm}, samples={\\nbsampleslargevalues}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1},  color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] { '+lfi[0][0]+' };\n%\n')
            fcost.write('\\addplot[id=costfunctionlvU,domain={'+str(self.tau)+'}:{\\xlimm}, samples={\\nbsampleslargevalues}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1},  color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] { '+lfi[0][1]+' };\n%\n')
            #
            curveu="{:>15}  {:>15}\n".format('x','y')
            curvel="{:>15}  {:>15}\n".format('x','y')
            for p in range(len(x)):
                interval=self.df_interval(x[p],ar,omega)
                curveu+="{:>15}  {:>15}\n".format(x[p],interval[1])      
                curvel+="{:>15}  {:>15}\n".format(x[p],interval[0])      
            vftableu[i].write(curveu)
            vftableu[i].close()
            vftablel[i].write(curvel)
            vftablel[i].close()
            vf.write('\\addplot[id=valuefunctionl,domain={\\xlim}:{\\xlimm}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] table{'+latexaccesspath+filename+'-vftablel-'+str(k)+'.tex};\n%\n')
            vf.write('\\addplot[id=valuefunctionl,domain={\\xlim}:{\\xlimm}, smooth,opacity=1.0,line width={\\curvewidth}, mark = none, mark options={scale=0.1}, color = {lightcol!'+str(float(int((nmax-k)*nstep*1000)/float(10)))+'!darkcol}] table{'+latexaccesspath+filename+'-vftableu'+str(k)+'.tex};\n%\n')
            i+=1
        fcost.close()
        vf.close()
        print 'Ok.'


            
        
def intervalsum(l1,l2):
    return([l1[0]+l2[0],l1[1]+l2[1]])    
def intervaldifference(l1,l2):
    return([l1[0]-l2[1],l1[1]-l2[0]])    
def intervalnosmaller(l1,l2):
    return(l1[0]>=l2[1])   
def intervalnogreater(l1,l2):
    return(l1[1]<=l2[0])   
        
class sin_bernstein(bernstein):
    scale=3
    def f(self,x):
        return sin(x*self.scale)
    def continuitymodulus(self,delta,tau):
        if (delta*self.scale<pi):
            return sin(min(pi/2,delta*self.scale))
        else:
            return 1-sin(min(3*pi/2,delta*self.scale))    
    

class sin_approximationinterval(approximationinterval):
    def approx_init(self,n,tau):
        self.bernstein_approximation=sin_bernstein(tau,n)
        self.korovkin_approximation=korovkin(n)

        
class xsquaredoverasquaredplusxsquared_approximationinterval(approximationinterval):
    def approx_init(self,n=None,tau=None):
        self.bernstein_approximation=xsquaredoverasquaredplusxsquared_bernstein(tau,n)
        self.korovkin_approximation=xsquaredoverasquaredplusxsquared_korovkin(n)
                
                                
class xsquaredoverasquaredplusxsquared_bernstein(bernstein):
    a=None
    def set_a(self,a):
        self.a=a
    def f(self,x):
        if (x>0):
            return 1/float(float(self.a/float(x))**2+1)
        else:
            return 0
    def continuitymodulus(self,delta,tau):
        hds=0.5*delta*delta
        asq=self.a*self.a
        xstar=max(0.5*delta,min(tau-0.5*delta,sqrt((hds-asq+2*sqrt(asq*asq+asq*hds+hds*hds))/float(3))))
        return(self.f(xstar+0.5*delta)-self.f(xstar-0.5*delta))

class xsquaredoverasquaredplusxsquared_korovkin(korovkin):
    a=None
    def set_a(self,a):
        self.a=a
    def f(self,x):
        if (x>0):
            return 1/float(float(self.a/float(x))**2+1)
        else:
            return 0
    def continuitymodulus(self,delta,tau):
        hds=0.5*delta*delta
        asq=self.a*self.a
        xstar=max(0.5*delta,min(tau-0.5*delta,sqrt((hds-asq+2*sqrt(asq*asq+asq*hds+hds*hds))/float(3))))
        return(self.f(xstar+0.5*delta)-self.f(xstar-0.5*delta))
    def compute_fourier_coefficients(self,tau):
        if self.a is None:
            print "Parameter a not yet defined."
        else:
            self.tau=tau
            self.laurent=[None for k in range(self.n+1)]
            self.kappa=[None for k in range(self.n+1)]
            self.firstterm=[None for k in range(self.n+1)]
            self.residues=[None for k in range(self.n+1)]
            self.fourier_beta_coef=[None for k in range(self.n+1)]
            self.unnamed=[None for k in range(self.n+1)]
            self.sumunnamed=[None for k in range(self.n+1)]
            factor=-2/float(self.tau)
            scaledfactor=self.a*factor
            scaledfactorsquared=-factor*factor*self.a*self.a
            theta=atan(self.tau/float(self.a))
            coshalftheta=cos(0.5*theta)
            sinhalftheta=sin(0.5*theta)
            rescoef=pow(self.a**2/float(self.a**2+self.tau**2),0.25)
            for k in range(self.n+1):
                # Laurent coefficient
                unnamed=[None for q in range(k+1)]
                sumunnamed=[None for q in range(k+1)]
                f=1.0                
                self.laurent[k]=[None for q in range(k+1)]
                for q in range(k+1):
                    self.laurent[k][q]=0
                    unnamed[q]=[None for j in range((k-q)/2+1)]
                    sumunnamed[q]=[0 for j in range((k-q)/2+1)]
                    if self.is_even(k):
                        value=-2*((-q)/2)
                        incr=0
                    else:
                        value=-2*((1-q)/2)
                        incr=1
                    for j in range(len(unnamed[q])):
                        unnamed[q][j]=[None for t in range(j+1)]
                        sumunnamed[q][j]=0
                        fs=1.0
                        for t in range(j+1):
                            if not t:
                                # initialisation of first element
                                if not (q or j):
                                    unnamed[0][0][0]=1
                                else:
                                    if not j:
                                        if self.is_even(q)^self.is_even(k):
                                            unnamed[q][0][0]=q+1
                                        else:
                                            unnamed[q][0][0]=1
                                    else:
                                        unnamed[q][j][0]=unnamed[q][j-1][0]*(value+2*j+incr)*(value+2*j+incr-1)/((value+2*j+incr-q)*(value+2*j+incr-q-1))                      
                            else:
                                # t=0 initialized
                                x=value+2*j+incr
                                qpdt=q+2*t
                                unnamed[q][j][t]=unnamed[q][j][t-1]*(x-qpdt+2)*(x-qpdt+1)/((qpdt)*(qpdt-1))#*scaledfactorsquared
                            # unnamed[q][j][t] ready
                            sumunnamed[q][j]+=unnamed[q][j][t]*fs
                            fs*=scaledfactorsquared
                        # sumunnamed[q][j] ready
                        self.laurent[k][q]+=sumunnamed[q][j]*self.nu[k][j+value/2]
                    if self.is_even(k):
                        self.laurent[k][q]*=f
                    else:
                        self.laurent[k][q]*=-f
                    f*=factor
                #print k,unnamed
                self.unnamed[k]=unnamed
                self.sumunnamed[k]=sumunnamed
                #
                # residues
                self.kappa[k]=[None for q in range(k/2+1)]
                for q in range(k/2+1):
                    if self.is_even(k):
                        realvector=[None for t in range(q+1)]
                        realvector[0]=1
                        for t in range(1,q+1):
                            realvector[t]=realvector[t-1]*(2*q-2*t+1)*(2*q-2*t+2)/((2*t)*(2*t-1))
                        imaginaryvector=[None for t in range(q)]
                        if q:
                            imaginaryvector[0]=2*q
                            for t in range(1,q):
                                imaginaryvector[t]=imaginaryvector[t-1]*(2*q-2*t)*(2*q-2*t+1)/((2*t+1)*(2*t))
                    else:
                        realvector=[None for t in range(q+1)]
                        realvector[0]=1
                        for t in range(1,q+1):
                            realvector[t]=realvector[t-1]*(2*q-2*t+2)*(2*q-2*t+3)/((2*t)*(2*t-1))
                        imaginaryvector=[None for t in range(q+1)]
                        imaginaryvector[0]=2*q+1
                        for t in range(1,q+1):
                            imaginaryvector[t]=imaginaryvector[t-1]*(2*q-2*t+1)*(2*q-2*t+2)/((2*t+1)*(2*t))
                    #print realvector
                    #print imaginaryvector
                    # realvector,imaginaryvector ready
                    fs=1
                    realterm=0
                    for t in range(len(realvector)):
                        realterm+=realvector[t]*fs
                        fs*=scaledfactorsquared
                    fs=1
                    imaginaryterm=0
                    for t in range(len(imaginaryvector)):
                        imaginaryterm+=imaginaryvector[t]*fs
                        fs*=scaledfactorsquared
                    if self.is_even(k):
                        self.kappa[k][q]=coshalftheta*realterm-sinhalftheta*(-scaledfactor)*imaginaryterm
                    else:
                        self.kappa[k][q]=-coshalftheta*realterm+sinhalftheta*(-scaledfactor)*imaginaryterm
                # self.kappa[k] ready
                #print self.kappa[k]
                self.residues[k]=0
                for q in range(k/2+1):
                    self.residues[k]+=self.kappa[k][q]*self.nu[k][q]
                self.residues[k]*=rescoef
                # computation of the first term
                self.firstterm[k]=0
                firsttermcoef=1.0
                for q in range(k+1):
                    self.firstterm[k]+=firsttermcoef*self.laurent[k][q]
                    firsttermcoef*=-(q+0.5)/float(q+1)*(-self.tau)#
                #
                self.fourier_beta_coef[k]=self.firstterm[k]-self.residues[k]
                if (k>0):
                    self.fourier_beta_coef[k]*=2                                               
            #
            self.fourier_coefs_ready=True
            #print self.firstterm
            #print self.residues
            #print self.fourier_beta_coef
                
class xsquaredoverasquaredplusxsquared_approximationinterval(approximationinterval):
    def approx_init(self,n,tau=None):
        self.a=None
        self.korovkin_approximation=xsquaredoverasquaredplusxsquared_korovkin(n)
    def f(self,x):
        return(self.korovkin_approximation.f(x))
    def set_a(self,a):
        self.a=a
        self.korovkin_approximation.set_a(a)
    def adjust_largevalues_bounds(self,tau):
        if not self.a is None:
            f=self.korovkin_approximation.f(tau)
            self.largevalues_interval=[realexponentialmixture([[f,0,0]]),realexponentialmixture([[1,0,0],[-(1-f)*exp(2*self.a*self.a*tau/float((1-f)*(self.a**2+tau**2)**2)*tau),0,2*self.a*self.a*tau/float((1-f)*(self.a**2+tau**2)**2)]])] # interval for large backlog values
    def latex_f(self):
        return '(x^2)/(('+str(self.a)+')^2+(x^2))'

    

                    
        
                   
