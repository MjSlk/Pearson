import numpy as np
from scipy.integrate import quad

def Q(x,a0,a1,a2,a3):  
    ### x: density bins
    ### a,a1,a2,a3: Pearson parameters
    return (x-a0)/(a1+a2*x+a3*x**2)


def moment(x,y,n):
    dx = np.diff(x)

    return np.sum(x[:-1]**n * y*dx)
    
    
def params_pearson(c,d,e,f):
##Fct to compute Pearson coefficients
###c,d,e,f: 4 first non-statistical moments of the data
    a0=c*(-13*d*f+8*e**2+9*d**3) + e*f + 12*(c**3)*f + 3*(d**2)*e - 20 *(c**2)*d*e
    a1=c*(e*f+ (d**2) *e) +(c**2) *(3*d*f -4*e**2) -4*(d**2)*f +3*d*(e**2)
    a2=c*(-7*d*f+2*e**2+ 3*d**3) +e*f +6*(c**3)*f+ 3*(d**2)*e -8*(c**2)*d*e
    a3=(c**2 )*(2*f-3*d**2)-2*d*f+3*e**2 -10*c*d*e +4 *(c**3)*e + 6*d**3
    DeN=c**2*(10*f-6*d**2)-10 *d*f +12*e**2 -32*c*d*e +8*(c**3)*e +18*d**3


    a0=a0/DeN
    a1=-a1/DeN
    a2=a2/DeN
    a3=-a3/DeN

    return a0,a1,a2,a3

def Pearson_fit(moments,delta,pdf):
    m1,m2,m3,m4 = moments
    x=delta[:-1] ##taking the bin starting values
    y=pdf
    logy=np.log(np.where(y==0,1e-7,y))   ##replacing 0 values with simulation resolution
    logder=np.gradient(logy,x)
    der=np.gradient(y,x)  ##derivative of the PDF of the data

    a0,a1,a2,a3= params_pearson(m1,m2,m3,m4)
    
    Delta=a2**2-4*a3*a1 ##discriminant of the denominator in Pearson differential equation
    x1=(-a2-np.sqrt(Delta))/(2*a3)
    x2=(-a2+np.sqrt(Delta))/(2*a3)

    # eps=(max(x1,x2)-min(x1,x2))/100
    eps=0.001
    Q_fit=Q(x, a0, a1, a2, a3) 
    
###Applying domain of applicability conditions
    if Delta>0:
        if a3>0:
            dPmin=min(x1,x2)
            dPmax=max(x1,x2)
            di,df= max(dPmin+eps, x.min()), min(dPmax-eps,x.max())


        else:
            if x.max()<=max(x1,x2) and x.min()<=min(x1,x2):
                dPmin=-np.inf
                dPmax=min(x1,x2)
                di,df=x.min(), min(x.max(),min(x1,x2)-eps)

            elif x.min()>=min(x1,x2) and x.max()>=max(x1,x2):
                dPmin=max(x1,x2)
                dPmax=np.inf
                di,df= max(x.min(),max(x1,x2)+eps),x.max()
            else:
                dPmin= np.nan
                dPmax= np.nan
                di,df=x.min(),x.max()
    else:
        if a3<=0:
            dPmin=-np.inf
            dPmax=np.inf
            di,df=x.min(),x.max()
        else:
            dPmin= np.nan
            dPmax= np.nan
            di,df=x.min(),x.max()
  
    x_fit=x[np.where(x>=di)] ##chosing bins which are inside the domain of applicability (practically df is always very large, no need to include it in the condition)
    y_=y[np.where(x>=di)]
    phi=np.zeros(len(x_fit))
    for i in range(len(phi)):
        phi[i]+=quad(Q,di,min(x_fit[i],df),args=(a0,a1,a2,a3))[0] ###Integratubg the differential equation
    #C=moment(x_fit,y_,0)/moment(x_fit, np.exp(phi), 0) ## Assuring the integrated fit and the data have the same normalistion within the domain of applicability
    x_norm = delta[np.where(x>=di)]
    C=moment(delta,y_,0)/moment(delta, np.exp(phi), 0) ## Assuring the integrated fit and the data have the same normalistion within the domain of applicability
    logy_fit=phi+np.log(C)
    y_fit=C*np.exp(phi) ##Taking the exponential (Pearson diff. equation is for the logarithmic derivative)
    #print(dPmin,dPmax)
    return x_fit, Q_fit, y_fit ####y_fit is Pearson fit (x_fit can be useful in the analysis if domain of applicability does not include all the density bins)
    
    
def D_KL(x,P,Q): ###fct to compute KL divergence between 2 continuous PDFs.
    ### x: density bins 
    ### P: data PDF
    ### Q: Fitted PDF
    dx = np.diff(x)
    return np.sum(P * np.log(P/Q)*dx)

