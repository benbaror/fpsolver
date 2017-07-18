from pylab import *
def get_Trr(Cm,M,g,x,rh2rg,N,e_flag=1.0):
    tm = sqrt(8)*sum(Cm*M)*x**(3./2.)/sum(g*M)
    Jc2 = rh2rg/32/x
    J2 =  linspace(1,Jc2,N)[:-1]
    tgr = 8/3.*J2
    ecc2 = 1-J2/Jc2
    if e_flag:
        trr = sum(abs(1/tm-1/tgr)*J2/(pi**2*0.25*ecc2))*abs(J2[1]-J2[0])/Jc2**2
    else:
        trr = sum(abs(1/tm-1/tgr)*J2)*abs(J2[1]-J2[0])/Jc2**2
    return trr

def get_Trr2(Cm,M,g,x,rh2rg,N,e_flag=1.0):
    tm = sqrt(8)*sum(Cm*M)*x**(3./2.)/sum(g*M)
    Jc2 = rh2rg/32/x
    Jc = sqrt(Jc2)
    J =  (linspace(1,Jc,N))[:-1]
    J2 = J**2
    tgr = 8/3.*J2
    ecc2 = 1-J2/Jc2
    if e_flag:
        trr = sum(abs(1/tm-1/tgr)*J2/(pi**2*0.25*ecc2))*abs(J[1]-J[0])/Jc**3
    else:
        trr = sum(abs(1/tm-1/tgr)*J2)*abs(J[1]-J[0])/Jc**3
    return trr

def get_trr(Cm,M,g,rh2rg,X,N,e_flag=1.0):
    trr = zeros(len(X))
    for i in xrange(len(X)):
        trr[i] = get_Trr(Cm,M,g[i,:],X[i],rh2rg,N,e_flag=e_flag)
    return trr

def get_trr2(Cm,M,g,rh2rg,X,N,e_flag=1.0):
    trr = zeros(len(X))
    for i in xrange(len(X)):
        trr[i] = get_Trr2(Cm,M,g[i,:],X[i],rh2rg,N,e_flag=e_flag)
    return trr




def get_Trr3(a,N):
    G =  6.6720e-8
    c = 2.99792458e10
    M0 = 1.9891e33
    pc = 3.08568e18
    yr = 31536000.0
    Mbh = 3e6
    Ms  = 1.
    Q   = Mbh/Ms
    rh = 2*sqrt(Mbh/3e6) 
    alpha = 1.75
    Nh = 4.*pi*4e4*rh**3/(3.-alpha)/sqrt(Mbh/3e6)/Ms
    Jc2 = G*M0*Mbh*a*pc
    J2 = linspace(0,Jc2,N)[:-1]
    ecc2 = 1-J2/Jc2
    P = 2*pi*(a*pc)**(3/2.)/sqrt(G*Mbh*M0)
#    J = sqrt(G*Mbh*M0*a*pc*(1-ecc**2))
    Jlso2 = (4*G*Mbh*M0/c)**2
    N = Nh*(a/2)**(3-1.75)
    TM = Q/N*P
    TGR = 8/3.*(J2/Jlso2)*P
#    TRR = (0.53)**-2*Q**2/N*P**2*abs(sqrt(8)/TM-1.0/TGR)/(pi**2*0.25*ecc2)
    TRR = Q**2/N*P**2*abs(sqrt(8)/TM-1.0/TGR)
#    TRR2 = TRR[argmin(abs(ecc-0.7))]/yr/1e6
#    ecc_r = rand(100000)
#    J_r = sqrt(G*Mbh*M0*a*pc*(1-ecc_r**2))
#    TGR_r = 8/3.*(J_r/Jlso)**2*P
#    TRR3 = mean((0.53)**-2*Q**2/N*P**2*abs(1/TM-1/TGR_r))/yr/1e6
    TNR = Q**2*P/N/log(Q)/yr/1e6
    TRR = sum(TRR[1:])*abs(J2[0]-J2[1])/Jc2/yr/1e6
    return TRR,TNR


def R_NR(g,x):
    Q = 4e6
    nh = 4e4
    rh = 1.9
    xD = 1e5
    alpha = Q**(-2)*nh*rh**3*xD
    q = 1.6*log(6/pi*Q)*alpha*x**(-5./2.)*g
    return 0.5*g*x**(5./2.)/alpha/log(Q)/(5.56+log(xD/4./x)/q)






D = loadtxt('../MS6/FP_last.out')
X = D[:,1]
g = D[:,[2,5,8,11]]
Cm = array([1.000E-01,1.000E+00,1.000E-02,1.000E-03])
M = array([6.000E-01,1.000E+00,1.400E+00,1.000E+01])

loglog(X,1/(1+1/get_trr(Cm,M,g,1e7,X,100000,e_flag=0)/log(4e6)),label='iso')
loglog(X,1/(1+1/get_trr(Cm,M,g,1e7,X,100000,e_flag=1)/log(4e6)),label='iso e dep')
loglog(X,1/(1+1/get_trr2(Cm,M,g,1e7,X,100000,e_flag=0)/log(4e6)),label='uni')
loglog(X,1/(1+1/get_trr2(Cm,M,g,1e7,X,100000,e_flag=1)/log(4e6)),label='uni e dep')
#loglog(1.9/2.0/a,1/(1+TNR/TRR),label='single mass no e dep')
