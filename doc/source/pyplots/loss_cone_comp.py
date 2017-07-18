import matplotlib.pyplot as plt
from numpy import *

def R_NR_BW77(g,x):
    Q = 4e6
    nh = 4e4
    rh = 1.9
    xD = 1e5
    alpha = Q**(-2)*nh*rh**3*xD
    q = 1.6*log(6/pi*Q)*alpha*x**(-5./2.)*g
    return 0.5*g*x**(5./2.)/alpha/log(Q)/(5.56+log(xD/4./x)/q)

def R_NR_HA06b(g,x):
    x_pinhole = 0.5
    R = zeros(len(x))
    for i in xrange(len(x)):
        if (x[i] < x_pinhole):
            R[i] = 0.0
        else:
            R[i] = g[i]**2/log(sqrt(1e7/32./x[i])) 
    return R

def R_RR(g,x):
    Jc2 = 1e7/32./x
    tm = g*sqrt(8.0)*x**(3.0/2.0)
    tgr_inv =  3./8.*log(Jc2)/Jc2
    trr = 1/tm+tgr_inv
    return g/trr/log(4e6)/1.98**2

D = loadtxt('/home/ben/MasterProject/Code/FP_solver/loss_cone_tests/SM_NR_HA06A/FP_last.out')
x = D[:,1]
g = D[:,2]
plt.loglog(x,R_NR_BW77(g,x),label='$R_{NR}^{BW77}$')
plt.loglog(x,R_NR_HA06b(g,x)*0.5,label='$R_{NR}^{HA06b}/2$ (for single mass)')
plt.loglog(x,R_NR_BW77(g,x)+R_RR(g,x),label='$R_{NR}^{BW77}+R_{RR}$')
plt.loglog(x,R_NR_HA06b(g,x)/2+R_RR(g,x),label='$R_{NR}^{HA06b}/2+R_{RR}$ (for single mass)')
leg = plt.legend(loc = 4)
for t in leg.get_texts():
    t.set_fontsize('small')
plt.xlabel('$X$')
plt.ylabel('sink term')
plt.show()
