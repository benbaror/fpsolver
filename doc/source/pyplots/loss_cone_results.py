import matplotlib.pyplot as plt
from numpy import *
path='/home/ben/MasterProject/Code/FP_solver/loss_cone_tests/'
dirs_NR = ['SM_NR_HA06A', # BW77 loss cone term
           'SM_NR_HA06b', # HA06b loss cone term multply by half with X_pinhole = 10
           'SM_NR_HA06b_xp_0.5'] # HA06b loss cone term multply by half with X_pinhole = 0.5
dirs_RR = ['SM_RR_HA06A2','SM_RR_HA06b','SM_RR_HA06b_xp_0.5']
QNR = []
for d in dirs_NR:
    QNR.append(loadtxt(path + d + '/Q_last.out'))

QRR = []
for d in dirs_RR:
    QRR.append(loadtxt(path + d + '/Q_last.out'))

gNR = []
for d in dirs_NR:
    gNR.append(loadtxt(path + d + '/FP_last.out'))

gRR = []
for d in dirs_RR:
    gRR.append(loadtxt(path + d + '/FP_last.out'))

nNR = []
for d in dirs_NR:
    nNR.append(loadtxt(path + d + '/density.out'))

nRR = []
for d in dirs_RR:
    nRR.append(loadtxt(path + d + '/density.out'))

NNR = []
for d in dirs_NR:
    NNR.append(loadtxt(path + d + '/number.out'))

NRR = []
for d in dirs_RR:
    NRR.append(loadtxt(path + d + '/number.out'))


labels = ['HA06a',
          'HA06b x_pinhole=10',
          'HA06b x_pinhole=0.5']
colors = ['r','b','g']
#############################################
plt.figure()
for i in xrange(len(QNR)):
    plt.loglog(QNR[i][:,1],QNR[i][:,5],label=labels[i],color=colors[i])    

for i in xrange(len(QRR)):
    plt.loglog(QRR[i][:,1],QRR[i][:,5],'--',color=colors[i])    

for i in xrange(len(QRR)):
    plt.loglog(QRR[i][:,1],QRR[i][:,5]/QNR[i][:,5],color=colors[i])    

leg = plt.legend(loc = 3)
for t in leg.get_texts():
    t.set_fontsize('small')
plt.xlim([2,1e3])
plt.ylim([1e-3,10])
plt.xlabel('$X$')
plt.ylabel('$Q$')
############################################
plt.figure()
for i in xrange(len(gNR)):
    plt.loglog(gNR[i][:,1],gNR[i][:,2],label=labels[i],color=colors[i])    

for i in xrange(len(gRR)):
    plt.loglog(gRR[i][:,1],gRR[i][:,2],'--',color=colors[i])    

leg = plt.legend(loc = 2)
for t in leg.get_texts():
    t.set_fontsize('small')
plt.xlim([1e-1,2e4])
plt.ylim([1,30])
plt.xlabel('$X$')
plt.ylabel('$g$')
################################################
plt.figure()
for i in xrange(len(nNR)):
    plt.loglog(nNR[i][:,0],nNR[i][:,1],label=labels[i],color=colors[i])    

for i in xrange(len(nRR)):
    plt.loglog(nRR[i][:,0],nRR[i][:,1],'--',color=colors[i])    

leg = plt.legend(loc = 3)
for t in leg.get_texts():
    t.set_fontsize('small')
plt.xlim([1e-4,2])
plt.ylim([2e4,2e12])
plt.xlabel('$r[pc]$')
plt.ylabel('$n[pc^{-3}]$')

################################################
plt.figure()
for i in xrange(len(NNR)):
    plt.loglog(NNR[i][:,0],NNR[i][:,1],label=labels[i],color=colors[i])    

for i in xrange(len(NRR)):
    plt.loglog(NRR[i][:,0],NRR[i][:,1],'--',color=colors[i])    

leg = plt.legend(loc = 2)
for t in leg.get_texts():
    t.set_fontsize('small')
plt.xlim([1e-4,2])
plt.ylim([1,4e6])
plt.xlabel('$r[pc]$')
plt.ylabel('$N$')

###############################3
plt.show()

