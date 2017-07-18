import matplotlib.pyplot as plt
from numpy import *
from FP import *    


path='/home/ben/MasterProject/Code/FP_solver/Runs/'
dirs_NR = ['MS_NR_AH09', # AH09 loss cone term  with X_pinhole = 0.5
           'MS_NR_AH09_E09']
dirs_RR = ['MS_RR_AH09',
           'MS_RR_AH09_E09']
labels = ['AH09',
          'AH09 with E09']
colors = ['r','b','g','k']
FP_NR = [FP(path+dirs_NR[i],labels[i],4) for i in xrange(len(dirs_NR))]
FP_RR = [FP(path+dirs_RR[i],labels[i],4) for i in xrange(len(dirs_NR))]
for i in xrange(len(FP_RR)):
    FP_RR[i].linestyles = ['--','--','--','--']
    FP_RR[i].labels = ['','','','']
    FP_NR[i].FP_RR = FP_RR[i]
#############################################
def Plot_Q(FP_NR,FP_RR,colors):
    """
    
    Arguments:
    - `FP_NR`:
    - `FP_RR`:
    """
    plt.figure()
    if (FP_NR.nM == 1):
        plt.loglog(FP_NR.Q_x,FP_NR.Q,label=FP_NR.label,color=colors[i])
        plt.loglog(FP_RR.Q_x,FP_RR.Q,'--',color=colors[i])    
        plt.loglog(FP_RR.Q_x,FP_RR.Q/FP_NR.Q,color=colors[i])    
    else:
        for i in xrange(FP_NR.nM):
            plt.loglog(FP_NR.Q_x,FP_NR.Q[:,i],label=FP_NR.label,color=colors[i])
            plt.loglog(FP_RR.Q_x,FP_RR.Q[:,i],'--',color=colors[i])    
            plt.loglog(FP_RR.Q_x,FP_RR.Q[:,i]/FP_NR.Q[:,i],color=colors[i])    

    leg = plt.legend(loc = 3)
    for t in leg.get_texts():
        t.set_fontsize('small')
    plt.xlim([2,1e3])
    #plt.ylim([1e-3,10])
    plt.xlabel('$X$')
    plt.ylabel('$Q$')
# ############################################

#Plot_Q(FP_NR[0],FP_RR[0],colors)
for i in xrange(len(FP_NR)):
    plt.figure()
    FP_NR[i].Plot_Q()
    FP_RR[i].Plot_Q()
    plt.figure()
    FP_NR[i].Plot_g()
    FP_RR[i].Plot_g()
    plt.figure()
    FP_NR[i].Plot_n()
    FP_RR[i].Plot_n()
    plt.figure()
    FP_NR[i].Plot_N()
    FP_RR[i].Plot_N()
    plt.show()
