import matplotlib.pyplot as plt
from numpy import *
from FP import *    


path='/home/ben/MasterProject/Code/FP_solver/Runs/'
dirs_NR = ['SM_NR_HA06a',
           'SM_NR_AH09', # AH09 loss cone term  with X_pinhole = 0.5
           'SM_NR_AH09_E09']
dirs_RR = ['SM_RR_HA06a',
           'SM_RR_AH09',
           'SM_RR_AH09_E09']
labels = ['HA06a',
          'AH09',
          'AH09 with E09']
colors = ['r','b','g','k']
FP_NR = [FP(path+dirs_NR[i],labels[i],1) for i in xrange(len(dirs_NR))]
FP_RR = [FP(path+dirs_RR[i],labels[i],1) for i in xrange(len(dirs_NR))]
for i in xrange(len(FP_RR)):
    FP_RR[i].linestyles = '--'
    FP_RR[i].labels = ''
    FP_NR[i].colors = colors[i]
    FP_RR[i].colors = colors[i]
    FP_RR[i].FP_NR = FP_NR[i]
#############################################
for i in xrange(len(FP_NR)):
    plt.figure(1)
    FP_NR[i].Plot_Q()
    FP_RR[i].Plot_Q()
    plt.figure(2)
    FP_NR[i].Plot_g()
    FP_RR[i].Plot_g()
    # plt.figure()
    # FP_NR[i].Plot_n()
    # FP_RR[i].Plot_n()
    # plt.figure()
    # FP_NR[i].Plot_N()
    # FP_RR[i].Plot_N()
    plt.show()
