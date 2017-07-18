import os
import numpy as np
####### Submit NR
name = 'SM_NR_AH09_E09'
command = 'python Submit_FP.py ' + name +  ' -E SM_NR_AH09_E09'
os.system(command)

name = lambda ARR,AGR : 'SM_RR_AH09_E09_ARR_AGR_%2.2f_%2.2f' % (ARR, AGR) 
command = lambda ARR,AGR: 'python Submit_FP.py ' + name(ARR,AGR) +  ' -E SM_RR_AH09_E09 --ARR '+ str(ARR) + ' --AGR ' + str(AGR)

#ARR = [0.05,0.1,0.2,0.4,0.6,0.8,1.0,2.0,4.0,8.0,10.0,50.0,100.0]
#AGR = [0.2,0.4,0.6,0.8,1.0,2.0,4.0,8.0,10.0]
path = '/home/ben/mp/Data/FP_Solver/Runs/'

ARR = np.logspace(-1,1,10)
AGR = 1/ARR
#AGR = [0.2,0.4,0.6,0.8,1.0,2.0,4.0,8.0,10.0]
#AGR = [1.0]
path = '/home/ben/mp/Data/FP_Solver/Runs/'


dir_list = os.listdir(path)
print len(dir_list)
for Arr in ARR:
    for Agr in AGR:
        temp_path = path+name(Arr,Agr)+'/'
        if name(Arr,Agr) in dir_list:
            file_list = os.listdir(temp_path)
            if 'FP.inp' not in file_list:
                print name(Arr,Agr)
                for f in file_list:
                    if 'src' in f:
                        for sf in os.listdir(temp_path+f):
                            os.remove(temp_path+f+'/'+sf)
                        os.removedirs(temp_path+f)
                    else:
                        os.remove(temp_path+f)
                #os.removedirs(temp_path)
                print command(Arr,Agr)
                os.system(command(Arr,Agr)+' r')
        else:
            os.system(command(Arr,Agr))

