# -*- coding: utf-8 -*-
import tempfile
import os
import pylab
from input_param import Input
import pickle
import subprocess
import shutil

def submit_job(input_file, job_name ,Report_name , binary,path):
    """Submit a jot using qsub"""
    path = path + '/'
    template_file = 'FP.job_Template'
    try:
        f = open(template_file,"r")
        allLines = f.readlines()
    finally:
        f.close()
    # Move the binary file Report directory.
    binary_new = path+job_name
    shutil.move(binary,binary_new)
    # Copy source files
    os.mkdir(path+'src')
    src_files = os.listdir('./')
    for sf in src_files:
        if sf.split('.')[-1] in ['F90','f90','py']:
            shutil.copy(sf,path+'src/')
    tempf = tempfile.NamedTemporaryFile(suffix='.job',delete=False)
    if (job_name == ''):
        job_name = exp_label+'_'+str(param)
    for eachLine in allLines:
        eachLine = eachLine.replace('@job_name', job_name.replace("/","_"), 1)
        eachLine = eachLine.replace('@input_file', input_file, 1)
        eachLine = eachLine.replace('@output_path', path, 1)
        eachLine = eachLine.replace('@binary', binary_new, 1)
        tempf.write(eachLine)
    tempf.close()
    #print a	
    command = 'qsub ' + tempf.name 
    os.system(command)
#    subprocess.Popen(command)
    os.unlink(tempf.name)



def submit_job_nml(Report_name, Exp_name, path, **keywords):
    """ Run a MPI job on SGE """
  # Make a namelist file.
    path = path+Report_name
    print keywords['binary']
    input_parmas=Input(Report_name, Exp_name,path=path, **keywords)
    if (input_parmas.namelist_text):
	if not os.path.exists(path):
            os.makedirs (path)
        input_file_name = path+'/'+Exp_name+'.input'
        input_parmas.write_namelist(input_file_name)
        submit_job(input_file_name, Exp_name, Report_name, keywords['binary'],path)
    pass

def make(PPA=''):  # Preprocessor arguments
    p1 = subprocess.Popen(['make','clean'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.output = p1.stdout.readlines()
    p1.error  = p1.stderr.readlines() 
    for l in p1.output:
        print l
    if not p1.error:
        p2 = subprocess.Popen(['make','PPA=' + PPA], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2.output = p2.stdout.readlines()
        p2.error  = p2.stderr.readlines() 
        for l in p2.output:
            print l
        if not p2.error:
            return 1
        else:
            for l in p2.error: 
                print l  
    else:
        for l in p1.error: 
            print l 
            return 0 
    pass

