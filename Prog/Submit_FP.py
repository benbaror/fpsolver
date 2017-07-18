# -*- coding: utf-8 -*-
import pylab as pl
import subprocess
import jobs as jobs
import os
import time
from optparse import OptionParser
def main():
	path = '/home/ben/mp/Data/FP_Solver/Runs/'
	Path = '/home/ben/mp/Data/FP_Solver/Runs/'
	usage = "usage: %prog report_name"
	parser = OptionParser(usage)
	parser.add_option("-r", "--redo", action="store_true", dest="redo",
                      help="Rerun to redo rejected jobs")
        parser.add_option("-S", "--Single_mass", action="store_true", 
                      dest="Single_mass", help="Use Single Mass")
	parser.add_option("--AGR",action="store", type="float", dest="AGR", default=1.0,
			  help='AGR')
	parser.add_option("--ARR",action="store", type="float", dest="ARR", default=1.0,
			  help='ARR')
        parser.add_option("--ANR",action="store", type="float", dest="ANR", default=1.0,
			  help='ANR')
	parser.add_option("-E", "--experiment",
			  action="store", type="string",
			  dest="Experiment", default='SM_NR_AH09',
			  help="Experiment ('SM_NR_HA06a','SM_RR_HA06a','SM_NR_AH09','SM_RR_AH09','SM_NR_AH09_E09','SM_RR_AH09_E09')")
	(options, args) = parser.parse_args()
	if len(args) != 1:
		parser.error("incorrect number of arguments")
	Report_name = args[0]
	if os.path.exists(path+Report_name):
		if not options.redo:
			parser.error("Error: " + Report_name + " already exsist")
	elif  options.redo:
		parser.error("Error: can't redo " + Report_name + " Repot dosen't exsist")
	else:
		os.makedirs (path+Report_name)
	print Report_name
	Experiment={
		#                  [Aeps, ANR, ARR, Flags,Mass spectrum]
		'SM_NR_HA06a'   :  [0.24,0.24,3.56,'',False],
		'SM_RR_HA06a'   :  [0.24,0.24,3.56,'-D do_RR',False],
		'SM_NR_AH09'    :  [0.24,0.24,3.56,'-D pinhole',False],
		'SM_RR_AH09'    :  [0.24,0.24,3.56,'-D do_RR -D pinhole',False],
		'SM_NR_AH09_E09':  [0.24,1/1.01**2,1/1.05**2,' -D pinhole',False],
		'SM_RR_AH09_E09':  [0.24,1/1.01**2,1/1.05**2,'-D do_RR -D pinhole',False],
		'MS_NR_AH09'    :  [0.24,0.24,3.56,'-D pinhole',True],
		'MS_RR_AH09'    :  [0.24,0.24,3.56,'-D do_RR -D pinhole',True],
		'MS_NR_AH09_E09'    :  [0.24,1/1.01**2,1/1.05**2,'-D pinhole',True],
		'MS_RR_AH09_E09'    :  [0.24,1/1.01**2,1/1.05**2,'-D do_RR -D pinhole',True],

		}
	job_name = options.Experiment
	Aeps = Experiment[options.Experiment][0]
	ANR = Experiment[options.Experiment][1]*options.ANR
	ARR = Experiment[options.Experiment][2]*options.ARR
	AGR = options.AGR
	PPA = ' '
	PPA = PPA + Experiment[options.Experiment][3]
	Mass_spectrum = Experiment[options.Experiment][4]
	if jobs.make(PPA = PPA):
		jobs.submit_job_nml(Report_name,job_name,Path,
				    binary = 'FPSolverMS',
				    ARR = ARR,
				    AGR = AGR,
				    ANR = ANR,
				    Mass_spectrum = Mass_spectrum,
				    Aeps = Aeps)
if __name__ == '__main__':
	main()
