# -*- coding: utf-8 -*-
import pylab

class Input(object):

  def __init__(self, Report_name, Exp_name, **keywords):
      self.output_dir = Report_name
      self.out_name = Exp_name
      self.input_params = keywords
      self.error = False
      self.set_param()
      self.get_namelist()
      if not self.error:
	  self.make_namelist()
      else:
	  self.namelist_text = False
  def get_namelist(self):
      try:
	namelist_template = open('FP.inp_template',"r")
      except:
	print "Error: 'FP.inp_template' dosn't exist"
	self.error = True
      else:
	self.namelist_text = namelist_template.readlines()
	namelist_template.close()
  
  def set_param(self):
      self.params = dict.fromkeys([
        'path',
        'ARR',
	'ANR',
        'AGR',
        'Mass_spectrum',
        'Single_mass',
        'NPDE',
        'Aeps',])
      
      self.defaults = {
        'ANR'          : 1.0,
        'ARR'          : 1.0,
        'AGR'          : 1.0,
        'Aeps'         : 1.0,
        'Mass_spectrum': False,
        'Single_mass'  : True,
        'NPDE'         : 1,
        }
      
      self.params.update(self.defaults, **self.input_params)
      if self.params['Mass_spectrum']:
        self.params['Mass_spectrum'] = ''
        self.params['Single_mass'] = '!'
        self.params['NPDE'] = 4
      else:
        self.params['Mass_spectrum'] = '!'
        self.params['Single_mass'] = ''
        self.params['NPDE'] = 1
      if None in self.params.values():
	  print "Error '" + self.params.keys()[self.params.values().index(None)] + "' has no value"
	  self.error = True
  
  def make_namelist(self):
      for param in self.params.keys():
	  for i in xrange(len(self.namelist_text)):
	      self.namelist_text[i] = self.namelist_text[i].replace('@'+param,str(self.params[param]))
  
  def write_namelist(self,file_name):
      try:
	namelist_file = open(file_name,"w")
      except:
	print "Error: "+file_name+" dosn't exist"
	self.error = True
      else:
        namelist_file.writelines(self.namelist_text)
	namelist_file.close()
