#! /usr/bin/env python

from ROOT import *
import sys
import script_utils as script_utils
import create_extrapolated_bckg as create_bckg
import subprocess as subp
import BDT_utils as BDT_ut
import BDT_file_handler as BDT_fh

def launch_analysis(bolo_name, analysis_type, FWHM_type, d_cut, d_overwrite, nsimu, use_ERA):

	"""Launch BDT simulations
	
	Detail:
		Detailed description

	Args:
		bolo_name     = (str) bolometer name
		analysis_type = (str) type of analysis (cuts on heat and ion)
		FWHM_type     = (str) which FWHM settings
		d_cut         = (dict) indicates the cuts (inf/sup) on heat and ion 
		nsimu         = (int) number of true event simulations
		overwrite     = (dict) indicates if standard files need to be overwritten
		use_ERA       = (0/1) indicates if use ERA heat only spectrum
		
	Returns:
		void

	Raises:
		AssertionError 
	"""



bolo_name     = "FID837"
analysis_type = "standard_resolution_no_ion_cut"
FWHM_type     = "standard_resolution"
d_cut         = {"ECinf": 1, "ECsup": 15, "EIinf": 0, "EIsup": 15}
d_overwrite   = {"Gamma": 1, "BetaPb": 1, "True":1, "Heatonly":1, "WIMP":1, "Extra":0}
nsimu         = 10
use_ERA       = 1

launch_analysis(bolo_name, analysis_type, FWHM_type, d_cut, d_overwrite, nsimu, use_ERA)