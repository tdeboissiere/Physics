#!/usr/bin/env python

from ROOT import TFile
import os, sys

"""
Collection of utilities used in the analysis scripts
"""



def COL(string_name, color):
	"""
    Attaches color prefixes/suffixes to given string

    Detail:

    Arguments:
    string_name (str) the string name
    color       (str) the desired color

    Outputs:

    returns a string enhanced with color  pref/suff
    """

	if color =='blue':
		pref = '\033[94m'
        
	if color =='header':
		pref = '\033[95m'
        
	if color =='green':
		pref = '\033[92m'

	if color =='warning':
		pref = '\033[93m'
        
	if color =='fail':
		pref = '\033[91m'
                
	suff = '\033[0m'
	return pref + string_name + suff

def print_utility(string_name):
	"""
	prints a string, highlighted with *****
	
	Detail:
	
	Arguments:
	string_name (str) the string name
	
	Outputs:
	
	Highlights a string print output with ****
	"""	

	print ''
	print '********* '+ string_name +' ************'
	print '' 

def create_directory(path_name):
	"""
    Creates a directory with name path_name

    Detail:

    Arguments:
    path_name (str) the full directory path name

    Outputs:

    Creates a directory at the given location
    """
	if not os.path.exists(path_name): 
		os.makedirs(path_name)
		print_utility(COL("Creating directory " +path_name, 'blue'))
		return path_name
	return path_name

def open_text_file(path_name, file_name, option):
	"""
    Opens a file specified by path + name for writing or reading

    Detail:

    Arguments:
    path_name (str) the  path name
    file_name (str) the file name
    option    (str)  "r" or "w" for reading or writing resp.

    Outputs:

    Creates the file path_name+file_name for reading or writing
    """

	newfile = open(path_name + file_name, option)
	if option =='r':
		print_utility(COL('Opening ' + path_name+file_name + ' for reading', 'blue'))
		return newfile
	if option =='w':
		print_utility(COL('Opening ' + path_name+file_name + ' for writing', 'blue'))
		return newfile
	else:
		print_utility(COL('Option ' + option + ' not valid. Exiting ', 'fail'))
		sys.exit()

def open_ROOT_file(ROOT_path_name, ROOT_file_name, option):
	"""
    Opens a ROOT file specified by path + name 

    Detail:

    Arguments:
    path_name (str) the  path name
    file_name (str) the file name
    option    (str)  "read" or "recreate" for reading or writing resp.

    Outputs:

    Creates the  ROOT file path_name+file_name
    """

	newfile = TFile(ROOT_path_name + ROOT_file_name, option)
	if option =='read':
		print_utility(COL('Opening ' + ROOT_path_name+ROOT_file_name + ' for reading', 'blue'))
		return newfile
	if option =='recreate':
		print_utility(COL('Opening ' + ROOT_path_name+ROOT_file_name + ' for writing', 'blue'))
		return newfile
	else:
		print_utility(COL('Option ' + option + ' not valid. Exiting ', 'fail'))
		sys.exit()

def get_howmany(bolo_name) :
	"""
	Returns the number of polarisations for a given bolo_name
	
	Detail:
	
	Arguments:
	bolo_name (str) the bolometer name
	
	Outputs:

	The number of polarisations
	
	"""
	#Define path for the file. Create the directory if it does not exist
	path_name  = create_directory('../Analyse_' + bolo_name + '/Text_files/')

	file_polar = open_text_file(path_name, bolo_name + "_all_polars_with_entries.txt", "r")
	list_lines = file_polar.readlines()
	#Ignore first line hence the -1
	print_utility('Found ' + str(len(list_lines) -1) + ' polar condition(s)')
	return len(list_lines) -1