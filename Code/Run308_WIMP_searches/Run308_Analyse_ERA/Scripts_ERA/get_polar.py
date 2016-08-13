#!/usr/bin/env python

from ROOT import TTree, TCut, TFile
from collections import namedtuple

import script_utils as script_utils


def get_polar(bolo_name, data_dir, tree_name ):
    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file bolo_name + "_all_polars_with_entries.txt", "w"
    it contains the polarisation name, the values of VVET and VOLT and the #entries

    """

    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)
    nEntries    = tree.GetEntries()
    
    Polars      = namedtuple('Polars', 'VVET VOLT')
    list_polars = []

    #Loop over tree to find all the (VVET VOLT) tuples
    for i in range(nEntries):
        tree.GetEntry(i)
        if Polars(tree.VVET, tree.VOLT) not in list_polars:
            list_polars.append(Polars(tree.VVET, tree.VOLT))

    #Define path for the file. Create the directory if it does not exist
    path_name=script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')

    #Write the results to a .txt file
    outputfile = script_utils.open_text_file(path_name, bolo_name + "_all_polars_with_entries.txt", "w")
    outputfile.write("POLAR_NAME" + "," + "POLAR_CUT_SELECTION" +"," +"NUMBER_ENTRIES" +","+"VFID"+","+"VET\n")
    for polar in list_polars:
        polar_cut     ="VVET=="+str(polar.VVET)+"&&VOLT=="+str(polar.VOLT)
        polar_entries =str(tree.GetEntries(polar_cut))
        polar_line=""
        #The following distrinction comes from a Jules sign convention
        if (polar.VVET>0 and polar_entries>0):
            polar_line=bolo_name + "+" +str(abs(polar.VOLT)/2.) + "V-" + str(polar.VVET) + "V" + "," + polar_cut+","+polar_entries+","+str(abs(polar.VOLT))+","+str(abs(polar.VOLT)/2.+polar.VVET)
            outputfile.write(polar_line + "\n")
        elif (polar.VVET <=0 and polar_entries>0):
            polar_line = bolo_name+ "+" + str(abs(polar.VOLT)/2.)+"V+"+ str(-polar.VVET) + "V" +","+polar_cut+","+polar_entries+","+str(abs(polar.VOLT))+","+str(abs(polar.VOLT)/2.+polar.VVET)
            outputfile.write(polar_line + "\n")
        
        #Print file content
        script_utils.print_utility(polar_line)

    outputfile.close()
