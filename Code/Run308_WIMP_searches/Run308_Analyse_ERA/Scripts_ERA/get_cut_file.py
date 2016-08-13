#!/usr/bin/env python

from ROOT import TTree, TFile, TH2F, TH1D, TGraph, TH1F
import script_utils as script_utils
import sys,os

def get_cut_file(bolo_name, data_dir, tree_name):
    """
    Creates a .txt file with the prelim cut list of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file '_prelim_TCuts.txt'
    it contains the prelim cuts fin TCut form for each bolometer and is updated when a new one is added
    """


    # Define path for the .txt file. Create the directory if it does not exist, then open file
    cut_path_name= script_utils.create_directory('../Cut_files/')  
    cut_file_name="standard_cut_values.txt" 
    file_cut="" 
    #Add an exception if the file does not exist
    try:
        file_cut = script_utils.open_text_file(cut_path_name, cut_file_name , "r")
    except IOError:
        script_utils.print_utility(script_utils.COL("No such file, use get_standard_cuts.py first","fail"))
        sys.exit()

    # Load the cut values. Start from line 1 to remove header.
    list_file_cut_lines=[line.rstrip().split(",") for line in file_cut.readlines()[1:]]
    bolo_TCut_line=bolo_name+","
    # Add a boolean flag to check if the bolo has its cuts in the file
    is_bolo_in_file=False
    for line in list_file_cut_lines:
        if bolo_name == line[0]:
            bolo_TCut_line +="FWIA<"+line[1]+"&&FWIB<"+line[2]+"&&FWIC<"+line[3]+"&&FWID<"+line[4]+"&&OWC1<"+line[5]+"&&OWC2<"+line[6]
            bolo_TCut_line +="&&(CHIA-RCIA)<"+line[7]+"&&(CHIB-RCIB)<"+line[8]+"&&(CHIC-RCIC)<"+line[9]+"&&(CHID-RCID)<"+line[10]+"&&XOC1<"+line[11]+"&&XOC2<"+line[12]
            bolo_TCut_line += "&&abs(EC1-EC2)<1" 
            bolo_TCut_line += "&&TBEF>0.5&&TAFT>0.5&&APAT>1&&MULT==1" 
            is_bolo_in_file=True

    #Add an error message if the bolo is no found
    if not is_bolo_in_file:
        script_utils.print_utility(script_utils.COL("Bolo no in the cut file. Verify process", "fail"))
        sys.exit()

    # Create file if it does no exist
    TCut_file_name='../Cut_files/TCuts.txt'

    #2 cases: the file exists or does not.
    if not os.path.isfile(TCut_file_name) :
        with open(TCut_file_name, 'w') as TCut_file: 
            TCut_file.write(bolo_TCut_line+ "&&KTH>0&&KTH<1\n")
            script_utils.print_utility(script_utils.COL('Creating file with line ' + bolo_TCut_line,'blue'))  

    else:
        # Read the whole file, save its contents.
        TCut_lines=[]
        with open(TCut_file_name, 'r') as TCut_file:
            TCut_lines=TCut_file.readlines()

        #Rewrite it with the appropriate content
        with open(TCut_file_name, 'w') as TCut_file:        
            # Update if the line already exists. Else, add it to the file
            is_bolo_in_line=False
            for line in TCut_lines:
                if bolo_name in line:
                    TCut_file.write(bolo_TCut_line+ "&&KTH>0&&KTH<1\n")
                    script_utils.print_utility(script_utils.COL('Updating line to ' + bolo_TCut_line,'blue'))
                    is_bolo_in_line=True
                else:
                    TCut_file.write(line)

            if not is_bolo_in_line:
                TCut_file.write(bolo_TCut_line+ "&&KTH>0&&KTH<1\n")
                script_utils.print_utility(script_utils.COL('Adding line ' + bolo_TCut_line,'blue'))

  
    
    # Create file if it does no exist
    TCut_file_name='../Cut_files/TCuts_noKTHcut.txt'


    #2 cases: the file exists or does not.
    if not os.path.isfile(TCut_file_name) :
        with open(TCut_file_name, 'w') as TCut_file: 
            TCut_file.write(bolo_TCut_line+'\n')
            script_utils.print_utility(script_utils.COL('Creating file with line ' + bolo_TCut_line,'blue'))  

    else:
        # Read the whole file, save its contents.
        TCut_lines=[]
        with open(TCut_file_name, 'r') as TCut_file:
            TCut_lines=TCut_file.readlines()

        #Rewrite it with the appropriate content
        with open(TCut_file_name, 'w') as TCut_file:        
            # Update if the line already exists. Else, add it to the file
            is_bolo_in_line=False
            for line in TCut_lines:
                if bolo_name in line:
                    TCut_file.write(bolo_TCut_line+'\n')
                    script_utils.print_utility(script_utils.COL('Updating line to ' + bolo_TCut_line,'blue'))
                    is_bolo_in_line=True
                else:
                    TCut_file.write(line)

            if not is_bolo_in_line:
                TCut_file.write(bolo_TCut_line+ '\n')
                script_utils.print_utility(script_utils.COL('Adding line ' + bolo_TCut_line,'blue'))

    # Create file if it does no exist
    TCut_file_name='../Cut_files/TCuts_forduration.txt'
    bolo_TCut_line=bolo_name+","
    for line in list_file_cut_lines:
        if bolo_name == line[0]:
            bolo_TCut_line +="FWIA<"+line[1]+"&&FWIB<"+line[2]+"&&FWIC<"+line[3]+"&&FWID<"+line[4]+"&&OWC1<"+line[5]+"&&OWC2<"+line[6]
    bolo_TCut_line+="&&KTH>0&&KTH<1&&APAT>1&&MULT==1"

    #2 cases: the file exists or does not.
    if not os.path.isfile(TCut_file_name) :
        with open(TCut_file_name, 'w') as TCut_file: 
            TCut_file.write(bolo_TCut_line+'\n')
            script_utils.print_utility(script_utils.COL('Creating file with line ' + bolo_TCut_line,'blue'))  

    else:
        # Read the whole file, save its contents.
        TCut_lines=[]
        with open(TCut_file_name, 'r') as TCut_file:
            TCut_lines=TCut_file.readlines()

        #Rewrite it with the appropriate content
        with open(TCut_file_name, 'w') as TCut_file:        
            # Update if the line already exists. Else, add it to the file
            is_bolo_in_line=False
            for line in TCut_lines:
                if bolo_name in line:
                    TCut_file.write(bolo_TCut_line+'\n')
                    script_utils.print_utility(script_utils.COL('Updating line to ' + bolo_TCut_line,'blue'))
                    is_bolo_in_line=True
                else:
                    TCut_file.write(line)

            if not is_bolo_in_line:
                TCut_file.write(bolo_TCut_line+ '\n')
                script_utils.print_utility(script_utils.COL('Adding line ' + bolo_TCut_line,'blue'))