#!/usr/bin/env python

from ROOT import TCanvas, TFile, TTree
import script_utils as script_utils
import os

def get_standard_cuts(bolo_name, data_dir, tree_name):

    """
    Creates a .txt file with the polar information of a given bolometer

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the tree directory (containing Jules n-tuple)
    tree_name (str) the tree name in the tree file

    Outputs:

    A .txt file "standard_cut_values".txt
    it contains the cut values on all the FWHM and CHI2 and is updated when a new one is added
    """
    
    #Get tree 
    file_tree   = TFile(data_dir+bolo_name+"_fond.root")
    tree        = file_tree.Get(tree_name)
    nEntries    = tree.GetEntries()

    #Initialize the start and end time 
    tree.GetEntry(0)
    all_tree_start_time =str(1E6*tree.UT1+tree.UT2)
    all_tree_end_time   =str(1E6*tree.UT1+tree.UT2)
    
    #Find the start and end time of the whole data tree
    for i in range(1, nEntries):
        tree.GetEntry(i)
        time= 1E6*tree.UT1+tree.UT2
        if time < all_tree_start_time:
            all_tree_start_time = str(time)
        if time > all_tree_end_time:
            all_tree_end_time = str(time)
    


    c_fwhm_ion=TCanvas("c_fwhm_ion", "c_fwhm_ion")
    c_fwhm_ion.Divide(2,2)

    c_fwhm_ion.cd(1)
    tree.Draw("FWIA:1E6*UT1+UT2>>fwhist1(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")
    c_fwhm_ion.cd(2)
    tree.Draw("FWIB:1E6*UT1+UT2>>fwhist2(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")
    c_fwhm_ion.cd(3)
    tree.Draw("FWIC:1E6*UT1+UT2>>fwhist3(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")
    c_fwhm_ion.cd(4)
    tree.Draw("FWID:1E6*UT1+UT2>>fwhist4(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")

    c_fwhm_heat=TCanvas("c_fwhm_heat", "c_fwhm_heat")
    c_fwhm_heat.Divide(1,2)

    c_fwhm_heat.cd(1)
    tree.Draw("OWC1:1E6*UT1+UT2>>fwist5(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")
    c_fwhm_heat.cd(2)
    tree.Draw("OWC2:1E6*UT1+UT2>>fwhist6(10000," + all_tree_start_time + "," + all_tree_end_time + ",1000,-2,10)")
    
    # Manually fill in the CHI2 cut after looking at the data.
    str_cut_fwhm = raw_input('Enter (FWIA,FWIB, FWIC, FWID, OWC1, OWC2) cut (a, b, c, d, e, f) format: ')

    fwhm_cut_line = str_cut_fwhm.split(",")
    bolo_TCut_line ="FWIA<"+fwhm_cut_line[0]+"&&FWIB<"+fwhm_cut_line[1]+"&&FWIC<"+fwhm_cut_line[2]+"&&FWID<"+fwhm_cut_line[3]+"&&OWC1<"+fwhm_cut_line[4]+"&&OWC2<"+fwhm_cut_line[5]
    bolo_TCut_line += "&&FWIA>0&&FWIB>0&&FWIC>0&&FWID>0&&OWC1>0&&OWC2>0 && APAT>1"

    print bolo_TCut_line

    c_chi2_heat_of_energ=TCanvas("c_chi2_heat_of_energ", "c_chi2_heat_of_energ")
    c_chi2_heat_of_energ.Divide(1,2)

    c_chi2_heat_of_energ.cd(1)
    tree.Draw("pow(10,XOC1):0.5*(EC1+EC2)>>chist5(10000,0,20,100,-2,10)", bolo_TCut_line)
    c_chi2_heat_of_energ.cd(2)
    tree.Draw("pow(10,XOC2):0.5*(EC1+EC2)>>chist6(10000,0,20,100,-2,10)", bolo_TCut_line)


    str_cut_chi2 = raw_input('Enter (CHIA-RCIA, CHIB-RCIB, CHIC-RCIC, CHID-RCID, XOC1, XOC2) cut (a, b, c, d, e, f) format: ')

    # Create file if it does no exist
    cut_file_name='../Cut_files/standard_cut_values.txt'

    #2 cases: the file exists or does not.
    if not os.path.isfile(cut_file_name) :
        with open(cut_file_name, 'w') as cut_file: 
            #Include header
            cut_file.write("Bolo_name,FWIA,FWIB,FWIC,FWID,OWC1,OWC2,CHIA,CHIB,CHIC,CHID,XOC1,XOC2" + "\n")
            bolo_cut_line =bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
            cut_file.write(bolo_cut_line+ '\n')
            script_utils.print_utility(script_utils.COL('Creating file with line ' + bolo_cut_line,'blue'))       

    else:
        # Read the whole file, save its contents.
        cut_lines=[]
        with open(cut_file_name, 'r') as cut_file:
            cut_lines=cut_file.readlines()

        #Rewrite it with the appropriate content
        with open(cut_file_name, 'w') as cut_file:
            #Include header
            cut_file.write("Bolo_name,FWIA,FWIB,FWIC,FWID,OWC1,OWC2,CHIA,CHIB,CHIC,CHID,XOC1,XOC2" + "\n")

            # Update if the line already exists. Else, add it to the file
            is_bolo_in_line=False
            for line in cut_lines[1:]:
                if bolo_name in line:
                    bolo_cut_line=bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
                    cut_file.write(bolo_cut_line+'\n')
                    script_utils.print_utility(script_utils.COL('Updating line to ' + bolo_cut_line,'blue'))
                    is_bolo_in_line=True
                else:
                    cut_file.write(line)

            if not is_bolo_in_line:
                bolo_cut_line=bolo_name+','+ str_cut_fwhm + ',' + str_cut_chi2
                cut_file.write(bolo_cut_line+ '\n')
                script_utils.print_utility(script_utils.COL('Adding line ' + bolo_cut_line,'blue'))












