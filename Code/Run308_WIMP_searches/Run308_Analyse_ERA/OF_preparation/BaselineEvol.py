#! /usr/bin/env python

import sys,os
import numpy as np
from ROOT import *

def get_baseline_evol(bolo_name,channel, interactif=0):

    """Add ERA heat to Jules NTP
    
    Detail:
        void

    Args:
        bolo_name    = (str) bolo name 
        channel_name = (str) channel name (Chal1 or Chal2)
        interactif   = (0/1) Set to 1 to see the plot
        
    Returns:
        void

    Raises:
        void
    """

    #Get the calibration coefficients
    calib_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff.txt"
    d_calib = {}
    assert(os.path.isfile(calib_file))
    with open(calib_file, "r") as fc:
        lines = fc.readlines()
        calibcoeff1, calibcoeff2 = lines[1].rstrip().split(",")[1], lines[1].rstrip().split(",")[2]  
        d_calib["Chal1"]=calibcoeff1
        d_calib["Chal2"]=calibcoeff2

    calibcoeff = float(d_calib[channel])

    list_periods = []
    with open("./Text_files/" + bolo_name + "/liste_periods.txt", "r") as fperiod:
        stuff = fperiod.readlines()
        list_periods = [elem.rstrip().split(" ") for elem in stuff[1:]]

    # print list_periods

    nbins=100 
    min_histo=-2
    max_histo=2

    class Gauss:
       def __call__( self, x, par ):
          return par[0]*TMath.Exp(-0.5*((x[0]-par[1])/par[2])**2)

    normal_law = TF1("f", Gauss(), -2,2,3)
    normal_law2 = TF1("f2", Gauss(), -2,2,3)
    normal_law.SetParameters(1,0,0.5)
    normal_law2.SetParameters(1,0,0.5)
    normal_law.FixParameter(1,0)
    normal_law.SetParLimits(2,0.01, 5)
    normal_law2.FixParameter(1,0)
    normal_law2.SetParLimits(2,0.01, 5)


    with open("./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_" + channel + ".txt", "w") as outputfile:
        for index, theperiod in enumerate(list_periods) :
            tinf=int(theperiod[0])
            tsup=int(theperiod[1])
            run=theperiod[2]
            
            fbasic=TFile("../Amp_files/" + bolo_name + "/basicntp_"+run+"_"+bolo_name+".root","READ")
            tbasic=fbasic.Get("basicntp_"+bolo_name)
            f=TFile("../Amp_files/" + bolo_name + "/tmpwiener_"+run+"_"+channel+"_"+bolo_name+".root","READ")
            t=f.Get("wienerntp_"+bolo_name+"_"+channel)

            Baseline=TH1F("Baseline","Baseline",nbins,min_histo,max_histo)

            for ievt in range(tbasic.GetEntries()) :
                tbasic.GetEntry(ievt)
                if tbasic.DateSec<tinf or tbasic.DateSec>tsup: continue
                t.GetEntry(ievt)
                cuts_part1 = tbasic.IsBoloTrigger == 0 and tbasic.Saturation ==0 and t.WienerZeroAmpl!=0
                cuts_part2 = t.WienerChi2<1.5 and t.WienerChi2>0.5 and t.WienerZeroAmpl <10 and eval("tbasic.PileUp"+channel+"==0")
                if ( cuts_part1 and cuts_part2) :
                    zeroampl_kev =t.WienerZeroAmpl/calibcoeff
                    Baseline.Fill(zeroampl_kev)

            # Calcul du fwhm a partir de l'histo Baseline:
            fwhm = 0
            if Baseline.GetEntries()==0 : 
        	   print run,tinf, tsup," : No data to compute baseline.."
        	   fwhm=0
            else :
                if interactif==1 :

                    Baseline.Fit("f","NQ", "", -2, 0)
                    Baseline.Fit("f2","NQ")

                    cc = TCanvas("cc", "cc")
                    Baseline.Draw()
                    normal_law.Draw("same")
                    normal_law.SetLineColor(kBlue)
                    normal_law2.Draw("same")
                    fwhm = 2.3548*normal_law.GetParameter(2)
                    fwhm2 = 2.3548*normal_law2.GetParameter(2)
                    gPad.Update()
                    print "Blue", fwhm, "Red", fwhm2

                    raw_input("Next..")
                    # gPad.SaveAs("./Figures/baseline_fit/"+channel+"/"+str(index)+".png")
                    del cc
                else :
                    Baseline.Fit("f","QN", "", -2, 0)
                    fwhm=2.35482*normal_law.GetParameter(2)

            # Le FWHM a partir de la liste des SigmaAmpl:
            outputfile.write(str(tinf)+","+str(tsup)+","+str(fwhm)+"\n")
            f.Close()
            fbasic.Close()
            del tbasic
            del t


def refine_baseline(bolo_name):

    """Recalibrate baseline given the refined calib coefficients
    
    Detail:
        void

    Args:
        bolo_name    = (str) bolo name 
        
    Returns:
        void

    Raises:
        void
    """

    #Get the old calibration coefficients
    calib_file = "./Text_files/" + bolo_name + "/" + bolo_name + "_calibcoeff_refined.txt"
    calibcoeff1, calibcoeff2=1,1
    assert(os.path.isfile(calib_file))
    with open(calib_file, "r") as fc:
        lines = fc.readlines()
        calibcorrection1, calibcorrection2 = float(lines[1].rstrip().split(",")[3]), float(lines[1].rstrip().split(",")[4]) 

    data_types = {"names": ("tinf", "tsup", "FWHM"), "formats": ("i", "i", "f")}
    arr_FWHM_EC1 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal1.txt", delimiter=",",  dtype=data_types)
    arr_FWHM_EC1["FWHM"] = arr_FWHM_EC1["FWHM"]*calibcorrection1

    data_types = {"names": ("tinf", "tsup", "FWHM"), "formats": ("i", "i", "f")}
    arr_FWHM_EC2 = np.loadtxt("./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal2.txt", delimiter=",",  dtype=data_types)
    arr_FWHM_EC2["FWHM"] = arr_FWHM_EC2["FWHM"]*calibcorrection2

    np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal1_refined.txt", arr_FWHM_EC1, delimiter = ",", fmt="%i,%i,%f")    
    np.savetxt("./Text_files/" + bolo_name + "/" + bolo_name + "_baseline_evol_Chal2_refined.txt", arr_FWHM_EC2, delimiter = ",", fmt="%i,%i,%f")