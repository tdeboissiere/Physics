#! /usr/bin/env python

################################################################################
### Routines to read Jules files mostly, for Run 308 only
################################################################################

import pickle,os,sys,glob,string,random
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d
from ROOT import *
from utils_python import *
import PyROOTPlots as PyRPl
from array import array

def compute_coeff(FWHM1, FWHM2):
    # Compute the coeff of the best estimator of 2 quantities with resolutions FWHM1 and FWHM2
    s1=FWHM1/2.3548
    s2=FWHM2/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)
    return w1

################################################################################

def compute_resolution(FWHM1, FWHM2) :
    # Compute the resolution of the best estimator of 2 quantities with resolutions FWHM1 and FWHM2
    s1=FWHM1/2.3548
    s2=FWHM2/2.3548

    w1=(s2*s2)/(s1*s1+s2*s2)

    res=np.sqrt(np.power(w1*s1,2)+np.power((1-w1)*s2,2))
    return 2.3548*res

################################################################################

def get_electronic_cut_line(bolo_name) :
    # Open a .txt file where a cut is stored to get electronic recoil events
    with open("../Text_files/elec_recoil_cuts.txt", "r") as f:
        list_lines = f.readlines()
        d_cuts = {}
        for cut_line in list_lines:
            line = cut_line.rstrip().split(",")
            d_cuts[line[0]] = line[1]

        try :
            return d_cuts[bolo_name]
        except KeyError :
            print "Bolo not in cut file"
            


################################################################################

def getbolonum_ana(bolo) :
    # liste numeros de bolos ana-run308 en dur
    theliste_nums = {"FID823":5,"FID824":6,"FID825":7,"FID810":8,"FID826":9,"FID827":10,"FID828":11,"FID817":12,"FID837":21,"FID838":22,"FID839":23,"FID831":24,"FID820":25,"FID821":26,"FID822":27,"FID807":28,"FID840":29,"FID841":30,"FID842":31,"FID832":32,"FID843":33,"FID844":34,"FID845":35,"FID846":36}
    return theliste_nums[bolo]

################################################################################

def getmac(bolo) :
    # liste numeros des macs par bolo en dur
    theliste_macs = {"FID823":1,"FID824":1,"FID825":1,"FID810":1,"FID826":1,"FID827":1,"FID828":1,"FID817":1,"FID837":2,"FID838":2,"FID839":2,"FID831":2,"FID820":3,"FID821":3,"FID822":3,"FID807":3,"FID840":3,"FID841":3,"FID842":3,"FID832":3,"FID843":3,"FID844":3,"FID845":3,"FID846":3}
    return theliste_macs[bolo]

################################################################################

def getbittrig_ana(bolo) : 
    # liste "bit trigger local" = variable TR dans ana
    # seulement pour les 9 bolos basse masse. Determination faite dans sandbox
    theliste_tr={"FID824":1,"FID825":2,"FID826":4,"FID827":5,"FID837":0,"FID838":1,"FID839":2,"FID841":5,"FID842":6}
    return theliste_tr[bolo]

################################################################################

def skim_periodlist(periodlist,listcuts=["kth","kthlow","fwf","fwc","fwvet"]) :
    w_outlier = ((periodlist['FWF']>0.1)&(periodlist['FWIA']>0.1)&(periodlist['FWIB']>0.1)&(periodlist['FWIC']>0.1)&(periodlist['FWID']>0.1)&(periodlist['FWC']>0.1)&(periodlist['Kth']>0.1)&(periodlist['Kth']<20.)&(periodlist['FWF']<9.9)&(periodlist['FWC']<9.9)&(periodlist['FWIA']<9.9)&(periodlist['FWIB']<9.9)&(periodlist['FWIC']<9.9)&(periodlist['FWID']<9.9)&(periodlist['Sat']==0))
    w_kth=(periodlist['Kth']<1.5)
    w_kthlow=(periodlist['Kth']>0.1+0.6*periodlist['TWC'])
    w_fwf=(periodlist['FWF']<0.7)
    w_fwc=(periodlist['FWC']<1.0)
    w_fwvet=((periodlist['FWIA']<1.5)&(periodlist['FWIC']<1.5))
    w_all=w_outlier
    if "kth" in listcuts : w_all=(w_all&w_kth)
    if "kthlow" in listcuts : w_all=(w_all&w_kthlow)
    if "fwf" in listcuts : w_all=(w_all&w_fwf)
    if "fwc" in listcuts : w_all=(w_all&w_fwc)
    if "fwvet" in listcuts : w_all=(w_all&w_fwvet)
    periods_cut=periodlist[w_all]
    return periods_cut

################################################################################

def get_chiion_cut(bolo,return_string=True,cut_only_positive=True) :
    # liste des parametres de cuts chi2 ion. En dur, seulement pour les 9 bolos basse masse
    # Cuts identiques a la note de Jules maj fin juillet 2015 (similaire Christoph) 
    # SAUF pour 827 (optimisation du cut faite dans boloperf_measurement.py)

    p_chii=[99.,99.,99.,99.]
    if bolo=="FID824" : p_chii=[0.34,0.50,0.44,0.49]
    if bolo=="FID825" : p_chii=[0.24,0.43,0.30,0.27]
    if bolo=="FID826" : p_chii=[0.25,0.26,0.24,0.27]
    if bolo=="FID827" : p_chii=[0.39,0.34,0.47,0.24] # CHANGE PAR RAPPORT A JULES
    # Valeur jules [0.58,0.48,0.69,0.24]
    if bolo=="FID837" : p_chii=[0.24,0.19,0.16,0.21]
    if bolo=="FID838" : p_chii=[0.19,0.20,0.17,0.22]
    if bolo=="FID839" : p_chii=[0.16,0.15,0.19,0.14]
    if bolo=="FID841" : p_chii=[0.21,0.22,0.25,0.22]
    if bolo=="FID842" : p_chii=[0.14,0.21,0.20,0.25]

    if not return_string : return p_chii
    if cut_only_positive : cut_chi_ion="( (CHIA-RCIA)<"+str(p_chii[0])+" && (CHIB-RCIB)<"+str(p_chii[1])+" && (CHIC-RCIC)<"+str(p_chii[2])+" && (CHID-RCID)<"+str(p_chii[3])+" )"
    else : cut_chi_ion="( (CHIA-RCIA)>-"+str(p_chii[0])+" && (CHIA-RCIA)<"+str(p_chii[0])+" && (CHIB-RCIB)>-"+str(p_chii[1])+" && (CHIB-RCIB)<"+str(p_chii[1])+" && (CHIC-RCIC)>-"+str(p_chii[2])+" && (CHIC-RCIC)<"+str(p_chii[2])+" && (CHID-RCID)>-"+str(p_chii[3])+" && (CHID-RCID)<"+str(p_chii[3])+" )"
    return cut_chi_ion

################################################################################

def get_chichal_cut(bolo,return_string=True,chi2type="all") :
    # liste des parametres de cuts chi2 chal. En dur, seulement pour les 9 bolos basse masse
    # Cuts identiques a la note de Jules maj fin juillet 2015.
    if chi2type not in ["freq","time","all"] : print "Bad chi2type!!"
    
    p_chic_f=[99.,99.,99.,99.] # chi2 cut en frequence
    if bolo=="FID824" : p_chic_f=[0.18,0.24]
    if bolo=="FID825" : p_chic_f=[0.09,0.09]
    if bolo=="FID826" : p_chic_f=[0.24,0.21]
    if bolo=="FID827" : p_chic_f=[0.11,0.10]
    if bolo=="FID837" : p_chic_f=[0.10,0.12]
    if bolo=="FID838" : p_chic_f=[0.09,0.13]
    if bolo=="FID839" : p_chic_f=[0.09,0.09]
    if bolo=="FID841" : p_chic_f=[0.10,0.22]
    if bolo=="FID842" : p_chic_f=[0.09,0.10]
    
    p_chic_t=[99.,99.,99.,99.] # chi2 cut en temps
    if bolo=="FID824" : p_chic_t=[0.31,0.28]
    if bolo=="FID825" : p_chic_t=[0.30,0.16]
    if bolo=="FID826" : p_chic_t=[0.26,0.33]
    if bolo=="FID827" : p_chic_t=[0.95,0.32]
    if bolo=="FID837" : p_chic_t=[0.30,0.30]
    if bolo=="FID838" : p_chic_t=[0.18,0.27]
    if bolo=="FID839" : p_chic_t=[0.16,0.17]
    if bolo=="FID841" : p_chic_t=[0.21,0.27]
    if bolo=="FID842" : p_chic_t=[0.29,0.34]
    
    if not return_string :
        if chi2type=="freq" : return p_chic_f
        if chi2type=="time" : return p_chic_t
        if chi2type=="all" : return [p_chic_f,p_chic_t]

    cut_chi_chal_f="( XOC1>-"+str(p_chic_f[0])+" && XOC1<"+str(p_chic_f[0])+" && XOC2>-"+str(p_chic_f[1])+" && XOC2<"+str(p_chic_f[1])+" )"
    cut_chi_chal_t="( (CHIC1-RCC1)>-"+str(p_chic_t[0])+" && (CHIC1-RCC1)<"+str(p_chic_t[0])+" && (CHIC2-RCC2)>-"+str(p_chic_t[1])+" && (CHIC2-RCC2)<"+str(p_chic_t[1])+" )"
    if chi2type=="freq" : cut_chi_chal=cut_chi_chal_f
    if chi2type=="time" : cut_chi_chal=cut_chi_chal_t
    if chi2type=="all" : cut_chi_chal="("+cut_chi_chal_f+" && "+cut_chi_chal_t+")"
    return cut_chi_chal

################################################################################

def get_heatbiplot_cut(bolo,return_string=True) :
    # Liste cuts dans le plan (EC1,EC2), faits si OWC2 est assez bonne
    # Optimisation faite dans dataselection.py
    p_cut=[99999.,99999.]
    if bolo=="FID824" : p_cut=[1.0,1.3]
    if bolo=="FID825" : p_cut=[1.8,1.5]
    if bolo=="FID826" : p_cut=[1.0,0.9]
    if bolo=="FID827" : p_cut=[1.2,1.3]
    if bolo=="FID837" : p_cut=[1.5,1.3]
    if bolo=="FID838" : p_cut=[2.5,2.0]
    if bolo=="FID839" : p_cut=[2.0,1.5]
    if bolo=="FID842" : p_cut=[2.5,2.0]
    # NB: pour 841 pas de cut
    if not return_string : return p_cut
    cut_heatbiplot="( (EC2-EC1)<"+str(p_cut[0])+" && (EC1-EC2)<"+str(p_cut[1])+" )"
    fullcut="( (OWC2>0.5+2*OWC1) || "+cut_heatbiplot+")"
    return fullcut
    
################################################################################

def get_ioncorr_params(bolo,filename) :
    if bolo=="" or filename=="" : print "Params pas ok pour get_ioncorr_params"
    ff=open(filename,'r')
    lines=ff.readlines()
    ff.close()
    toto=[l for l in lines if bolo in l]
    if len(toto)!=1 : print "Pbl with",bolo,"for file",file
    ioncorr_pars=[float(x) for x in (toto[0]).split()[1:]]
    return (ioncorr_pars[0],ioncorr_pars[1])
    
################################################################################

def get_ho_rate_params(bolo,file) :
    if bolo=="" or file=="" : print "Params pas ok pour get_ho_rate_params"
    ff=open(file,'r')
    lines=ff.readlines()
    ff.close()
    toto=[l for l in lines if bolo in l]
    if len(toto)!=1 : print "Pbl with",bolo,"for file",file
    ho_pars=[float(x) for x in (toto[0]).split()[1:]]
    return ho_pars

################################################################################

def get_ho_1kev_rate(bolo,tablefile,days) :
    toto=np.loadtxt(tablefile)
    xtab=[0]
    xtab.extend(toto[:,0])
    xtab.append(600) # on "etend" le array pour pas avoir pbl interp au bord
    ytab=[toto[0,1]]
    ytab.extend(toto[:,1])
    ytab.append(toto[-1,1])
    ratefunc=interp1d(xtab,ytab)
    return ratefunc(days)

################################################################################

def get_fiducial_mass(bolo,file,return_fraction=False) :
    if bolo=="" or file=="" : print "Params pas ok pour get_fiducial_mass"
    ff=open(file,'r')
    lines=ff.readlines()
    ff.close()
    toto=[l for l in lines if bolo in l]
    if len(toto)!=1 : print "Pbl with",bolo,"for file",file
    m_full=float((toto[0]).split()[1])
    m_fid=float((toto[0]).split()[2])
    if return_fraction : return m_fid/m_full
    else : return m_fid*1.e-3 # en kg


################################################################################

def get_skimmed_axion(bolo_name,with_chiion=True,with_chichal=True, multiple="single",runtype="fond",ion_topology="All",chichaltype="all",chiion_onlypos=True,with_heatbiplotcut=True) :
    # Retourne un TChain de fichiers anaj rootifies, et avec un EventList correspondant a un ensemble de cuts
    # Cut par defaut = selection des periodes low mass ; cuts optionnels: chi2, topologie ion, multi
    # ATTENTION: retourne la chaine entiere, la selection n'etant que l'EventList associe.
    # Il faut faire ch.CopyTree() pour avoir l'arbre reellement coupe.
    
    if multiple not in ["empty","single","mult","all"] : print "Keyword multiple not ok..."
    if ion_topology not in ["All","NoSurf","FidStrict","Surf1","Surf2","Triple1","Triple2"] : print "Keyword ion_topology not ok..."

    tree,tree_file = PyRPl.open_ROOT_object("../ROOT_files/CCLyon/unblind_skim_fond_" + bolo_name + ".root", "data")
    listname=randomname()
    thelist=TEventList(listname,listname)

    # General "cSkim" cut a la Silvia, except for MULT
    cut_general="(KTH<1.5 && KTH>(0.1+0.6*TWC) && FWF<0.7 && FWC<1 && FWIA<1.5 && FWIC<1.5 && SAT==0 && FWF>0.1 && FWIA>0.1 && FWIB>0.1 && FWIC>0.1 && FWID>0.1 && FWI>0.1 && FWC>0.1 && FWIB<9.9 && FWID<9.9 && FWI<9.9)"
    # Pour 826 on jete les donnees a 8V par un cut sur JOUR
    if bolo_name=="FID826" : cut_general="("+cut_general+" && (JOUR>73.35) )"
    
    cut_multi="(MULT==1)"
    if multiple=="empty" : cut_multi="(MULT==0)"
    if multiple=="mult" : cut_multi="(MULT>1)"
    if multiple=="all" : cut_multi="(MULT>=1)"
    
    cut_chi_ion = get_chiion_cut(bolo_name,cut_only_positive=chiion_onlypos) if with_chiion else "1"
    cut_chi_chal = get_chichal_cut(bolo_name,chi2type=chichaltype) if with_chichal else "1"
    cut_heatbiplot = get_heatbiplot_cut(bolo_name) if with_heatbiplotcut else "1"

    cut_ion_topology="1"
    if ion_topology=="NoSurf" : # "Loose fiducial", 5 sigma cut on vetos (same as blinding cut)
        cut_ion_topology="(EIA<5.0*FWIA/2.35 && EIC<5.0*FWIC/2.35)"
    if ion_topology=="FidStrict" : # Complementary to surface cuts (2.7 sigma on vetos)
        cut_ion_topology="(EIA<2.7*FWIA/2.35 && EIC<2.7*FWIC/2.35)"
    if ion_topology=="Surf1" : # Silvia def du 09/07
        cut_ion_topology="(EIA>2.7*FWIA/2.35 && EIB>2.7*FWIB/2.35 && EIC<2.7*FWIC/2.35 && EID<2.7*FWID/2.35)"
    if ion_topology=="Surf2" : # Silvia def du 09/07
        cut_ion_topology="(EIA<2.7*FWIA/2.35 && EIB<2.7*FWIB/2.35 && EIC>2.7*FWIC/2.35 && EID>2.7*FWID/2.35)"
    if ion_topology=="Triple1" : # Cut at 5 sigma, orthogonal to fiducial cut : avoid double counting of 10keV with slight charge sharing
        cut_ion_topology="(EIA>5.0*FWIA/2.35 && EIB>5.0*FWIB/2.35 && EIC<5.0*FWIC/2.35 && EID>5.0*FWID/2.35)"
    if ion_topology=="Triple2" : # 
        cut_ion_topology="(EIA<5.0*FWIA/2.35 && EIB>5.0*FWIB/2.35 && EIC>5.0*FWIC/2.35 && EID>5.0*FWID/2.35)"

    full_cut = cut_general+" && "+cut_chi_ion+" && "+cut_chi_chal+" && "+cut_multi+" && "+cut_ion_topology+ " && "+cut_heatbiplot
    tree.Draw(">>" + listname, full_cut)
    num_events = thelist.GetN()

    #Store axion tree in final file 
    fout = TFile("../ROOT_files/Axion/" + bolo_name + "_skimmed_axion.root", "recreate")

    #Create skimmed tree with cuts
    t_axion = TTree("t_axion", "t_axion")

    # cSpecifiy variables to add in new tree
    list_var_names = ["EC","EFID", "EIA", "EIB", "EIC", "EID",
                       "FWC", "FWF", "FWIA", "FWIB", "FWIC", "FWID", "KTH", "JOUR"]
    list_var_arr = [array("f", [0.]) for i in range(len(list_var_names))]

    for name, arr in zip(list_var_names, list_var_arr) :
        tree.SetBranchAddress(name, arr)
        t_axion.Branch(name, arr, name + "/F")

    #Add combined heat and ionisation variable 
    Ebest = array("f", [0.])
    FWbest = array("f", [0.])
    t_axion.Branch("Ebest", Ebest, "Ebest/F")
    t_axion.Branch("FWbest", FWbest, "FWbest/F")

    for k in range(num_events):
        counter = thelist.GetEntry(k)
        tree.GetEntry(counter)

        coeff_EC = compute_coeff(tree.FWC, tree.FWF)
        FWbest[0] = compute_resolution(tree.FWC,tree.FWF)
        Ebest[0] = coeff_EC*tree.EC + (1-coeff_EC)*tree.EFID
        t_axion.Fill()

    t_axion.Write()
    fout.Close()