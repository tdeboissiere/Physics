from ROOT import *


fheat = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Heatonly/ROOT_files/FID837_heatonly_tree.root", "read")
fheat_above2 = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Heatonly/ROOT_files/FID837_heatonly_tree_above2.root", "read")
fheat_below2 = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Heatonly/ROOT_files/FID837_heatonly_tree_below2.root", "read")

fheat_above2_true = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Heatonly/ROOT_files/FID837_heatonly_above2_from_true_events_tree.root", "read")
fheat_below2_true = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Heatonly/ROOT_files/FID837_heatonly_below2_from_true_events_tree.root", "read")

fFidGamma = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Gamma/ROOT_files/FID837_Gamma_fid_tree.root", "read")
fS1Gamma = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Gamma/ROOT_files/FID837_Gamma_surf1_tree.root", "read")
fS2Gamma = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Gamma/ROOT_files/FID837_Gamma_surf2_tree.root", "read")

fS1Beta = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Beta_and_Pb/ROOT_files/FID837_Beta_surf1_tree.root", "read")
fS2Beta = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Beta_and_Pb/ROOT_files/FID837_Beta_surf2_tree.root", "read")

fS1Pb = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Beta_and_Pb/ROOT_files/FID837_Pb_surf1_tree.root", "read")
fS2Pb = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/Beta_and_Pb/ROOT_files/FID837_Pb_surf2_tree.root", "read")

ftrue = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/True_events/ROOT_files/FID837_true_events_tree.root", "read")

fWIMP = TFile("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_FID837/standard_resolution/WIMP/Events/FID837_WIMP_mass_10_tree.root", "read")

theat = fheat.Get("t_new")
theat_above2 = fheat_above2.Get("t_new")
theat_below2 = fheat_below2.Get("t_new")

theat_above2_true = fheat_above2_true.Get("t_new")
theat_below2_true = fheat_below2_true.Get("t_new")

tFidGamma = fFidGamma.Get("t_new")
tS1Gamma= fS1Gamma.Get("t_new")
tS2Gamma= fS2Gamma.Get("t_new")

tS1Beta= fS1Beta.Get("t_new")
tS2Beta= fS2Beta.Get("t_new")

tS1Pb= fS1Pb.Get("t_new")
tS2Pb= fS2Pb.Get("t_new")

ttrue= ftrue.Get("t_new")

tWIMP= fWIMP.Get("t1")


# list_tree = [theat, theat_above2, theat_above2_true, theat_below2, theat_below2_true, tFidGamma, tS1Gamma, tS2Gamma, tS1Beta, tS2Beta, tS1Pb, tS2Pb, tWIMP, ttrue]
list_tree = [theat, tFidGamma, tS1Gamma, tS2Gamma, tS1Beta, tS2Beta, tS1Pb, tS2Pb, tWIMP, ttrue]

for elem in list_tree:
	print elem.GetEntries()

# for elem in list_tree:
# 	cc= TCanvas("cc", "cc")
# 	elem.Draw("0.5*(EC1+EC2)>>hist(200,-1,20)")
# 	raw_input()
# 	del cc
# 	del hist


for elem in list_tree:
	cc= TCanvas("cc", "cc")
	elem.Draw("0.5*(EIB+EID):0.5*(EC1+EC2)>>hist(200,-1,20, 200, -1, 20)")
	hist.SetContour(40)
	hist.Draw("cont4z")
	raw_input()
	del cc
	del hist