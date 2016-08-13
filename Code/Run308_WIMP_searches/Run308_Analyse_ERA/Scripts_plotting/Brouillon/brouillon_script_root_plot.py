    # ######################
    # # ROOT plot
    # #####################
    # hheatonly= TH2F("heatonly", "heatonly", 100, 0, 100, 100, -0.2, 1.2)
    # hfid= TH2F("fid", "fid", 100, 0, 100, 100, -0.2, 1.2)    
    # hvet= TH2F("vet", "vet", 100, 0, 100, 100, -0.2, 1.2)

    # tree.Project("fid", "Q:ER", standard_cuts + fidcut)
    # tree.Project("vet", "Q:ER", standard_cuts + vetcut)
    # tree.Project("heatonly", "Q:ER", standard_cuts + heatcut)
    
    # PyRPl.process_TH2(hheatonly, X_title = "Recoil energy (keV)", Y_title = "Ionisation yield", color = kRed, marker_style = 20, marker_size = 0.6)
    # PyRPl.process_TH2(hfid, X_title = "Recoil energy (keV)", Y_title = "Ionisation yield", color = kBlue, marker_style = 20, marker_size = 0.6)
    # PyRPl.process_TH2(hvet, X_title = "Recoil energy (keV)", Y_title = "Ionisation yield", color = kGreen+3, marker_style = 20, marker_size = 0.6)

    # cc=TCanvas("cc","cc")
    # hvet.Draw() 
    # hfid.Draw("same") 
    # hheatonly.Draw("same") 

    # raw_input()
    # cc.Print(figure_path + bolo_name + "_Qplot.png")

    #Other root plot

    # hheatonly= TH2F("heatonly", "heatonly", 100, -1, 10, 100, -1, 10)
    # hNTD= TH2F("NTD", "NTD", 100, -1, 10, 100, -1, 10)
    # tree.Project("NTD", "EC2:EC1", cut_line + "")
    # tree.Project("heatonly", "EC2:EC1", cut_line + "&&XOC1<0.114&&XOC2<0.114&&abs(EC1-EC2)<1")
    
    # PyRPl.process_TH2(hheatonly, X_title = "EC1 (keV)", Y_title = "EC2 (keV)", color = kRed)
    # PyRPl.process_TH2(hNTD, X_title = "EC1 (keV)", Y_title = "EC2 (keV)", color = kBlack)

    # cc=TCanvas("cc","cc")
    # hNTD.Draw("")   
    # hheatonly.Draw("same")   

    # raw_input()
    # figure_path = script_utils.create_directory("./Figures/" + bolo_name + "/")
    # cc.Print(figure_path + bolo_name + "_heatonly_NTD.png")