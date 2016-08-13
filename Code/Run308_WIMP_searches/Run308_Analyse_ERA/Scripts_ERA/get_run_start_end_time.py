from ROOT import TH1F, TFile

import script_utils as script_utils



def get_run_start_end_time(bolo_name, data_dir, tree_name ):
    """
    Creates a .txt file with the start and end time information for all polars of a given bolo

    Detail:

    Arguments:
    bolo_name (str) the bolometer name
    data_dir  (str) the data directory (containing Jules n-tuple)
    tree_name (str) the tree name in the data file

    Outputs:

    A .txt file bolo_name + "_polar_start_and_end_time_prelim.txt", "w"
    it contains the polarisation name, the start and end time, for each polar

    """    
    num_polar  = script_utils.get_howmany(bolo_name)
    
    #Load the data
    file_tree  = TFile(data_dir+bolo_name+"_fond.root")
    tree       = file_tree.Get(tree_name)
    nEntries   = tree.GetEntries()
    
    #Initialize the start and end time 
    tree.GetEntry(0)
    all_tree_start_time =1E6*tree.UT1 + tree.UT2
    all_tree_end_time   =1E6*tree.UT1 + tree.UT2
    
    #Find the start and end time of the whole data tree
    for i in range(1,nEntries):
        tree.GetEntry(i)
        time= 1E6*tree.UT1 + tree.UT2
        if time < all_tree_start_time:
            all_tree_start_time = float(time)
        if time > all_tree_end_time:
            all_tree_end_time = float(time)

    #use a fine binning for the histograms that will be used to
    #find the start and end time 
    bin = 2*int(all_tree_end_time-  all_tree_start_time)

    #Define path for the file. Create the directory if it does not exist
    path_name=script_utils.create_directory('../Analyse_' + bolo_name + '/Text_files/')    
    
    #Create output file
    file_polar = script_utils.open_text_file(path_name, bolo_name + "_all_polars_with_entries.txt", "r")
    
    #Read the file lines
    list_lines_polar = [elem.rstrip().split(",") for elem in file_polar.readlines()]
    
    #Define and Fill in the polar cuts and an array of titles for the histograms
    list_cut_polar, list_name_polar = [], []
    for i in range(1,len(list_lines_polar)):
        list_cut_polar.append(list_lines_polar[i][1])
        list_name_polar.append(list_lines_polar[i][0])

    #Project the data to match the polar cut
    list_hist=[TH1F() for i in range(num_polar)]
    htemp= TH1F("temp","temp" ,bin,all_tree_start_time,all_tree_end_time) #Temporary histogram to load the histogram array
    for k in range(num_polar):
        print list_cut_polar[k]
        tree.Project("temp","1E6*UT1+UT2", list_cut_polar[k])
        list_hist[k] = htemp.Clone("h"+str(k))

    #Define and fill a list of starting and ending indices for each polar condition
    #The start index is defined as the first non zero entry
    #The end index is defined as the last non zero entry
    list_start_index=[]
    list_end_index=[]

    for k in range(num_polar):
        counter=0
        while (list_hist[k].GetBinContent(counter) ==0):
            counter+=1
        list_start_index.append(counter)

    for k in range(num_polar):
        counter=int(bin)
        while (list_hist[k].GetBinContent(counter) ==0):
            counter+=1
        list_end_index.append(counter)  

    #Create output file
    outputfile = script_utils.open_text_file(path_name , bolo_name +"_polar_start_and_end_time_prelim.txt", "w")
    for k in range(num_polar):
        start_and_end_line=str(list_name_polar[k]) + ","+ str(list_hist[k].GetBinCenter(list_start_index[k]))+ "," + str(list_hist[k].GetBinCenter(list_end_index[k]))
        outputfile.write(start_and_end_line + "\n")

        #Print file_content
        script_utils.print_utility(start_and_end_line) 
    outputfile.close()