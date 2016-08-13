from os import listdir
from os.path import isfile, join
from ROOT import *



def fusion_test_files(bolo_name, data_dir="../Fond/", tree_name="data"):

 
    #List of files in data_dir
	list_files_all = [ f for f in listdir(data_dir) if isfile(join(data_dir,f)) ]        

	list_files_bolo =sorted([file_bolo for file_bolo in list_files_all if bolo_name.lower() in file_bolo])
	chain=TChain(tree_name,"eionbis")
	for partition in list_files_bolo:
		print partition
		chain.AddFile(data_dir+partition)

	fusion_file_path = '../Fond_ERA_merged/'
	fusion_file_name = fusion_file_path+bolo_name.upper()+"_fond.root"
	# fusion_file      = TFile(fusion_file_name,"recreate")
	chain.Merge(fusion_file_name)
	# fusion_file.Close()    
	    

