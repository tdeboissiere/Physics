#This script copies files from a directory to another directory. Can give a pattern to avoid
#It uses a bottom up approach to do that.


import os

def has_file_pattern(file_name, list_pattern_to_avoid):

    #Check if any pattern of list_pattern_to_avoid is in file_name
    #Returns False if a pattern to be avoided is in the file_name
    for pattern in list_pattern_to_avoid:
        if pattern in file_name:
            return True
    
    return False


def delete_files(source_dir,  list_pattern_to_avoid):





    for root, subFolders, files in os.walk(source_dir,topdown=True):
        folder_path=root[len(source_dir):]
        if (len(folder_path)!=0 and folder_path[-1] != '/' ): folder_path+='/'
        for file_name in files:
            if has_file_pattern(file_name, list_pattern_to_avoid):
                os.remove(source_dir+folder_path+file_name)




# source_dir='/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_FID837_only/Analyse_FID837/'
# target_dir='/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_FID837_only/'
# delete_files(source_dir, target_dir, ['.eps' ])

