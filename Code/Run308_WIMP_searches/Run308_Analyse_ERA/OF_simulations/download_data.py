#Combine all bolo data to build a spectrum

import os
import paramiko
import getpass
import script_utils as script_utils

def download_ANA_data(bolo_name, lyon_ANA_dir):

	"""Quick description
	
	Detail:
		Download ANA files

	Args:
		bolo_name = (str) bolo name 
		lyon_ANA_dir  = (str) directory from which ANA files are downloaded

		
	Returns:
		void

	Raises:
		void
	"""

	data_dir = "../Fond/"

	ssh = paramiko.SSHClient() 
	ssh.load_host_keys(os.path.expanduser(os.path.join('~', '.ssh', 'known_hosts')))
	pwd = getpass.getpass('Password please:   ')  # command line prompt without echo
	ssh.connect('ccage.in2p3.fr', username='tmain', password=pwd)
	stdin, stdout, stderr = ssh.exec_command("cd " + lyon_ANA_dir + " && ls *"+ bolo_name[3:]+ "*")

	#list containing the file names
	liste=stdout.read().splitlines()	

	#print  stderr.read.splitlines()
	sftp = ssh.open_sftp()
	for bolo_file in liste:
	    if ('-fond-' in bolo_file and 'root' in bolo_file): 
			script_utils.print_utility(script_utils.COL("Adding file " + bolo_file, "blue"))
			sftp.get(lyon_ANA_dir+bolo_file,data_dir+bolo_file)
	sftp.close()
	ssh.close()


def download_spectra_files(bolo_name, lyon_dir):

	"""Quick description
	
	Detail:
		Download ERA Spectra files

	Args:
		bolo_name = (str) bolo name 
		lyon_dir = (str) directory from which files are downloaded

		
	Returns:
		void

	Raises:
		void
	"""

	amp_dir = lyon_dir + "Spectra/"
	data_dir = "../Spectra_files/" + bolo_name + "/"

	ssh = paramiko.SSHClient() 
	ssh.load_host_keys(os.path.expanduser(os.path.join('~', '.ssh', 'known_hosts')))
	pwd = getpass.getpass('Password please:   ')  # command line prompt without echo
	ssh.connect('ccage.in2p3.fr', username='tmain', password=pwd)
	stdin, stdout, stderr = ssh.exec_command("cd " + amp_dir + " && ls *"+ bolo_name+ "*")

	#list containing the file names
	liste=stdout.read().splitlines()	

    #Create directory for data if they don't exist
	script_utils.create_directory(data_dir)

	#print  stderr.read.splitlines()
	sftp = ssh.open_sftp()
	for bolo_file in liste:
		script_utils.print_utility(script_utils.COL("Adding file " + bolo_file, "blue"))
		sftp.get(amp_dir+bolo_file,data_dir+bolo_file)
	sftp.close()
	ssh.close()


def download_txt_files(bolo_name, lyon_dir):

	"""Quick description
	
	Detail:
		Download txt files from bolo processing

	Args:
		bolo_name = (str) bolo name 
		lyon_dir = (str) directory from which files are downloaded

		
	Returns:
		void

	Raises:
		void
	"""

	data_dir = "./Text_files/" + bolo_name + "/"

	ssh = paramiko.SSHClient() 
	ssh.load_host_keys(os.path.expanduser(os.path.join('~', '.ssh', 'known_hosts')))
	pwd = getpass.getpass('Password please:   ')  # command line prompt without echo
	ssh.connect('ccage.in2p3.fr', username='tmain', password=pwd)
	stdin, stdout, stderr = ssh.exec_command("cd " + lyon_dir + " && ls *.txt*")

	#list containing the file names
	liste=stdout.read().splitlines()	

    #Create directory for data if they don't exist
	script_utils.create_directory(data_dir)

	#print  stderr.read.splitlines()
	sftp = ssh.open_sftp()
	for bolo_file in liste:
		script_utils.print_utility(script_utils.COL("Adding file " + bolo_file, "blue"))
		sftp.get(lyon_dir+bolo_file,data_dir+bolo_file)
	sftp.close()
	ssh.close()