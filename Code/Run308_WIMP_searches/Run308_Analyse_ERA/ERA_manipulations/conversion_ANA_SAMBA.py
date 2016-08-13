#Definitions that will help you navigate the  different IDs.


#The following two auxiliary functions convert lower case letters into numbers and viceversa
#1 to 'a', 2 to 'b', etc.

def int2chr(n):
	return chr(96+n)

def chr2int(l):
	return ord(l)-96


#The following dictionary gives the mac number for each bolo.
#Please check if the numbers are correct, they might change.

bolo_mac = {
"FID807" : 3,
"FID810" : 1,
"FID817" : 1,
"FID820" : 3,
"FID821" : 3,
"FID822" : 3,
"FID823" : 1,
"FID824" : 1,
"FID825" : 1,
"FID826" : 1,
"FID827" : 1,
"FID828" : 1,
"FID831" : 2,
"FID832" : 3,
"FID837" : 2,
"FID838" : 2,
"FID839" : 2,
"FID840" : 3,
"FID841" : 3,
"FID842" : 3,
"FID843" : 3,
"FID844" : 3,
"FID845" : 3,
"FID846" : 3
}

#the following function makes two dictionaries from two text files with the same number of lines.

def make_dict(file1, file2):

	[list1, list2] = [[line.strip() for line in open(file1)], [line.strip() for line in open(file2)]]
	dict1={}
	dict2={}
	for i in range(len(list1)):
		dict1[list1[i]] = list2[i]
		dict2[list2[i]] = list1[i]
	return (dict1, dict2)



#Conversion between ANA and Samba. Note that ANA has less information than Samba concerning the date, so this might not stay relevant for ever.

#The following function takes an ANA run number and a machine number, and gives back the samba run number.

def ana2samba(run_ana):
	mac = run_ana / 100000	#Mac number
	r0 = run_ana % 100000
	B = r0 / 10000			#Month
	
	r1 = r0 % 10000
	CC= r1 / 100		
	DD= r1 % 100
	
	A = "o"		#Year, hardcoded here	

	run_samba = A + int2chr(B+5) + str(CC).zfill(2) + int2chr(mac) + str(DD).zfill(3)

	return run_samba


#The following function takes a samba run number, and gives back the ANA run number and machine number.
	
def samba2ana(run_samba):
	A = run_samba[0]			#Year, note it doesn't appear later
	if A != "o":
		print "this run is not from the year this function was programmed for"
	B = chr2int(run_samba[1])	#Month
	
	CC = int(run_samba[2:4])
	mac = chr2int(run_samba[4])	#Mac number
	DD = int(run_samba[5:8])	#While DD technically has three digits, the first one is always zero
	
	run_ana = mac*100000 + (B-5)*10000 + CC*100 + DD
	
	return run_ana
