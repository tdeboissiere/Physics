
def ReadGlobalParams(file):
    params=dict()
    f=open(file,"r")
    lines=f.readlines()
    for theline in lines :
        theline=theline.strip() # degage le "\n" a la fin
        toto=theline.split(' ')
        toto=[x for x in toto if x!='']
        if len(toto)==3 and toto[1]=="=" and (toto[0])[0]!="#":
            params[toto[0]] = toto[2]#.rstrip("\n")
        if len(toto)>0 :
            if toto[0].find("=")!=-1 : print "Warning, error in param file syntax?"
    return params
