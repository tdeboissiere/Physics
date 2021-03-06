
import csv

def MacLetter(mac) :
    if len(mac)!=2 :
        print "MacLetter error: wrong mac name ",mac
    # Recette magique pour avoir le lien chiffre-lettre:
    letterlist=[chr(i) for i in xrange(ord('a'), ord('j')+1)]
    return letterlist[int(mac[1])-1]

def GetConsigne_Log(fichierlog,voie,version) :
    # recupere la polar ou le D2 pour une voie donnee
    # a partir du log (pour les vieux fichiers samba)
    # TODO: aussi pour les nouveaux fichiers: cross-check des polars appliquees..
    consigne=0
    f=open(fichierlog,"r")
    detect_samba=voie[voie.rfind(' ')+1:]
    type_bolo="GeNTD" # on en a besoin pour lire les polars...
    if detect_samba[0:2]=="ID" or detect_samba[0:3]=="FID" :
        type_bolo="InterDigit"
    typevoie_samba="NONE"
    if voie.find("centre") != -1 :
        typevoie_samba="centre"
    elif voie.find("garde") != -1 :
        typevoie_samba="garde"
    elif voie.find("chaleur") != -1 :
        typevoie_samba="chaleur"
    else : print "Pbl type de voie"
    lettre_annee=fichierlog[fichierlog.rfind("/")+1:fichierlog.rfind("/")+2]
    pattern_cherche=""
    # On reprend les "oldies" crasseuses du code EdwRawDB.cc
    # ca doit marcher pour l'essentiel (ID/GeNTD, Run 10/12...)
    # gestion polar consigne:
    #    // Change in polar consigne starting run 11.
    #  // Before run 11: centre = polar-voie3, garde = polar-ionis
    #  // After run 11: centre = polar-ionis, garde = polar-voie3
    # TODO = TESTER CA SUR DES ID RUN10/12 !!!!
    if (typevoie_samba=="chaleur") : pattern_cherche="d2 ="
    elif (float(version)<9 and type_bolo=="GeNTD") : pattern_cherche="polar-ionis ="
    elif (float(version)>=9 and type_bolo=="GeNTD") : pattern_cherche="polar-centre ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="centre" and float(version)<9 and lettre_annee<="i") :
        pattern_cherche="polar-voie3 ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="centre" and float(version)<9 and lettre_annee>"i") :
        pattern_cherche="polar-ionis ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="centre" and float(version)>=9) :
        pattern_cherche="polar-centre ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="garde" and float(version)<9 and lettre_annee<="i") :
        pattern_cherche="polar-ionis ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="garde" and float(version)<9 and lettre_annee>"i") :
        pattern_cherche="polar-voie3 ="
    elif (type_bolo=="InterDigit" and typevoie_samba=="garde" and float(version)>=9) :
        pattern_cherche="polar-garde ="
    ok=1
    while ok :
        line=f.readline()
        if line.find("(detecteur: "+detect_samba+")") != -1 :
            ok=0
            while line.find("initialisation terminee") == -1 :
                line=f.readline()
                if (line.find(pattern_cherche)!=-1) :
                    toto=line[line.find("=")+2:]
                    consigne=float(toto)
        if line.find("Creation du fichier") > 0 :
            ok=0
    f.close()        
    return consigne

def GetConsigne_Partition(partition,voie) :
    # Recupere la polar ou le D2 pour une voie donnee.
    # a partir d'un fichier _000
    consigne=0
    f=open(partition,"r")
    detect_samba=voie[voie.rfind(' ')+1:]
    typevoie_samba="NONE"
    if voie.find("centre") != -1 :
        typevoie_samba="centre"
    elif voie.find("garde") != -1 :
        typevoie_samba="garde"
    elif voie.find("chaleur") != -1 :
        typevoie_samba="chaleur"
    else : print "Pbl type de voie"
    ok=1
    while ok :
        line=f.readline()
        if line.find("Detecteur "+detect_samba) != -1 :
            ok=0
            line=f.readline()
            while line.find("Detecteur") == -1 :
                line=f.readline()
                if ( (typevoie_samba=="centre" and line.find(" polar-centre :=")!=-1) or (typevoie_samba=="garde" and line.find(" polar-garde :=")!=-1) or (typevoie_samba=="chaleur" and line.find(" d2 :=")!=-1) ) :
                    toto=line[line.find(":=")+3:line.rfind("}")-1]
                    consigne=float(toto)
        if line.find("Entete de run") > 0 :
            ok=0
    f.close()
    return consigne


def GetConsigne_CSV(csvfile,voie) :
    # Recupere la polar ou le D2 pour une voie donnee.
    # A partir du fichier setup.csv
    consigne=-10000
    reader = csv.reader(open(csvfile, 'r'),delimiter=';')
    row=reader.next()
    row=[x.strip() for x in row if x != '']
    # souvent il y a un ";" au debut, en fonction de la version de samba...
    if row[1]!="detecteur" or row[3]!="capteur" or row[8]!="reglage" :
        print "Warning, pbl in CSV file format.. :",row
    for row in reader :
        row=[x.strip() for x in row if x!= '']
        if 'NEANT' in row : continue # nouveau samba : voies en doublon dans le csv avec une ligne "neant"...
        detecteur=row[1]
        capteur=row[3]
        if (capteur+" "+detecteur) == voie :
            if capteur.find("chal")!= -1 :
                if 'D2' not in row : print "D2 not found in",voie,"- file:",csvfile
                pos = row.index('D2')
                consigne=int(row[pos+1])
            elif capteur.find("ionis")!=-1 :
                lettre=capteur[len(capteur)-1]
                if 'polar-'+lettre not in row : print "polar not found in",voie,"- file:",csvfile
                pos = row.index('polar-'+lettre)
                consigne=float(row[pos+2])
    if consigne==-10000 : print "Warning: not found consigne in CSV file for",voie," - file:",csvfile
    return consigne


def GetConsigne(sambadir,run,voie) :
    coefficient=0
    partition=sambadir+run+"/"+run+"_000"
    f=open(partition,"r")
    line=f.readline()
    if line != "* Archive SAMBA\n" :
        print "Est-ce bien un fichier samba?"
    line=f.readline()
    version=line[line.rfind(" ")+1:line.rfind("\n")]
    f.close()
    if float(version) >= 9.29 :
        coefficient=GetConsigne_CSV(sambadir+run+"/"+run+"_setup.csv",voie)
    elif version=="9.11" :
        coefficient=GetConsigne_Partition(sambadir+run+"/"+run+"_000",voie)
    elif version=="8.15" :
        coefficient=GetConsigne_Log(sambadir+run+"/"+run+"_log",voie,version)
    else :
        print "Version de samba pas testee pour lire consigne!",version
    return coefficient

def GetRunType_Partition(partition) :
    f=open(partition,"r")
    line=f.readline()
    if line != "* Archive SAMBA\n" :
        print "Est-ce bien un fichier samba?"
    line=f.readline()
    version=line[line.rfind(" ")+1:line.rfind("\n")]
    liste_versions_ok=["8.15","9.11","9.29","9.31","9.32","9.33","9.34","9.35"]
    if version not in liste_versions_ok :
        print "Version de samba pas testee pour lire type de run!",version
    while line.find("Entete de run") == -1 :
        line=f.readline()
    while line.find("Entete d'evenement") == -1  :
        line=f.readline()
        if line[0:6]=="Type =" :
            toto=(line[line.find("=")+1:line.find("#")]).strip()
            # Ajout Run300: type de run "b133"... (bug samba?)
            if toto=="ba133" or toto=="b133" : runtype="Gamma"
            elif toto=="co60" : runtype="Neutron"
            elif toto=="fond" : runtype="Bckgd"
            # A completer... !!! (en utilisant "condition")
            else : runtype="Undefined"
            break
    else :
        print "Run type not found"
        runtype="Undefined"
    f.close()
    return runtype
