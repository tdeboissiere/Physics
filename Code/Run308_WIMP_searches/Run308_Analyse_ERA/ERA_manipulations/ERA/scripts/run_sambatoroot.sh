#!/bin/sh

# Parametres
eradir=$1
anadir=$2
sambadir=$3
bolo=$4
run=$5
overwrite=0

## Compilation en local
#mkdir ${TMPDIR}/ERA
#cp -r ${eradir}/* ${TMPDIR}/ERA/.
#cd ${TMPDIR}/ERA
#make clean
#make

## Copie fichiers de donnees
# Reproduction arborescence anadir
mkdir ${TMPDIR}/AnaDir
cp ${anadir}/*.txt ${TMPDIR}/AnaDir/.
mkdir ${TMPDIR}/AnaDir/${bolo}
cd ${TMPDIR}/AnaDir/${bolo}
cp ${anadir}/${bolo}/*.txt .
mkdir Amplitudes Traces Spectra Figures
# SambaToRoot: copie donnees samba
mkdir ${TMPDIR}/SambaDir
cp -rv ${sambadir}/${run} ${TMPDIR}/SambaDir/.

## Fichiers de params local
paramfile=${TMPDIR}/params.txt
echo "AnaDir = "${TMPDIR}"/AnaDir" >> ${paramfile}
echo "SambaDir = "${TMPDIR}"/SambaDir" >> ${paramfile}
echo "BoloName = "${bolo} >> ${paramfile}
echo "RunName = "${run} >>  ${paramfile}
echo "OverWrite = "${overwrite} >> ${paramfile}
more ${paramfile}

## Lancement programme
echo "** Etape: SambaToRoot.exe"
${eradir}/bin/SambaToRoot.exe ${paramfile}
echo "** Done."

## Copie fichiers de sortie
cp ${TMPDIR}/AnaDir/${bolo}/Traces/${run}* ${anadir}/${bolo}/Traces/.
