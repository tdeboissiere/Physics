#!/bin/sh

# Parametres
eradir=$1
anadir=$2
bolo=$3
run=$4
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
# BuildNoiseSpectra : copie fichiers traces
cp ${anadir}/${bolo}/Traces/${run}*.root ${TMPDIR}/AnaDir/${bolo}/Traces/.

## Fichiers de params local
paramfile=${TMPDIR}/params.txt
echo "AnaDir = "${TMPDIR}"/AnaDir" >> ${paramfile}
echo "BoloName = "${bolo} >> ${paramfile}
echo "RunName = "${run} >> ${paramfile}
echo "OverWrite = "${overwrite} >> ${paramfile}
more ${paramfile}

## Lancement programme
echo "** Etape: BuildNoiseSpectra.exe"
${eradir}/bin/BuildNoiseSpectra.exe ${paramfile}
echo "** Done."

## Copie fichiers de sortie
cp ${TMPDIR}/AnaDir/${bolo}/Spectra/spectra_${run}_${bolo}.root ${anadir}/${bolo}/Spectra/.

