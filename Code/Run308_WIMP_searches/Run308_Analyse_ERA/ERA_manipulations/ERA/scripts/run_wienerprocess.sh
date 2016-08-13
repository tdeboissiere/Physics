#!/bin/sh

# Parametres
eradir=$1
anadir=$2
bolo=$3
run=$4
overwrite=1
#listevoies=(Chal1 Chal2 Col1 Col2) # en dur pour l'instant
listevoies=(Chal1 Chal2 Col1 Col2 Vet1 Vet2) # en dur pour l'instant
# aussi implicitement en dur: marche uniquement pour fichiers "tmpwiener_voie"
# et pour une seule voie par sous-job
# (eg. pas ok pour lancer WienerProcess sur plusieurs voies ni sur ampls ntds)
voie=${listevoies[${SGE_TASK_ID}-1]}

## Copie fichiers de donnees
# Reproduction arborescence anadir
mkdir ${TMPDIR}/AnaDir
cp ${anadir}/*.txt ${TMPDIR}/AnaDir/.
mkdir ${TMPDIR}/AnaDir/${bolo}
cd ${TMPDIR}/AnaDir/${bolo}
cp ${anadir}/${bolo}/*.txt .
mkdir Amplitudes Traces Spectra Figures
# WienerProcess : copie fichiers traces et spectres
cp ${anadir}/${bolo}/Traces/${run}*.root ${TMPDIR}/AnaDir/${bolo}/Traces/.
cp ${anadir}/${bolo}/Spectra/spectra_${run}_${bolo}.root ${TMPDIR}/AnaDir/${bolo}/Spectra/.
#cp ${anadir}/${bolo}/Amplitudes/*_${run}_*_${bolo}.root ${TMPDIR}/AnaDir/${bolo}/Amplitudes/.

## Fichiers de params local
paramfile=${TMPDIR}/params.txt
echo "AnaDir = "${TMPDIR}"/AnaDir" >> ${paramfile}
echo "BoloName = "${bolo} >> ${paramfile}
echo "OverWrite = "${overwrite} >> ${paramfile}
echo "RunName = "${run} >> ${paramfile}
echo "ChannelSelection = "${voie} >> ${paramfile}
more ${paramfile}

## Lancement programme
echo "** Etape: WienerProcess.exe"
${eradir}/bin/WienerProcess.exe ${paramfile}
echo "** Done."

## Copie fichiers de sortie
cp ${TMPDIR}/AnaDir/${bolo}/Amplitudes/tmpwiener_${run}_${voie}_${bolo}.root ${anadir}/${bolo}/Amplitudes/.
