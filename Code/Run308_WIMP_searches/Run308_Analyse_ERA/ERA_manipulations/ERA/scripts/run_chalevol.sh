#!/bin/sh

eradir=$1
anadir=$2
bolo=$3
voie=$4

overwrite=0
interactif=0
export PYTHONPATH=${PYTHONPATH}:${eradir}

paramfile=${TMPDIR}/params.txt
echo "eradir = "${eradir} >> ${paramfile}
echo "anadir = "${anadir} >> ${paramfile}
echo "bolo = "${bolo} >> ${paramfile}
echo "overwrite = "${overwrite} >> ${paramfile}
echo "interactif = "${interactif} >> ${paramfile}
echo "voie = "${voie} >> ${paramfile}

more ${paramfile}
echo "** Etape: GainChalEvol.py"
${eradir}/python/GainChalEvol.py ${paramfile}
echo "** Etape: BaselineEvol.py"
${eradir}/python/BaselineEvol.py ${paramfile}
echo "** Done."

