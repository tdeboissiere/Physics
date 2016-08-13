#!/bin/sh

eradir=$1
anadir=$2
bolo=$3

overwrite=0
export PYTHONPATH=${PYTHONPATH}:${eradir}

paramfile=${TMPDIR}/params.txt
echo "eradir = "${eradir} >> ${paramfile}
echo "anadir = "${anadir} >> ${paramfile}
echo "bolo = "${bolo} >> ${paramfile}
echo "overwrite = "${overwrite} >> ${paramfile}

more ${paramfile}
echo "** Etape: BuildPeriodList.py"
${eradir}/python/BuildPeriodList.py ${paramfile}
echo "** Done."

