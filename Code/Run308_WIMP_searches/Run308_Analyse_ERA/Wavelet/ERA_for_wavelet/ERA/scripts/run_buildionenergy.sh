#!/bin/sh

eradir=$1
anadir=$2
bolo=$3

overwrite=0
export PYTHONPATH=${PYTHONPATH}:${eradir}

paramfile_py=${TMPDIR}/params_py.txt
echo "eradir = "${eradir} >> ${paramfile_py}
echo "anadir = "${anadir} >> ${paramfile_py}
echo "bolo = "${bolo} >> ${paramfile_py}
echo "overwrite = "${overwrite} >> ${paramfile_py}
more ${paramfile_py}

paramfile_era=${TMPDIR}/params_era.txt
echo "AnaDir = "${anadir} >> ${paramfile_era}
echo "BoloName = "${bolo} >> ${paramfile_era}
echo "OverWrite = "${overwrite} >> ${paramfile_era}
more ${paramfile_era}


echo "** Etape: BuildIonEnergy.py"
${eradir}/python/BuildIonEnergy.py ${paramfile_py}
echo "** Etape: ClassifyEvent.exe"
${eradir}/bin/ClassifyEvent.exe ${paramfile_era}
echo "** Done."

