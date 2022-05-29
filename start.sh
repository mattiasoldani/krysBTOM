#!/bin/bash

# preliminary settings ############################################################################

FIRSTEXEC=1  # if 1 (0), first cycle is executed automatically at startup

RUNSTRL=1  # filename part w/ run number, left index
RUNSTRR=9  # filename part w/ run number, right index
	
TOSERVER=1  # if 1 (0), data are (not) sent to server on EOS
SERVERPATH="$HOME/eos_sshfs_vm/www/krysBTOM"  # web server location on EOS (should be mounted via sshfs)
SYNCBACKSERVER=1  # if 1, sync website back from server after sending data (subject to TOSERVER)
	
###################################################################################################

echo "---------------------------"
echo "krysBTOM: here we go again!"
echo "---------------------------"

# read cycle parameters & input data path from configuration files
# DATAPATH: path with data
# TCYCLE: interval time between checks for new files in minutes
# WAITAFTERCYCLE: time to wait after each (positive) cycle in seconds
# NFILES: nr. of files per cycle
# EXCLUDECOMPILER: if 1 (0), recompilation of the analysis software is (not) performed
# NRUN: number of run to be processed; if 0, process the latest (last spill of latest run is always processed)
read DATAPATH < 'inFileDataPath.conf'
	# note: the format of this file has to be checked carefully: only 1 line, 1 entry separated by a space
read TCYCLE WAITAFTERCYCLE NFILES EXCLUDECOMPILER NRUN < 'confCycle.conf'
	# note: the format of this file has to be checked carefully: only 1 line, 3 entries separated by a space

if [ $FIRSTEXEC -eq 1 ] ; then
	source cycle.sh $RUNSTRL $RUNSTRR $TOSERVER $SYNCBACKSERVER $SERVERPATH
else
	echo "waiting for new data..."
	echo "---"
fi
	
while [ true ] ; do
	if [ ! -z "$(find $DATAPATH -mmin $TCYCLE)" ] ; then
		source cycle.sh $RUNSTRL $RUNSTRR $TOSERVER $SYNCBACKSERVER $SERVERPATH
		sleep $WAITAFTERCYCLE
	else
		NOW="$(date +'%Y %b %d, %X')"
		STRNODATA="no new data @ $NOW, sticking to older files..."
		echo $STRNODATA
		echo $STRNODATA > "./outtext/infoCycle.conf"
		source talkToServer.sh $TOSERVER $SYNCBACKSERVER $SERVERPATH
		echo "---"
		sleep 1
	fi
done

