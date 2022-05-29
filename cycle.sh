#!/bin/bash

# data processing cycle
# 5 arguments needed: 
# - left & right character position limit for the sun-identifying substring in file names, 
# - boolean for enable copy to server, 
# - boolean for enable sync from server, 
# - server location (either remote or mounted locally)
RUNSTRL="$1"
RUNSTRR="$2"
TOSERVER="$3"
SYNCBACKSERVER="$4"
SERVERPATH="$5"

# read cycle parameters & input data path from configuration files
# DATAPATH: path with data
# TCYCLE: interval time between checks for new files in minutes
# WAITAFTERCYCLE: time to wait after each (positive) cycle in seconds
# NFILES: nr. of files per cycle
# EXCLUDECOMPILER: if 0 (1), recompilation of the analysis software is (not) performed
# NRUN: number of run to be processed; if 0, process the latest (last spill of latest run is always processed)
read DATAPATH < 'inFileDataPath.conf'
	# note: the format of this file has to be checked carefully: only 1 line, 1 entry separated by a space
read TCYCLE WAITAFTERCYCLE NFILES EXCLUDECOMPILER NRUN < 'confCycle.conf'
	# note: the format of this file has to be checked carefully: only 1 line, 3 entries separated by a space

# datetime @ cycle start --> writing it to text file
NOW="$(date +'%Y %b %d, %X')"
printf "welcome to the cycle! current datetime: %s\n" "$NOW"

# updating text files with cycle parameters & datetime
echo "last cycle was @ $NOW" > "./outtext/infoCycle.conf"
echo "checking for new files every $TCYCLE min, $NFILES files per cycle" >> "./outtext/infoCycle.conf"

# <<comment # uncomment this (& other "comment" below) & create a lsFiles.conf manually to use the latter
# list of files to process
echo "will work with up to $NFILES files..."
if [ $NRUN -eq 0 ] ; then 
	LSFILESID=$(ls $DATAPATH | tail -1 | cut -c$RUNSTRL-$RUNSTRR)
	ls $DATAPATH | grep $LSFILESID | tail -$NFILES > "./outtext/lsFiles.conf"
	ls $DATAPATH | grep $LSFILESID | tail -1 >> "./outtext/lsFiles.conf"
else
	STRWARNINGMANUAL="WARNING: manual selection mode is on, run $NRUN selected!"
	echo $STRWARNINGMANUAL
	echo $STRWARNINGMANUAL >> "./outtext/infoCycle.conf"
	LSFILESIDLAST=$(ls $DATAPATH | tail -1 | cut -c$RUNSTRL-$RUNSTRL)
	ls $DATAPATH | grep $NRUN | tail -$NFILES > "./outtext/lsFiles.conf"
	ls $DATAPATH | grep $LSFILESIDLAST | tail -1 >> "./outtext/lsFiles.conf"
fi
# comment # uncomment this (& other "<<comment" above) & create a lsFiles.conf manually to use the latter

# also write the latter to text file
rm -f ./outtext/lsFilesAdv.conf ; touch ./outtext/lsFilesAdv.conf
I=1
NLINESTOT=0
while read FILE0 ; do
	NLINES0=$(wc -l < "$DATAPATH/$FILE0")
	if [ $I -eq $(wc -l "./outtext/lsFiles.conf" | awk '{ print $1 }') ] ; then
		echo "-----" >> "./outtext/lsFilesAdv.conf"
		echo "compared to the latest-spill file:" >> "./outtext/lsFilesAdv.conf"
		NLINESLAST=$NLINES0
	fi
	printf "%0*d) %s: %5d lines; modified @ %s\n" 3 $I $FILE0 $NLINES0 "$(date -r $DATAPATH$FILE0 +'%Y %b %d, %X')" >> "./outtext/lsFilesAdv.conf"
	let "I+=1"
	let "NLINESTOT+=NLINES0"
done < "./outtext/lsFiles.conf"
if [ $NLINESLAST -eq 0 ] ; then
	STRWARNINGLAST0="WARNING: 0 lines in last spill --> no last-spill histograms available!"
	echo $STRWARNINGLAST0
	echo $STRWARNINGLAST0 >> "./outtext/infoCycle.conf"
fi
if [ $NLINESTOT -eq 0 ] ; then
	STRWARNINGTOT0="WARNING: 0 lines in total --> no analysis at all is executed!"
	echo $STRWARNINGTOT0
	echo $STRWARNINGTOT0 >> "./outtext/infoCycle.conf"
else
	echo "total nr. of lines read: $NLINESTOT" >> "./outtext/infoCycle.conf"
fi
	
# remove outdated output files
rm -f ./outhist/*
rm -f ./outtree/*
rm -f ./outpdf/*
rm -f ./outfigs/*
touch ./outhist/.keep
touch ./outtree/.keep
touch ./outpdf/.keep
touch ./outfigs/.keep

# recompile analysis software to avoid system mismatch after syncing
if [ $EXCLUDECOMPILER -eq 0 ] ; then
	rm -f krysBTOM_ana.out
	echo "recompiling analysis software (krysBTOM_ana.cpp)..."
	g++ krysBTOM_ana.cpp -o krysBTOM_ana.out `root-config --cflags --libs`
else
	echo "recompilation of the analysis software (krysBTOM_ana.cpp) not performed"
fi

# data analysis -> output histogram (& maybe tree files)
if [ $NLINESTOT -ne 0 ] ; then
	./krysBTOM_ana.out
fi

# talk to server (copy output forth, sync website back)
source talkToServer.sh $TOSERVER $SYNCBACKSERVER $SERVERPATH

echo "done! now waiting for new data..."
echo "-----"