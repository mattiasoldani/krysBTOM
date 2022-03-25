#!/bin/bash

# copy all output file to web server & (if requested) sync website back
# 3 arguments needed: boolean for enable copy to server, boolean for enable sync from server, server location (either remote or mounted locally)
TOSERVER="$1"
SYNCBACKSERVER="$2"
SERVERPATH="$3"
 
echo "talking to server @ $SERVERPATH..." 
if [ $TOSERVER -eq 1 ] ; then
	if [ -d $SERVERPATH ] && [ ! -z "$(ls $SERVERPATH)" ] ; then
		echo "copying all files in outtext, outpdf, outfigs, outhist files to server"
		# note: tree files are never sent, regardless of their generation
		rm -f $SERVERPATH/outdata/outtext/* ; cp -r "outtext" "$SERVERPATH/outdata"
		rm -f $SERVERPATH/outdata/outpdf/* ; cp -r "outpdf" "$SERVERPATH/outdata"
		rm -f $SERVERPATH/outdata/outfigs/* ; cp -r "outfigs" "$SERVERPATH/outdata"
		rm -f $SERVERPATH/outdata/outhist/* ; cp -r "outhist" "$SERVERPATH/outdata"
		if [ $SYNCBACKSERVER -eq 1 ] ; then
			source syncBackServer.sh $SERVERPATH
		else
			echo "chose not to sync website back from server"
		fi
	else
		echo "server not found (directory inexistent or empty!)"
	fi
else
	echo "no, chose not to talk to server"
fi
