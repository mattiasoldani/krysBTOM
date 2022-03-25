#!/bin/bash

# delete all data from server (in order to prepare it for a new beamtest)
# one argument needed: server location (either remote or mounted locally)

SERVERPATH="$1"

echo "cleaning server up from all latest output..."
rm -f $SERVERPATH/outdata/outtext/*
rm -f $SERVERPATH/outdata/outpdf/*
rm -f $SERVERPATH/outdata/outfigs/*
rm -f $SERVERPATH/outdata/outhist/*

echo "also, creating empty/dummy info files to be displayed in the website..."
touch $SERVERPATH/outdata/outtext/lsFilesAdv.conf
echo "OM empty; ready for the next beamtest!" > $SERVERPATH/outdata/outtext/infoCycle.conf

echo "done"
