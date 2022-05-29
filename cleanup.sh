#!/bin/bash

# delete all data from server (in order to prepare it for a new beamtest)
# one argument needed: server location (either remote or mounted locally)

SERVERPATH="$1"

echo "cleaning server up from all latest output..."
rm -f $SERVERPATH/outtext/*
rm -f $SERVERPATH/outpdf/*
rm -f $SERVERPATH/outfigs/*
rm -f $SERVERPATH/outhist/*
rm -f $SERVERPATH/outtree/*

echo "also, creating empty/dummy info files to be displayed in the website..."
touch $SERVERPATH/outtext/lsFilesAdv.conf
echo "OM empty; ready for the next beamtest!" > $SERVERPATH/outtext/infoCycle.conf

echo "done"
