#!/bin/bash

# after cycle, sync website from web server to local folder
# one argument needed: server location (either remote or mounted locally)

echo "syncing website back from server"
rsync -arzu -delete $1 "./website"
touch ./website/.keep