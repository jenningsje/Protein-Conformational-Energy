#!/bin/bash

echo This bash script will extract all of the cif files from the protein databank into the current working directory
echo afterwards it will compute the conformational entropies pertaining to each protein within that archive.
echo shall we proceed with the program? (y/n)

read response

yes_response=(Y y Yes yes)
no_response=(N n No no)

count = 0

while [[count = 0]]

    if [response in yes_response]
    then

        ./batch_download.sh -f list_file.txt -c
        chmod +x batch_download.sh

        find . 

        python Protein-Cnfrm-Energy.py pwd

        echo You will find the database containing the equations
        echo pertaining to every protein within the protein databank
        echo within the following directory:
        pwd

        count = 1

    elif [response in no_response]
    then

        echo end program

        count = 1

    else
    then
        echo please give me a yes or no response