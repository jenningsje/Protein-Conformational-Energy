#!/bin/bash

atoms=( ["C"]=1 ["N"]=2 ["O"]=3 ["ZN"]=4)

acids=( ["ALA"]=1, ["ARG"]=2, ["ASN"]=3, ["ASP"]=4, ["CYS"]=5, ["GLU"]=6, ["GLN"]=7, ["GLY"]=8, ["HIS"]=9, ["LIE"]=10, ["LEU"]=11, ["LYS"]=12, ["MET"]=13, ["PHE"]=14, ["PHO"]=15, ["SER"]=16, ["THR"]=17, ["TRP"]=18, ["TYR"]=19, ["VAL"]=20, ["HOH"]=21)

for i in 'grep -l " ${acids[$i]}" protein_coordinate_database.csv'
do
    sed -i "s/${i}/${acids[$i]}/g" $i;
done
