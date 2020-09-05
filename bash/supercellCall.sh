#!/bin/bash

# i="$HOME/proj/proMin/minerals/supercell/actinolite/1981/Actinolite_1981.cif"
# s='2x2x2'
i=$1
s=$2
o=$3
# e=$4
# echo $i
# echo $s




# grep -m 1 return first match 
err=$(timeout 2 supercell -i $i -s $s -c no -m -d | grep -m 1 '^ERROR:')

# -n returns True not empty / -z returns True if empty
if [[ -n $err ]];
then 
    comb=`echo $err`
    printf "error,$s,%s,%s\n" "$i" "$comb" #>> $e
    exit()
fi


# grep -m 1 return first match 
warn=$(timeout 2 supercell -i $i -s $s -c no -m -d | grep -m 1 '^WARN')

# -n returns True not empty / -z returns True if empty
if [[ -n $warn ]];
then 
    comb=`echo $warn | rev | cut -d' ' -f1 | rev`
    printf "No,$s,%s,%s\n" "$i" "$comb" >> $o
    exit()
fi

# grep "The total number of combinations is 1"
comb=$(timeout 2 supercell -i $i -s $s -c no -m -d | grep -m 1 '^The total number of combinations is')

if [[ -n $comb ]];
then 
    comb=`echo $comb | rev | cut -d' ' -f1 | rev`
    printf "Yes,$s,%s,%s\n" "$i" "$comb" >> $o
    exit()
fi
