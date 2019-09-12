#!/bin/bash 

FILE=$1
while read LINE; do
    STR=`samweb -e sbnd locate-file $LINE`
    STR=${STR/(*)}
    STR=${STR/enstore:}
    STR="$STR/"
    STR=$STR$LINE
    echo $STR >> "complete_path.txt"
done < $FILE
