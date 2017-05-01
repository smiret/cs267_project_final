#!/bin/bash

#RUN_PATH=./parallel/Code/exec
#COMMAND="$RUN_PATH/parallel.exe -i $RUN_PATH/input.txt"
RUN_PATH=./serial/Code/exec
COMMAND="$RUN_PATH/groupF.exe $RUN_PATH/input.txt"
OUTFILE=serial_time.txt

for size in 8 10 12 14 16 18 20
do
    sed -i "2s/.*/${size}/" $RUN_PATH/input.txt
    echo "Size, ${size}" >> $OUTFILE
    $COMMAND >> $OUTFILE
done