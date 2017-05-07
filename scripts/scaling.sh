#!/bin/bash

#RUN_PATH=./parallel/Code/exec
#COMMAND="$RUN_PATH/parallel.exe -i $RUN_PATH/input.txt"
#RUN_PATH=./serial/Code/exec
#COMMAND="$RUN_PATH/groupF.exe -i $RUN_PATH/input.txt"
RUN_PATH=./eigen/Code/exec
COMMAND="$RUN_PATH/parallel.exe -i $RUN_PATH/input.txt"
OUTFILE=eigen_scaling.txt

for nprocs in 2 4 6 8 10 12 14 16
do    
    echo "Procs, ${nprocs}" >> $OUTFILE
    $COMMAND -n $nprocs >> $OUTFILE
done
