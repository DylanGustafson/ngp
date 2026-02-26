#!/bin/bash
start=$1
stop=$2
chunk=$3
block=$4
logfile=$5

# Get initial g(p) maximum from end of gp-max file
gp_max_file="output/gp-max.csv"
if [ ! -f $gp_max_file ]; then
    gp_max="3"
else
    gp_max=$(tail -n 1 $gp_max_file | cut -d ',' -f 2)
fi
gp_fallback=$gp_max

# Write header to log file
touch $logfile
echo >> $logfile
echo $(date "+%F %T") >> $logfile
echo "Running from $start to $stop, with 1000 q/N bins."
echo "Thread chunk size: $chunk. Writing every $block. Initial gp_max: $gp_max" >> $logfile

tic=$(date +%s)

# Loop through blocks and execute ngp-bin
for ((i=$start; i<$stop; i+=$block)); do
    cmd="./ngp-bin c g=$gp_max $i $block $chunk"
    echo "Executing: $cmd" >> $logfile
    gp_max=$($cmd 2>> $logfile)

    if [ $? -ne 0 ]; then
        gp_max=$gp_fallback
    fi
done

# Print total time and exit
toc=$(date +%s)
echo "Finished! Total time: $(($toc-$tic)) seconds" >> $logfile
echo $(date "+%F %T") >> $logfile
