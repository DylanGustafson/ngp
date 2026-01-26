#!/bin/bash

start=$1
stop=$2
chunk=$3
class=$4
block=$5
outdir=$6
logfile=$7

touch $logfile
echo >> $logfile
echo $(date "+%F %T") >> $logfile
echo "Running from $start to $stop, with $class classes. Thread chunk size: $chunk. Writing every $block." >> $logfile

tic=$(date +%s)

for ((i=$start; i<$stop; i+=$block)); do

    if [ $block -lt "1000000000" ]; then
        from=$(printf "%014d" $i)
        to=$(printf "%014d" $(($i+$block)))
    elif [ $block -lt "1000000000000" ]; then
        from="$(printf "%06d" ${i:0:$((${#i} - 9))})e9"
        to_long=$(($i+$block))
        to="$(printf "%06d" ${to_long:0:$((${#to_long} - 9))})e9"
    else
        from="$(printf "%06d" ${i:0:$((${#i} - 12))})e12"
        to_long=$(($i+$block))
        to="$(printf "%06d" ${to_long:0:$((${#to_long} - 12))})e12"
    fi

    outfile="$outdir/$from-to-$to.csv"
    ./ngp-bin c $i $block $chunk $class > $outfile 2>> $logfile

done

toc=$(date +%s)

echo
echo $(date "+%F %T") >> $logfile
echo "Finished! Total time: $(($toc-$tic)) seconds" >> $logfile
