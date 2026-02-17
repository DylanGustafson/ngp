#!/bin/bash

start=$1
stop=$2
chunk=$3
block=$4
outdir=$5
logfile=$6

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

touch $logfile
echo >> $logfile
echo $(date "+%F %T") >> $logfile
echo "Running from $start to $stop, with 1000 q/N bins. Thread chunk size: $chunk. Writing every $block." >> $logfile

tic=$(date +%s)

for ((i=$start; i<$stop; i+=$block)); do

    if [ $block -lt "1000000000" ]; then
        from=$(printf "%014d" $i)
        to=$(printf "%014d" $(($i+$block)))
    elif [ $block -lt "1000000000000" ]; then
        if [ $i -eq 0 ]; then
            from="000000e9"
        else
            from="$(printf "%06d" ${i:0:$((${#i} - 9))})e9"
        fi

        to_long=$(($i+$block))
        to="$(printf "%06d" ${to_long:0:$((${#to_long} - 9))})e9"
    else
        if [ $i -eq 0 ]; then
            from="0000000e12"
        else
            from="$(printf "%07d" ${i:0:$((${#i} - 12))})e12"
        fi

        to_long=$(($i+$block))
        to="$(printf "%07d" ${to_long:0:$((${#to_long} - 12))})e12"
    fi

    outfile="$outdir/$from-to-$to.csv"
    echo "./ngp-bin c $i $block $chunk" >> $logfile
    ./ngp-bin c $i $block $chunk > $outfile 2>> $logfile

done

toc=$(date +%s)

echo
echo $(date "+%F %T") >> $logfile
echo "Finished! Total time: $(($toc-$tic)) seconds" >> $logfile
