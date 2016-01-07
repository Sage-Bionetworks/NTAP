#!/bin/bash

for file in `ls /scratch/NTAP/rnaseq/hhao_123771/FASTQ/*.gz`
do
cmd="fastqc $file -o ./"
echo $cmd
$cmd &
done
