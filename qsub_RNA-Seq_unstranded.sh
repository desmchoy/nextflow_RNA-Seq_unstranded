#!/bin/bash

task_name=RNA-Seq_unstranded

if [ ! -d "results_FastQC" ] ; then
	mkdir results_FastQC
fi




while read line ; do
	
	qsub -v sample=$line -o "$(pwd -P)"/log."${task_name}"."${line}".output -e "$(pwd -P)"/log."${task_name}"."${line}".error sge_"${task_name}".batch

done < "sample.list"
