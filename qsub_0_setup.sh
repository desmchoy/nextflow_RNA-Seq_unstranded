#!/bin/bash


ls ./RAW/*.gz | sed 's/\.\/RAW\///g;s/_.*_001//g;s/\.fastq\.gz//g' | sort | uniq > sample.list
