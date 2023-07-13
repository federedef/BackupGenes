#!/usr/bin/env bash

#cut -f 5 Big_Papi | sort
grep -v -w 'All 6\|no *' Big_Papi > probando
grep -v -w 'All 6\|no *' Big_Papi |  awk '{FS="\t";OFS="\t"}{if ( $6 >= 0.05 ) print $1,$2}' | sort -k1 | uniq -c | grep -w "6" | sed "s/6//1"
#grep -w 'All 6\|no *' Big_Papi | awk '{FS="\t";OFS="\t"}{if ( $6 ==1) print $1,$2}'
