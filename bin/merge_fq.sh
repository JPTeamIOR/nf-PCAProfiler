#!/bin/bash

f1=$1
f2=$2
outf=$3

zcat $f1 | awk -v F2=${f2} '
    { print;
    if ( NR %4 == 0 ) {
        for ( i=0 ; i < 4 ; i++ ) {
            "zcat "F2 | getline sline ;
            print sline
        }
    }
    }
    END { close("zcat "F2) } ' > $outf
