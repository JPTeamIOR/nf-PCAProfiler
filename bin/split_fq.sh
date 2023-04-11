#!/bin/bash

infile=$1
f1=$2
f2=$3


awk -v F2=${f2} '
    { print ;
    if ( NR %4 == 0 ) {
        for ( i=0 ; i < 4 ; i++ ) {
            getline
            print >> F2
        }
    }
    }
' $infile > $f1
