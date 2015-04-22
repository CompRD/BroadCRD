#!/bin/bash
qdir=/wga/scr4/qgrind
if (( $# > 0 ))
then
    color=$1
else
    color=""
fi

while :
do
    cols=`stty -a | head -1  | awk -F\; '{print $3}' | awk '{print $2}'`
#    let cols=cols-1
    clear
    cat $qdir/crd?_*.status $qdir/crd??_*.status | cut -b1-$cols | (
    while read line
    do
        stat=`echo "$line" | awk '{print $4}'`
        if [[ "$color" ]]
        then
            RED='\033[31m'
            GREEN='\033[32m'
            RESET='\033[39m'
            if [[ "$stat" == "RUN" ]]
            then
                echo -en "$GREEN"
            else
                echo -en "$RESET"
            fi
            echo $line
        else
            echo $line
        fi
    done
    if [[ "$color" ]]
    then
        echo -en "$RESET"
    fi
    )
    sleep 5
done
