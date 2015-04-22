#!/bin/bash 

base_dir=`pwd`

for (( FF=0 ; FF < 107 ; ++FF )) ; do
    mkdir ${FF}
    cd ${FF}
    for (( KK=21 ; KK < 152 ; KK+=10 )) ; do
        mkdir ${KK}
        cd ${KK}
        ~/wga/src/BroadCRD/other_assemblers/velvet/fosmid.sh ~/wga/velvet/software161 ${FF} ${KK} 400 50 > run_log
        cd ..
    done
    cd ..
done


