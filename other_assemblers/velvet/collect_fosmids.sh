#!/bin/bash 

base_dir=`pwd`
fosmids_fasta=${base_dir}/fosmids.fasta

rm -rf ${fosmids_fasta}
#for (( NN=0 ; NN < $3 ; ++NN )) ; do
  for (( FF=$1 ; FF < $2 ; ++FF )) ; do
  #    mkdir ${FF}
      cd ${FF}


      maxPen=100000000
      bestK=""
      for (( KK=21 ; KK < 152 ; KK+=10 )) ; do
  #        mkdir ${KK}
          cd ${KK}
          python ~/wga/src/BroadCRD/other_assemblers/velvet/contigs2scaffold.py output/contigs.noN.fa ref F${FF}_K${KK} > ${base_dir}/${FF}_${KK}.fasta
          grep '>' ${base_dir}/${FF}_${KK}.fasta
          Pen=`grep '>' ${base_dir}/${FF}_${KK}.fasta | sed -e 's/_/ /g' | awk '{print $NF}' `

          NBases=`grep '>' ${base_dir}/${FF}_${KK}.fasta | sed -e 's/_/ /g' | awk '{print $3}' `
          RefLength=`grep '>' ${base_dir}/${FF}_${KK}.fasta | sed -e 's/_/ /g' | awk '{print $4}' `
          if [[ ("${Pen}" -lt "${maxPen}")  ]] ;then
              maxPen=${Pen}
              bestK=${KK}
          fi
          cd ..
      done

      echo "best K is ${bestK}"

      cat ${base_dir}/${FF}_${bestK}.fasta >> ${fosmids_fasta}
      cd ..
  done
#done

