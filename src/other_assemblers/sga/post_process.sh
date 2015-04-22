#! /bin/bash


# make dot file from assembly
zcat primary-graph.asqg.gz | /wga/scr4/assemblers/sga_standalone/sga/src/bin/sga-asqg2dot.pl > primary-graph.dot


# try to map variants back to contigs
MakeLookupTable SOURCE=primary-contigs.fa

QueryLookupTable SEQS=primary-variants.fa  K=12 L=primary-contigs.fa.lookup \
                 TARGET_NAMING="from_record" QUERY_NAMING="from_record"\
                 VISUAL=False PARSEABLE=True SMITH_WAT=True NH=True \
                 | grep QUERY \
                 | awk '{for(i=1;i<=NF;++i) if(i!=1 && i!=8 && i!=9 ) printf("%s\t",$i);printf("\n")}'  \
                 | sed -e 's@\(.*variant-[0-9]\+\)/\(.*$\)@\1\t\2@' \
                 > tmp.tmp

python ~/wga/src/BroadCRD/other_assemblers/sga/variants2contigs.py tmp.tmp > variants_contigs.query 
rm tmp.tmp
