#!/bin/bash

if [ ! $# -eq 3 ]; then
  echo "merge 2 vcf files with possible overlaping regions" 1>& 2
  echo "Usage: `basename $0` file1 file2 output" 1>& 2
  exit 1
fi

f1=$1
f2=$2
f3=$3

if [ ! -f $f1 ];then
  echo $f1 "doesn't exist" 1>& 2
  exit 1
fi
if [ ! -f $f2 ];then
  echo $f2 "doesn't exist" 1>& 2
  exit 1
fi
if [ -f $f3 ];then
  echo $f3 "already exists" 1>& 2
  exit 1
fi

exec_list=("tabix" "bgzip" "vcf-isec" "vcf-concat" "vcf-sort")
for e in "${exec_list[@]}";do
  command -v $e > /dev/null 2>&1 || { echo >&2 "$e cannnot be found"; exit 1;}
done

z1=${f1}.gz
z2=${f2}.gz

bgzip -c $f1 > $z1
tabix $z1
bgzip -c $f2 > $z2
tabix $z2

vcf-isec $z1 $z2 | bgzip -c > i.vcf.gz
vcf-isec -c $z1 $z2 | bgzip -c > i1c.vcf.gz
vcf-isec -c $z2 $z1 | bgzip -c > i2c.vcf.gz

vcf-concat i.vcf.gz i2c.vcf.gz i1c.vcf.gz | vcf-sort >  $f3

bgzip -d -c i.vcf.gz > i.vcf

rm $z1 $z2 i.vcf.gz i1c.vcf.gz i2c.vcf.gz
