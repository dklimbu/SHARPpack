#!/bin/bash -I

root2bin=~/Desktop/MyDocs/SHARPpack/bin
exe=${root2bin}/sharp.x

echo JOB STARTED AT:
 date

k1=5
k2=35
dk=1

for (( i=1; $i<=$(bc<<<"($k2-$k1)/$dk+1"); i++ )); 
do
  echo $(bc<<<"$k1+($i-1)*$dk") > pkin

  ${exe} < pkin  

done

echo JOB FINISHED AT:
 date
echo ALL DONE!!!

