#!/bin/bash

arq=$1
sed '11,24!d' Fe.nml >> ${arq} 
line_ini=$(grep -in "center_band_s_up" ${arq} | awk -F: '{print $1}')
line_fin=$((${line_ini}+11)); count=1 ;
for i in $(seq ${line_ini} ${line_fin}) ; do
v=$(sed "${i}"'!d' ${arq} | awk '{print $3}')
if [ $(echo "${count} <= 3" | bc) -eq 1 ] ; then
value=$(sed "${count}"'!d' uppar | awk '{print $1}') ;
elif [ $(echo "${count} <= 6" | bc) -eq 1 ] ; then
k=$((${count}-3)) ;
value=$(sed "${k}"'!d' dwpar | awk '{print $1}') ;
elif [ $(echo "${count} <= 9" | bc) -eq 1 ] ; then
k=$((${count}-6)) ;
value=$(sed "${k}"'!d' uppar | awk '{print $2}') ;
else
k=$((${count}-9)) ;
value=$(sed "${k}"'!d' dwpar | awk '{print $2}') ; fi
sed -i "${i}"'s/'"${v}"'/'"${value}"'/g' ${arq}
count=$((${count}+1))
done
