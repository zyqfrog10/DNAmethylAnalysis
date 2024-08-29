#!/bin/bash
# compute methylation rate from pUC19 methylated CpG

dir=$1 # project directory

echo -e 'sample_name\tCpG_methylation_rate' > ${dir}/CpG_methylation_rate.txt

for file in `find ${dir} -name "bismark_summary_report.txt"`
do 
    sampleName=$(basename $(dirname $file))
    rate=`sed '1d' ${file} |awk '{printf "%.4f", $10/($10+$11)}'`
    echo -e "${sampleName}\t${rate}" >> ${dir}/CpG_methylation_rate.txt
done
