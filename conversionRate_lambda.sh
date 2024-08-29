#!/bin/bash
# compute conversion rate from lambda unmethylated Cs

dir=$1 # project directory

echo -e 'sample_name\tC_conversion_rate' > ${dir}/C_conversion_rate.txt

for file in `find ${dir} -name "bismark_summary_report.txt"`
do 
    sampleName=$(basename $(dirname $file))
    rate=`sed '1d' ${file} |awk '{printf "%.4f", ($11+$13+$15)/$9}'`
    echo -e "${sampleName}\t${rate}" >> ${dir}/C_conversion_rate.txt
done
