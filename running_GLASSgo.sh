#! /bin/bash

%% novel miRNAs from NCI-60
n=`grep ">" /home/mlee/Borellia_sRNA/GLASSgo/miRNA/novel_miRNA_NCI60.fasta | wc -l`

for i in $(eval echo "{1..$n}")
do
echo $i
line_start=`expr $i \* 2 - 1`
line_end=`expr $line_start + 1`
sed -n "${line_start},${line_end}p" /home/mlee/Borellia_sRNA/GLASSgo/miRNA/novel_miRNA_NCI60.fasta > query
fname="novel_miRNA_NCI60_"
fname+=`echo $i`
python3.6 /home/mlee/Borellia_sRNA/GLASSgo/GLASSgo_1_5_2/GLASSgo_1_5_2.py -d /home/mlee/blastdb/nt/nt -i ./query -t 2  > `echo $fname`
done
