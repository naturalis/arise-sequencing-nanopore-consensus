#!/bin/bash

# calculate number of reads and median of read lengths for all fastq files in folder
# medians seem to increase with a decreasing number of reads

# requirement: 'nanopore' environment (conda)
# ----
# conda activate nanopore
# conda install --yes -c conda-forge -c bioconda medaka==0.11.5 openblas==0.3.3 spoa racon minimap2
# pip install NGSpeciesID
# conda install -c bioconda prinseq seqtk cutadapt
# ----

# usage: ont_consensus.sh $1 $2 $3
# ----
#   $1 = name_code_index.txt (index [tab] replacemnent header)
#   e.g.    CCTCCAACCGCTG Andrena_spinigera___RMNH.INS.1092535___CCTCCAACCGCTG
#   $2 = forward primer (without the index)
#   $3 = revcomp of reverse primer (without the revcomp_index)
#   e.g. index-forward_primer-------amplicon-------revcomp_rev_primer-revcomp_index
# ----

# eval statement prevents conda activate from erroring out
eval "$(conda shell.bash hook)"
conda activate nanopore

# prevent files from being overwritten
[ -d out_length_filtered ] && { printf "out_length_filterd EXISTS !!! \nplease remove to continue\n"; exit 1; }
mkdir -p out_length_filtered
[ -f tmp.txt ] && { printf "tmp.txt !!! \nplease remove to continue\n"; exit 1; }
[ -f tmp2.txt ] && { printf "tmp2.txt EXISTS !!! \nplease remove to continue\n"; exit 1; }
[ -d out_length_filtered/out_sum ] && { printf "out_length_filtered/out_sum EXISTS !!! \nplease remove to continue\n"; exit 1; }
mkdir -p out_length_filtered/out_sum

# create associative array based on $1
declare -A index_array
while read index name
do
    index_array[$index]=$name
done < $1

# calculate number of reads and median of lengths
printf "File\t#_reads\tmedian(length)\tmin\tmax\n" >> tmp.txt
for i in $(ls -1 *.fastq)
do
    seq_count=$(wc -l < $i | sed 's/$/\/4/' | bc)
    seq_median=$(awk 'NR % 4 == 2' $i | awk '{print length($0)}' | sort -n |
     awk '{a[i++]=$0;s+=$0}END{print int((a[int(i/2)]+a[int((i-1)/2)])/2)}')
    printf "$i\t$seq_count\t$seq_median\n" >> tmp.txt
done

# calculate min_len and max_len (set to 50 nt below and above median)
printf "\n--- less_than_100_reads ---\n"
tail -n +2 tmp.txt | awk '$2 < 100' | sort -nk2 | cat -n  | awk 'BEGIN{print "_","file","#_reads","median(length)"}1'| column -t
printf "\n--- 100_or_more_reads ---\n"
tail -n +2 tmp.txt | awk '$2 >= 100' | awk '{print $0,($3-25),($3+25)}' | sort -nk2 | cat -n  | awk 'BEGIN{print "_","file","#_reads","median(length)","min","max"}1'| column -t
echo ""
tail -n +2 tmp.txt | awk '$2 >= 100' | awk '{print $0,($3-25),($3+25)}' | sort -nk2 > tmp2.txt

# use prinseq to filter on length
cat tmp2.txt | while read line
do
    index=$(echo "$line" | awk '{print $1}')
    min=$(echo "$line" | awk '{print $(NF-1)}')
    max=$(echo "$line" | awk '{print $NF}')
    printf "Prinseq filtering $index in length range $min - $max ..."
    prinseq-lite.pl -fastq $index -min_len $min -max_len $max -out_bad null -out_format 3 >> /dev/null 2>&1
    printf "  \t\tdone\n"
done

# cleanup temporary files
rm tmp.txt tmp2.txt

# move output to directory and rename files
mv *prinseq_good* out_length_filtered/
cd out_length_filtered
for i in $(ls -1 *.fastq)
do
    mv "$i" `echo "$i" | sed -e 's/\(.\{13\}\).*.\(.\{5\}\)/\1.\2/'`
done

# run NGSpeciesID on each fastq
echo ""
for i in $(ls -1 *.fastq)
do
    printf "NGSpeciesID processing $i ..."
    NGSpeciesID --ont --fastq $i --outfolder "$i".out --consensus --medaka >> /dev/null 2>&1
    printf "  \t\t\t\tdone\n"
done

find . -type f -name "*.fasta" > out_sum/tmp_01.txt                         # tmp_01.txt = list of paths to all fasta files
cat out_sum/tmp_01.txt | egrep -v "medaka" > out_sum/tmp_02.txt             # tmp_02.txt = list of paths to correct fasta files
cat out_sum/tmp_02.txt | cut -c 3-15 > out_sum/tmp_03.txt                   # tmp_03.txt = barcodes
cat out_sum/tmp_02.txt | awk -F"_" '{print "_"$3}' > out_sum/tmp_04.txt     # tmp_04.txt = extensions
paste -d "" out_sum/tmp_03.txt out_sum/tmp_04.txt > out_sum/tmp_05.txt      # tmp_05.txt = new names
awk '{print "out_sum/"$0}' out_sum/tmp_05.txt > out_sum/tmp_06.txt          # tmp_06.txt = add out location
paste -d " " out_sum/tmp_02.txt out_sum/tmp_06.txt > out_sum/tmp_07.txt     # tmp_07.txt = location destination
awk '{print "cp "$0}' out_sum/tmp_07.txt > out_sum/tmp_08.txt               # tmp_08.txt = copy statement

# echo ""
# cat out_sum/tmp_08.txt

# copy files to out
cat out_sum/tmp_08.txt | while read -r i
do
    $i
done

# remove the tmp files
rm out_sum/*.txt

# remove the output folders from ngspid.sh
# note: turned off by default, because ngspid.sh takes computation time
# and this script collects only part of the output
# ---
# rm *.out

mv out_sum ..
cd ../out_sum

# index2header
# create replacement loop
for fasta in $(ls -1 | grep "fasta$")
do
    fasta_index=$(echo $fasta | awk -F"_" '{print $1}')
    fasta_extension=$(echo $fasta | awk -F"_" '{print "_"$2}')
    consensus_number=$(echo $fasta_extension | sed 's/\.fasta//g')
# replace the fasta consensus headers
    sed 's/^>consensus.*reads/>consensus_reads/g' "$fasta" |
    sed "s/consensus/${index_array[$fasta_index]}$consensus_number/g" > "$fasta".tmp && mv "$fasta".tmp "$fasta"
# replace the fasta filename
    mv $fasta ${index_array[$fasta_index]}$fasta_extension
done

cat *.fasta > ../out_untrimmed.fasta
cd ..

echo ""
# trim sequences
cutadapt -a "$2"..."$3" --rc -o out_trimmed.fasta out_untrimmed.fasta

conda deactivate
