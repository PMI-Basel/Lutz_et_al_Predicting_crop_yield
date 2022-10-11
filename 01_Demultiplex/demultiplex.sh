#!/bin/bash

#SBATCH --qos=1day
#SBATCH --time=24:00:00
#SBATCH --mem=20g
#SBATCH --output=demultiplex.out
#SBATCH --error=demultiplex.error
#SBATCH --job-name=demultiplex
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules
module load foss/2018b #interpreters
module load SMRT-Link/9.0.0.92188-cli-tools-only
module load SAMtools/1.9-foss-2018b
module load picard/2.9.2
module load cutadapt/2.10-foss-2018b-Python-3.6.6

#raw data must be called: name.subreads.bam (name can't have a ".")
all_raw_data=$(ls ../00_raw_data/)

# running time notification
echo "start script"

for raw_data in $all_raw_data; do

    # get name
    name=$(echo ${raw_data} | cut -f1 -d ".")

    #running time notification
    echo ${name} "started"

    #convert raw data to circular consensus sequences (package SMRT)
    ccs ../00_raw_data/${raw_data} ${name}.ccs.bam --min-passes 5

    #running time notification
    echo ${raw_data} "converted to ccs"

    #convert ccs fastq to bam
    #useful if you start directly with ccs file in fastq format
    # java -jar $EBROOTPICARD/picard.jar FastqToSam \
    #     F1=../00_raw_data/${name}.ccs.fastq \
    #     O=${name}.ccs.bam \
    #     SM=sample001 \
    #     RG=rg0013

    #demultiplex the ccs file (package SMRT)
    #barcodes must be in one file (forward and reverse) in a file called barcode.fasta
    #barcode.fasta must be in a folder called ${name}
    #forward primers must start with "F"
    #reverse primers must start with "R"
    #no IUPAC bases allowed
    rm -r bam 2> /dev/null
    mkdir bam 2> /dev/null
    lima ${name}.ccs.bam ${name}/barcode.fasta bam/${name}.xml --different --split-bam-named --peek-guess --ccs

    #convert bam file to fastq file
    mkdir out 2> /dev/null
    rm -r out/${name} 2> /dev/null
    mkdir out/${name}
    cd bam
    ls *.bam > bam_names.txt
    while read infile; do
      outfile=$(echo ${infile} | cut -f3 -d ".") #adapt this part
      #samtools bam2fq ${infile} > ../out/${name}/${outfile}.uncut.fastq
      samtools bam2fq ${infile} > ../out/${name}/${outfile}.fastq
    done < bam_names.txt
    rm bam_names.txt
    cd ..
    rm -r bam

    #delete F-F and R-R combinations
    #forward primers must start with "F"
    #reverse primers must start with "R"
    cd out/${name}
    rm *--F* R* #F--F or R--R combinations
    rename "-" "" *
    cd ../..

    #running time notification
    echo ${raw_data} "demultiplexed"

    #remove primers
    #primers must be in one file (forward and reverse) in a file called primer.fasta
    #primer.fasta must be in a folder called ${name}
    #the first sequence is the forward primer (seq on 2nd line in file)
    #the second sequence is the reverse primer (seq on 4th line in file)
    #IUPAC bases allowed

    # cd out
    # ls ${name} > filenames.txt
    # while read infile; do
    #   outname=$(echo ${infile} | cut -f1 -d ".") #filename without .uncut.fastq
    #   forward=$(sed '2q;d' ../${name}/primer.fasta)
    #   reverse=$(sed '4q;d' ../${name}/primer.fasta)
    #   cutadapt -g 'X'${forward} -o temp.fastq ${name}/${infile} #"X" = non-internal adapters
    #   cutadapt -g 'X'${reverse} -o ${name}/${outname}.fastq temp.fastq #"X" = non-internal adapters
    #   rm temp.fastq
    #
    # done < filenames.txt
    # rm filenames.txt
    #
    # mkdir ${name}/uncut
    # mv ${name}/*uncut.fastq ${name}/uncut
    # cd ..

done

# running time notification
echo "end script"
echo "the following data has been processed:"
echo ${all_raw_data}
