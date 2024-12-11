#!/bin/bash


##to run any chunk of this script, run "bash bashName.sh commandname", (commandname = demux for first) at the command line

## Make sure to read all comments and make the necessary edits to your own saved copies of each piece of code used

##Download python, biopython, bwa, samtools - just google and follow instructions to download to your personal computer

##your reads will already be demultiplexed by illumina barcodes (i12/i6/i19/i5), but we still need to demultiplex by custom fwd primers. 
##This script will take illumina barcode fasta (for example i5) and spit out up to six demultiplexed output fastas by fwd primer
## change the argument for each fwd primer (e.g. --f1 f1) to whatever sample it is (e.g. --f1 Newman_D5)
 if [ $1 == "demux" ]; then
 	mkdir demultiplexed
 	touch demux_ultimeta.txt
 	for i in ./i5*.fastq
 	do python3 demultiplexerV2.py --input $i --f1 JW5507_P0 --f2 i5f2 --f3 i5f3 --f4 i5f4 --f5 i5f5 --f6 JW5513_P8;
 	done
 	for i in ./i6*.fastq
 	do python3 demultiplexerV2.py --input $i --f1 JW5507_P8STREP --f2 JW5508_P0 --f3 i6f3 --f4 i6f4 --f5 i6f5 --f6 i6f6;
 	done
 	for i in ./i7*.fastq
 	do python3 demultiplexerV2.py --input $i --f1 i7f1 --f2 JW5508_P8STREP --f3 JW5513_P0 --f4 i7f4 --f5 i7f5 --f6 i7f6;
 	done
 	for i in ./i12*.fastq
	do python3 demultiplexerV2.py --input $i --f1 i12f1 --f2 i12f2 --f3 JW5513_P8STREP --f4 JW5507_P8 --f5 i12f5 --f6 i12f6;
	done
	for i in ./i19*.fastq
	do python3 demultiplexerV2.py --input $i --f1 i19f1 --f2 i19f2 --f3 i19f3 --f4 i19f4 --f5 JW5508_P8 --f6 i19f6;
	done
	cat ./demultiplexed/*meta.txt >> ./demux_ultimeta.txt
	
fi

#if [ $1 == "convert" ]; then
#	mkdir converted
#	for i in ./*.fastq; 
#	do python fastqToFasta.py --input $i;
#	done
#fi

##This script will loop through all of your demultiplexed fastas and pull out spacer sequences
# Every file in the directory ./demultiplexed will be parsed
if [ $1 == "parse" ]; then
	mkdir parsed
	touch parse_ultimeta.txt
	for i in ./demultiplexed/*.fa; 
	do python3 Spacer_parserV2.py --input $i;
	echo $i parsed
	done
	cat ./parsed/*.meta.txt >> ./parse_ultimeta.txt
fi

	 
##Align your spacers to reference. Not included in this loop- first making a BWA reference. This command is "bwa index ref.fa". You will also need 
##to download/install bwa 	 
## I had to align each spacers.fa file to a different reference genome, so i had to run an edited version of this short 
# script many times, once for each reference genome. An example of one of mine is below the shell of it
if [ $1 == "YourAlignmentCommand" ]; then
	for i in OnlyTheFilesYouWantToMap*.spacers.fa; 
	do 
	/Your/Path/To/bwa/bwa/ aln -n 0 /Your/Path/To/ReferenceFasta.fa $i >./align/$i.sai
	/Your/Path/To/bwa/bwa samse /Your/Path/To/ReferenceFasta.fa ./align/$i.sai $i > ./align/$i.sam
	done
fi	
	
 	#My Example
 if [ $1 == "alignNewman" ]; then
 	mkdir align
 	for i in ./parsed/*Newman*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	/home/modelllab/NAS/NGSraw/nkeith/tools/bwa/bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/Newman.fa $i > ./align/$filename2.sai
 	/home/modelllab/NAS/NGSraw/nkeith/tools/bwa/bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/Newman.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignNCTC" ]; then
 	mkdir align
 	for i in ./parsed/NCTC*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignTB4" ]; then
 	mkdir align
 	for i in ./parsed/TB4*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/TB4.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/TB4.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignRN_3" ]; then
 	mkdir align
 	for i in ./parsed/RN_phi11_phi12_phi13*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa ./align/$filename2.sai $i > ./align/$filename2.sam
	done
 fi
 
 if [ $1 == "alignNM4" ]; then
 	mkdir align
 	for i in ./parsed/RN_phiNM4*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM4.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM4.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignPhi12" ]; then
 	mkdir align
 	for i in ./parsed/RNphi12*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi12.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi12.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignPhi11" ]; then
 	mkdir align
 	for i in ./parsed/RN_phi11_P*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
  if [ $1 == "align450Phi11" ]; then
 	mkdir align
 	for i in ./parsed/RN450_phi11_P*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignPhi13" ]; then
 	mkdir align
 	for i in ./parsed/RN_phi13_P*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi13.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi13.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi
 
 if [ $1 == "alignNM2" ]; then
 	mkdir align
 	for i in ./parsed/RN_phiNM2*.spacers.fa; 
 	do 
 	filename="${i##*/}"
 	filename2="${filename%%.*}"
 	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM2.fa $i >./align/$filename2.sai
 	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM2.fa ./align/$filename2.sai $i > ./align/$filename2.sam
 	done
 fi

if [ $1 == "alignNM1" ]; then
	mkdir align
	for i in ./parsed/RN_phiNM1*.spacers.fa; 
	do 
	filename="${i##*/}"
	filename2="${filename%%.*}"
	bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM1.fa $i >./align/$filename2.sai
	bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM1.fa ./align/$filename2.sai $i > ./align/$filename2.sam
	done
fi

if [ $1 == "alignSF" ]; then
	mkdir align
	for i in ./parsed/*.spacers.fa; 
	do 
	filename="${i##*/}"
	filename2="${filename%%.*}"
 	/home/modelllab/NAS/NGSraw/nkeith/tools/bwa/bwa aln -n 0 /home/modelllab/NAS/NGSraw/nkeith/REF/SF370.fa $i > ./align/$filename2.sai
 	/home/modelllab/NAS/NGSraw/nkeith/tools/bwa/bwa samse /home/modelllab/NAS/NGSraw/nkeith/REF/SF370.fa ./align/$filename2.sai $i > ./align/$filename2.sam
	done
fi


##samtools is for sorting and indexing your alignment files. Alternatively, a python script using pySam has a lot of built in 
##utilities for playing around with alignment files. 

##putting your sam alignments into indexed bam form. You will need to download samtools. The base command is to keep file 
##extensions from getting out of order 
# Replace /Users/nicholaskeith/.../samtools with your path to samtools
if [ $1 == "samtools" ]; then
	for i in ./align/*.sam; 
	do 
	base=${i%%.*}
	filename="${i##*/}"
	filename2="${filename%%.*}"
	echo $filename2
	samtools view -b -o ./align/$filename2.bam $i
	samtools sort -o ./align/$filename2.sorted.bam ./align/$filename2.bam 
	samtools index ./align/$filename2.sorted.bam;
	done
fi

##just to output number of reads aligned and assigned to either plasmid or chromosome
# Replace /Users/nicholaskeith/.../samtools with your path to samtools
if [ $1 == "stats" ]; then
	touch bamStats.txt
	for i in ./align/*.sorted.bam;
	do 
	filename="${i##*/}"
	filename2="${filename%%.*}"
	echo $filename2 >> bamStats.txt
	samtools idxstats $i >> bamStats.txt;
	samtools flagstat $i >> bamStats.txt;
	echo >> bamStats.txt
	done
fi

##convert alignment to bed file for plotting in R
# Replace /Users/nicholaskeith/.../bedtools with your path to bedtools
if [ $1 == "bed" ]; then
	mkdir bed
	for i in ./align/JW*.sorted.bam;
	do
	filename="${i##*/}"
	filename2="${filename%%.*}"
	bedtools bamtobed -i $i &>./bed/$filename2.bed;
	done
fi

	
# A script to determine the number of canonical pams and spacers for a given set of reads. 
# You'll have to go into pamStats.py and change a lot of stuff for each reads file
# you want to do this for, sorry i don't know how to code very well
if [ $1 == "pamStats" ]; then
	mkdir pamStats
	touch pam_ultimeta.txt
	
# 	for i in ./bed/*Newman*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/Newman.fa --input $i -p1Name phiNM4 -p1Start 320753 -p1End 363941 -p2Name phiNM2 -p2Start 1098876 -p2End 1142020 -p3Name phiNM1 -p3Start 1980129 -p3End 2023256 -p4Name phiNM3 -p4Start 2088889 -p4End 2132221 
# 	done
# 	
# 	for i in ./bed/*NM4*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM4.fa --input $i -p1Name phiNM4 -p1Start 316258 -p1End 359446 
# 	done
# 	
# 	for i in ./bed/*NM2*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM2.fa --input $i -p1Name phiNM2 -p1Start 1042168 -p1End 1085312 
# 	done
# 	
# 	for i in ./bed/*NM1*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phiNM1.fa --input $i -p1Name phiNM1 -p1Start 1877364 -p1End 1920491 
# 	done
# 	
# 	for i in ./bed/RN_phi11_P*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa --input $i -p1Name phi11 -p1Start 1879324 -p1End 1922927 
# 	done
# 	
# 	for i in ./bed/RNphi12_*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi12.fa --input $i -p1Name phi12 -p1Start 1462550 -p1End  1508553
# 	done
# 	
# 	for i in ./bed/RN_phi13_P*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi13.fa --input $i -p1Name phi13 -p1Start 1942325 -p1End  1985050
# 	done
# 	
# 	for i in ./bed/RN_phi11_phi12_phi13*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa --input $i -p1Name phi12 -p1Start 1464415 -p1End 1508488 -p2Name phi11 -p2Start 1923408 -p2End 1967011 -p3Name phi13 -p3Start 2033113 -p3End 2074574
# 	done
# 	
# 	for i in ./bed/*NCTC*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325.fa --input $i -p1Name phi12 -p1Start 1464415 -p1End 1508488 -p2Name phi11 -p2Start 1923408 -p2End 1967011 -p3Name phi13 -p3Start 2033113 -p3End 2074574
# 	done
# 	
# 	for i in ./bed/RN450*.bed
# 	do
# 	echo $i
# 	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/NCTC8325_cured_phi11.fa --input $i -p1Name phi11 -p1Start 1879324 -p1End 1922927
# 	done
	
	for i in ./bed/JW*.bed
	do
	echo $i
	python ./pamStats.py --ref /home/modelllab/NAS/NGSraw/nkeith/REF/SF370.fa --input $i -p1Name 370.1 -p1Start 529587 -p1End 570505 -p2Name 370.2 -p2Start 778520 -p2End 821005 -p3Name 370.3 -p3Start 1189120 -p3End 1222649 -p4Name 370.4 -p4Start 1773339 -p4End 1786888 
	done
	
	cat ./pamStats/*Meta.txt >> ./pam_ultimeta.txt
	
fi

if [ $1 == "LBNstats" ]; then
	mkdir LBNStats
	touch LBN_ultimeta.txt
	
	for i in ./LBN/finalNumbers/*.tsv
	do
	echo $i
	python ./LBNstats.py  --input $i 
	done
	
	cat ./LBNStats/*LBNMeta.txt >>./LBN_ultimeta.txt
	
fi

##Optional step to remove reads that map to multiple parts of the reference from dataset
# Replace /Users/nicholaskeith/.../samtools with your path to samtools
if [ $1 == "rmmulti" ]; then
	for i in *.sorted.bam;
	do
	base=${i%%.*}
	/Users/nicholaskeith/seq/tools/samtools-1.9/samtools view -b -q 5 $i > $base.rmdup.bam
	done
fi
