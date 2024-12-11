#!/bin/env python

# To run this code, run "python pamStats.py --ref [path to reference] --input [path to spacer bed file]" and optionally -p if prophages are present 
# outputs bed file with added NGG column, y if that read is NGG PAM, n if not, unk if unknown
# stats output is printed in console
# only uses first fasta sequence in fasta file (if multiple)
# only analyzes first chromosome listed in bed file (if multiple)
## Make sure the following packages have been downloaded and can be imported (pandas and Biopython)
import pandas as pd
from Bio import SeqIO
from Bio import Seq
import argparse

#adding reference fasta, spacer bed file, and prophage arguments
parser = argparse.ArgumentParser(description='input reference fasta and spacer bed file, output number of canonical (NGG) PAMs')
parser.add_argument('--ref', dest= "ref", help="path of reference fasta file")
parser.add_argument('--input', dest= "bed", help="path of input bed file")

parser.add_argument("-p1Name", dest = "p1N", help="add prophage name")
parser.add_argument("-p1Start", dest="p1S", help="add prophage start nt")
parser.add_argument("-p1End", dest="p1E", help="add prophage end nt")

parser.add_argument("-p2Name", dest = "p2N", help="add prophage name")
parser.add_argument("-p2Start", dest="p2S", help="add prophage start nt")
parser.add_argument("-p2End", dest="p2E", help="add prophage end nt")

parser.add_argument("-p3Name", dest = "p3N", help="add prophage name")
parser.add_argument("-p3Start", dest="p3S", help="add prophage start nt")
parser.add_argument("-p3End", dest="p3E", help="add prophage end nt")

parser.add_argument("-p4Name", dest = "p4N", help="add prophage name")
parser.add_argument("-p4Start", dest="p4S", help="add prophage start nt")
parser.add_argument("-p4End", dest="p4E", help="add prophage end nt")

args = parser.parse_args()

fileName = str(str(args.bed).split("/")[-1])
fileName2 = str(str(fileName).split(".")[0])
meta = open(str("./pamStats/")+str(fileName2)+ '.pamMeta.txt','w')


print("")
anyProphage = False

if args.p1N is not None:
	prophageNameList = []
	prophageStartList = []
	prophageEndList = []
	prophageCanonList = []
	prophageNonCanonList = []
	prophageNameList.append(str(args.p1N))
	prophageStartList.append(int(args.p1S))
	prophageEndList.append(int(args.p1E))
	prophageCanonList.append(0)	
	prophageNonCanonList.append(0)
	numPro=1
	anyProphage = True

if args.p2N is not None:
	prophageNameList.append(str(args.p2N))
	prophageStartList.append(int(args.p2S))
	prophageEndList.append(int(args.p2E))
	prophageCanonList.append(0)	
	prophageNonCanonList.append(0)
	numPro=2
	anyProphage = True

if args.p3N is not None:
	prophageNameList.append(str(args.p3N))
	prophageStartList.append(int(args.p3S))
	prophageEndList.append(int(args.p3E))
	prophageCanonList.append(0)	
	prophageNonCanonList.append(0)
	numPro=3
	anyProphage = True

if args.p4N is not None:
	prophageNameList.append(str(args.p4N))
	prophageStartList.append(int(args.p4S))
	prophageEndList.append(int(args.p4E))
	prophageCanonList.append(0)	
	prophageNonCanonList.append(0)
	numPro=4
	anyProphage = True



counter = 0

#only pulls first FASTA record aka reference sequence from record file (not pGG32ultracr)
#fastaRecord = next(SeqIO.parse(args.ref, "fasta"))
fastaRecord = list(SeqIO.parse(args.ref, "fasta"))

print("Analyzing "+ str(args.bed))


fastaNum = 0


#fastaNum = counter
fasta = fastaRecord[fastaNum].seq
print("Using " + fastaRecord[counter].id + " as reference.")


#imports in your BED file of reads
wholebed = pd.read_csv(args.bed, delimiter = "\t", names = ["chrom", "start", "stop", "name", "score", "strand"])
#to only read spacers on first chromosome in bed file
myGenome = fastaRecord[fastaNum].id
print("Using " + str(myGenome) +"...")

bed = wholebed.drop(wholebed[(wholebed.chrom != myGenome)].index)



#add columns
nggList = []
prophageList = []
spacerList=[]
spacer = ""

pam = str()

# Number of reads (rows) in your bed file
maxBed = len(bed.index)

# Reads
total_reads = 0
canonGenome = 0
noncanonGenome = 0
canonPro = 0
noncanonPro = 0
ngg = False

# If the read is in a prophage region
pro = False

#current row of the bed file
rw = 0

# starts at row(read)0 of bed file, ends when you're no longer reading lines from the first chromosome in the bed file 
# or at last bed file row
# pulls in the start/end/strand of each read (make sure your bed file matches the intended columns)
# if the start and end of a read is within a prophage region, pro designated as true
# sets current phage as the phage spacer is found in 
# if + strand then downstream pam must be gg or its noncanon, 
# otherwise must be cc upstream of start because on - strand, and fasta is top strand only
# if within a prophage region, added as canon or noncanon for prophage, 
# also adds read to separate list of individual prophage spacer counts, if multiple prophages are present
# added functionality to include the full spacer from 5' to 3', so that the terminal nucleotide is the one for last base normalization (if doing)

print("\n" + "Determining PAMs...")



while (rw != maxBed) and (bed.iloc[rw,0] == myGenome):

	total_reads +=1 
	chrStart = bed.iloc[rw,1]
	chrEnd = bed.iloc[rw, 2]
	strand = bed.iloc[rw, 5]
# to return 5' to 3' spacer sequence
	spacer = fasta[chrStart:chrEnd]
	if strand == '-':
		spacerOpp = Seq.Seq(str(spacer))
		spacer = spacerOpp.reverse_complement()

	spacerList.append(str(spacer))
	
	if anyProphage:
		for i in range(numPro):
			if chrStart > int(prophageStartList[i]) and chrEnd < int(prophageEndList[i]):
				pro = True
				currentPro = i
				prophageList.append(str(prophageNameList[currentPro]))
				break
			else:
				pro = False
		if not pro:
			prophageList.append('NA')
	if strand == '+':
		pam = fasta[chrEnd +1: chrEnd + 3]
		if pam.lower() == "gg":
			nggList.append('y')
			if pro:
				canonPro += 1
				prophageCanonList[currentPro] += 1
			else:
				canonGenome +=1
		else:
			nggList.append('n')
			if pro:
				noncanonPro += 1
				prophageNonCanonList[currentPro] += 1
			else:
				noncanonGenome += 1
	else:
		pam = fasta[chrStart - 3 : chrStart - 1]
		if pam.lower() == "cc":
			nggList.append('y')
			if pro:
				canonPro += 1
				prophageCanonList[currentPro] += 1
			else:
				canonGenome += 1
		else:
			nggList.append('n')
			if pro:
				noncanonPro += 1
				prophageNonCanonList[currentPro] += 1
			else:
				noncanonGenome +=1
	if (rw % 100000) == 0:
		print (str(rw))
	rw += 1


bed['NGG'] = nggList
bed['spacer'] = spacerList
if anyProphage:
	bed['Phage'] = prophageList

# writing our outputs!	
meta.write("\t"+ args.bed.rsplit("/", 1)[-1])
meta.write("\nTotal Spacers:\t" + str(total_reads))
if total_reads != 0:
	perCanTot = (float(canonGenome) + float(canonPro)) / float(total_reads) * 100
	meta.write('\nPercent NGG Spacers from Total=\t' + "{:.2f}".format(perCanTot) + "%")

meta.write("\n\nFrom Bacterial Genome:")
meta.write("\nNGG Spacers:\t" + str(canonGenome))
meta.write("\nNon-NGG Spacers:\t" + str(noncanonGenome))
if canonGenome and noncanonGenome !=0:
	perCanGen = float(canonGenome)/float(canonGenome + noncanonGenome) * 100
	meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanGen) + '%')

if anyProphage:
	meta.write("\n\nFrom Prophage Regions:")
	meta.write("\nNGG Spacers:\t" + str(canonPro))
	meta.write("\nNon-NGG Spacers:\t" + str(noncanonPro))
	if canonPro and noncanonPro !=0:
		perCanPro = float(canonPro)/float(canonPro + noncanonPro) * 100
		meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanPro) + "%")
		changingNumber = float(canonPro + noncanonPro) / float(total_reads) * 100 
		meta.write("\nPercentage of spacers from prophages out of total:\t" + "{:.2f}".format(changingNumber) + "%")
	
	if numPro > 1:
		for i in range(numPro):
			meta.write("\n\nFrom prophage " + str(prophageNameList[i]) + ":")
			meta.write("\nNGG Spacers:\t" + str(prophageCanonList[i]))
			meta.write("\nNon-NGG Spacers:\t" + str(prophageNonCanonList[i]))
			if (float(prophageCanonList[i]) + float(prophageNonCanonList[i])) != 0:
				changingNumber = float(prophageCanonList[i] / (float(prophageCanonList[i]) + float(prophageNonCanonList[i]))) * 100
				meta.write("\nPercent NGG:\t" + "{:.2f}".format(changingNumber) + "%")
				changingNumber = float(prophageCanonList[i] + prophageNonCanonList[i])/float(canonPro + noncanonPro) *100
				meta.write("\nPercentage of prophage spacers specific for " + str(prophageNameList[i]) +":\t" "{:.2f}".format(changingNumber) + "%")



fileName3 = str("./pamStats/"+ str(fileName2)+'.tsv')

print("Saving as "+ str(fileName3))
bed.to_csv(str(fileName3), sep='\t', index=False)
meta.write("\n")
meta.close()

print("")
print("")



