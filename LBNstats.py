import pandas as pd
from Bio import SeqIO
from Bio import Seq
import argparse

#adding reference fasta, spacer bed file, and prophage arguments
parser = argparse.ArgumentParser(description='input reference fasta and spacer bed file, output number of canonical (NGG) PAMs')
parser.add_argument('--ref', dest= "ref", help="path of reference fasta file")
parser.add_argument('--input', dest= "bed", help="path of input bed file")

args = parser.parse_args()

fileName = str(str(args.bed).split("/")[-1])
fileName2 = str(str(fileName).split(".")[0])




#imports in your BED file of reads
bed = pd.read_csv(args.bed, delimiter = "\t", names = ["source", "chrom", "start", "stop", "strand", "NGG", "spacer", "Phage", "unique", "normed", "newTotal", "RPM"])

# Number of reads (rows) in your bed file
maxBed = len(bed.index)

print("\n" + "Determining PAMs...")

#declaring shit

rw = 1
currentSource = bed.iloc[rw,0]
anyProphage = False
canonPro = 0
noncanonPro = 0
canonGenome = 0
noncanonGenome = 0
total_reads = 0

#for each line in file
while rw < maxBed:

	#while using the same sample 
	if bed.iloc[rw,0] == currentSource:
	
		#row 5 = NGG y or n, row 10 = reads normalized to RPM, row 7 is either NA or a phage name
		total_reads = total_reads + float(bed.iloc[rw,11])
		if bed.iloc[rw,5] == 'y':
			if str(bed.iloc[rw,7]) == "nan":
				canonGenome = canonGenome + float(bed.iloc[rw,11])
			else:
				canonPro = canonPro + float(bed.iloc[rw,11])
				anyProphage = True
		else:
			if str(bed.iloc[rw,7]) == "nan":
				noncanonGenome = noncanonGenome + float(bed.iloc[rw,11])
			else:
				noncanonPro = noncanonPro + float(bed.iloc[rw,11])
				anyProphage = True
		rw +=1
	else: 
		fileName = str(str(bed.iloc[rw-1,0]).split("/")[-1])
		fileName2 = str(str(fileName).split(".")[0])
		meta = open(str("./LBNStats/")+str(fileName2)+ '.pamLBNMeta.txt','w')
		meta.write("\t"+ fileName2)
		meta.write("\nTotal Spacers:\t" + "{:.2f}".format(total_reads))
		if total_reads != 0:
			perCanTot = (float(canonGenome) + float(canonPro)) / float(total_reads) * 100
			meta.write('\nPercent NGG Spacers from Total=\t' + "{:.2f}".format(perCanTot) + "%")

		meta.write("\n\nFrom Bacterial Genome:")
		meta.write("\nNGG Spacers:\t" + "{:.2f}".format(canonGenome))
		meta.write("\nNon-NGG Spacers:\t" + "{:.2f}".format(noncanonGenome))
		if canonGenome and noncanonGenome !=0:
			perCanGen = float(canonGenome)/float(canonGenome + noncanonGenome) * 100
			meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanGen) + '%')

		if anyProphage:
			meta.write("\n\nFrom Prophage Regions:")
			meta.write("\nNGG Spacers:\t" + "{:.2f}".format(canonPro))
			meta.write("\nNon-NGG Spacers:\t" + "{:.2f}".format(noncanonPro))
			if canonPro and noncanonPro !=0:
				perCanPro = float(canonPro)/float(canonPro + noncanonPro) * 100
				meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanPro) + "%")
				changingNumber = float(canonPro + noncanonPro) / float(total_reads) * 100 
				meta.write("\nPercentage of spacers from prophages out of total:\t" + "{:.2f}".format(changingNumber) + "%")
		meta.write("\n\n")
		meta.close()
		currentSource = bed.iloc[rw,0]
		anyProphage = False
		canonPro = 0
		noncanonPro = 0
		canonGenome = 0
		noncanonGenome = 0
		total_reads = 0

				


fileName = str(str(bed.iloc[rw-1,0]).split("/")[-1])
fileName2 = str(str(fileName).split(".")[0])
meta = open(str("./LBNStats/")+str(fileName2)+ '.pamLBNMeta.txt','w')
meta.write("\t"+ fileName2)
meta.write("\nTotal Spacers:\t" + "{:.2f}".format(total_reads))
if total_reads != 0:
	perCanTot = (float(canonGenome) + float(canonPro)) / float(total_reads) * 100
	meta.write('\nPercent NGG Spacers from Total=\t' + "{:.2f}".format(perCanTot) + "%")

meta.write("\n\nFrom Bacterial Genome:")
meta.write("\nNGG Spacers:\t" + "{:.2f}".format(canonGenome))
meta.write("\nNon-NGG Spacers:\t" + "{:.2f}".format(noncanonGenome))
if canonGenome and noncanonGenome !=0:
	perCanGen = float(canonGenome)/float(canonGenome + noncanonGenome) * 100
	meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanGen) + '%')

if anyProphage:
	meta.write("\n\nFrom Prophage Regions:")
	meta.write("\nNGG Spacers:\t" + "{:.2f}".format(canonPro))
	meta.write("\nNon-NGG Spacers:\t" + "{:.2f}".format(noncanonPro))
	if canonPro and noncanonPro !=0:
		perCanPro = float(canonPro)/float(canonPro + noncanonPro) * 100
		meta.write('\nPercent NGG:\t' + "{:.2f}".format(perCanPro) + "%")
		changingNumber = float(canonPro + noncanonPro) / float(total_reads) * 100 
		meta.write("\nPercentage of spacers from prophages out of total:\t" + "{:.2f}".format(changingNumber) + "%")
meta.write("\n\n")
meta.close()


	
				
			

