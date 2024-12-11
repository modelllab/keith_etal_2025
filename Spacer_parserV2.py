from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import argparse

parser = argparse.ArgumentParser(description='input spacer aquisition sequencing reads, output parsed spacers')
parser.add_argument('--input', dest= "input", help="name of input file")

args = parser.parse_args()

input = open(args.input, "r")

outputName = str(args.input).split("/")
outputName2 = str(outputName[-1]).split(".")
outputName3 = str(outputName2[0])

output = open("./parsed/"+ outputName3 +".spacers.fa", "w")
badspacers = open("./parsed/"+ outputName3 +".badspacers.txt", "w")
meta = open("./parsed/"+ outputName3 +".meta.txt", "w")
multi = open("./parsed/"+ outputName3 +".multi.txt", "w")


repB = "CTCTAAAA"
repE = "GTTTTGGG"
repES = "CTCTAAAACCTCGTAG"

totalreads = 0
totalspacers = 0
totalunadapted = 0
totalmulti = 0
totalbadspacers = 0

def my_rev_complement(seq):
	return Seq(seq).reverse_complement()

for record in SeqIO.parse(input,'fasta'):
        totalreads += 1
        if((totalreads % 100000) == 0):
                print(totalreads)
        posBf = nt_search(str(record.seq), repB)
        if len(posBf) > 1:
                posEf = nt_search(str(record.seq), repE)
                posESf = nt_search(str(record.seq), repES)
                if len(posEf) > 1:
                		if (len(posEf) > 2) and (len(posBf) > 2):
                				totalmulti += 1
                				spacer = str(record.seq)[posBf[1]+9:posEf[1]]
                				spacer_rev = my_rev_complement(spacer)
                				spacer2 = str(record.seq)[posBf[2]+9:posEf[2]]
                				spacer2_rev = my_rev_complement(spacer2)
                				multi.write(">" + str(totalmulti) + "\n" + str(spacer_rev) + "\t" + str(spacer2_rev) + "\n")
                		else:
                        		spacer = str(record.seq)[posBf[1]+9:posEf[1]]
                        		spacer_rev = my_rev_complement(spacer)
                        		if (len(spacer)>5):
                        				totalspacers += 1
                        				output.write(">" + str(totalspacers) + "\n" + str(spacer_rev) + "\n")
                        		else: 
                        				totalbadspacers +=1
                        				badspacers.write(str(record.seq) + "\n")
                elif (len(posESf) > 1) and (posESf[1] < 40):
                		totalunadapted +=1
                else:
                		totalbadspacers += 1
                		badspacers.write(str(record.seq) + "\n")
        else:
                totalbadspacers += 1
                badspacers.write(str(record.seq) + "\n")

meta.write("\t" + outputName3 + "\n")                        
meta.write("total reads:\t" + str(totalreads) + "\n")
meta.write("total spacers:\t" + str(totalspacers) + "\n")
meta.write("total undadapted:\t" + str(totalunadapted) +"\n")
meta.write("total multi-acquisition:\t" + str(totalmulti)+"\n")
meta.write("total bad spacers:\t" + str(totalbadspacers) + "\n")

meta.close()
input.close()
output.close()
multi.close()
badspacers.close()