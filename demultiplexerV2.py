
  
# Demultiplex internal indexes
# change output names to whatever the 'f' barcodes are for each i12
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import argparse
import sys


parser = argparse.ArgumentParser(description='input iPrimer raw NGS reads, output demuxed fwd reads')
parser.add_argument('--input', dest= "input", help="name of input file")
parser.add_argument('--convert', dest = "convert", help="convert from fastq to fasta only, with typed suffix")
parser.add_argument('--f1', dest = "f1", help="what you want to name the f1 file")
parser.add_argument('--f2', dest = "f2", help="what you want to name the f2 file")
parser.add_argument('--f3', dest = "f3", help="what you want to name the f3 file")
parser.add_argument('--f4', dest = "f4", help="what you want to name the f4 file")
parser.add_argument('--f5', dest = "f5", help="what you want to name the f5 file")
parser.add_argument('--f6', dest = "f6", help="what you want to name the f6 file")

args = parser.parse_args()

namerino = str(args.input).split("_")
prefix = str(namerino[0])

input = open(args.input,"r")

trimFile = input.name + ".trimFA"
trimmed = SeqIO.convert(input, "fastq", trimFile, "fasta")

records = SeqIO.parse(trimFile,'fasta')

if args.convert is not None:
	SeqIO.write(trimmed,  str("./") +str(input.name)+ ".fa", "fasta")
	sys.exit("files converted from fastq to fasta")



meta = open(str("./demultiplexed/")+str(prefix)+ 'meta.txt','w')



bar1 = 'ATCA'
bar2 = 'TTAG'
bar3 = 'ACTT'
bar4 = 'GATC'
bar5 = 'TAGC'
bar6 = 'GGCT'

total_reads = 0
output1_reads = 0
output2_reads = 0
output3_reads = 0
output4_reads = 0
output5_reads = 0
output6_reads = 0
noBarcode = 0


out1 = []
out2 = []
out3 = []
out4 = []
out5 = []
out6 = []
outNoBar = []

print("\nDemultiplexing " + str(input.name))
print("\nReads Demultiplexed:")

for r in records:
        total_reads += 1
        if total_reads %50000 == 0:
        	print(total_reads)
        seq = str(r.seq)
        barcode = seq[0:4]
        if barcode == bar1:
                output1_reads += 1
                out1.append(r)
        elif barcode == bar2:
                output2_reads += 1
                out2.append(r)
        elif barcode == bar3:
                output3_reads += 1
                out3.append(r)
        elif barcode == bar4:
                output4_reads += 1
                out4.append(r)
        elif barcode == bar5:
                output5_reads += 1
                out5.append(r)
        elif barcode == bar6:
                output6_reads += 1
                out6.append(r)
        else:
        		noBarcode += 1
        		outNoBar.append(r)



SeqIO.write(out1,  str("./demultiplexed/")+str(args.f1) + ".fa", "fasta")

if args.f2 is not None:
	SeqIO.write(out2, str("./demultiplexed/")+ str(args.f2) + ".fa", "fasta")
if args.f3 is not None:
	SeqIO.write(out3,  str("./demultiplexed/")+str(args.f3) + ".fa", "fasta")
if args.f4 is not None:
	SeqIO.write(out4,  str("./demultiplexed/")+str(args.f4) + ".fa", "fasta")
if args.f5 is not None:
	SeqIO.write(out5,  str("./demultiplexed/")+str(args.f5) + ".fa", "fasta")
if args.f6 is not None:
	SeqIO.write(out6, str("./demultiplexed/")+str(args.f6) + ".fa", "fasta")
	
SeqIO.write(outNoBar,str("./demultiplexed/")+ str(prefix) +"_noBar.fa", "fasta")



print(str(args.input) + " demultiplexed")

meta.write("\t" + str(input.name) + "\n")                        
meta.write("total reads:\t" + str(total_reads) + "\n")
meta.write("total "+ str(args.f1) +":\t" + str(output1_reads) + "\n")
if args.f2 is not None:
	meta.write("total "+ str(args.f2)+":\t" + str(output2_reads) + "\n")
if args.f3 is not None:
	meta.write("total "+ str(args.f3)+":\t" + str(output3_reads) + "\n")
if args.f4 is not None:
	meta.write("total "+ str(args.f4)+":\t" + str(output4_reads) + "\n")
if args.f5 is not None:
	meta.write("total "+ str(args.f5)+":\t" + str(output5_reads) + "\n")
if args.f6 is not None:
	meta.write("total "+ str(args.f6)+":\t" + str(output6_reads) + "\n")

meta.write("total no barcode:\t"+ str(noBarcode) + "\n")

                                
input.close()
meta.close()

print("Outputs in ./demultiplexed/")

