import vcfpy
import argparse
from Bio.Seq import Seq


parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-vcf','--input_vcf', metavar='', help='input vcf file' ,default='', type =str,nargs=1)
parser.add_argument('-rc','--rev_com', action='store_true', help='reverse complement')
parser.add_argument('-mt','--mt_code', action='store_true', help='use fungal mt code for translating')

#parser.add_argument('-r','--number_of_replicates', metavar='',default=[1], help='how many times to run the admixture analysis for each K (default=1)', nargs='*',type=int)
#parser.add_argument('-anc','--ancestor', metavar='',default='', help='list of outgroups', nargs=1,type=str)
parser.add_argument('-o','--out', default='', metavar='',help='stem for output file', nargs= 1,type=str)
#parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)



arguments = parser.parse_args()




# Open vcf file, this will read in the header
reader = vcfpy.Reader.from_path(arguments.input_vcf[0])

# Build and print header
header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names




#for i in range(len(reader.header.samples.names)):
#	print (i)
#		if sample == reader.header.samples.names[i]:

# initialize fasta file
OUT={}


for record in reader.header.samples.names:
	OUT.update({str(record):""})


#print (OUT)


for record in reader:
	geno = [call.data.get('GT') or './.' for call in record.calls]
	list_geno = list(geno)

	n_alt = len(record.ALT)


#	print(n_alt)

#		if len(record.ALT[0].value) > 1:
#	if len(record.REF) > 1:

	for z in range(len(geno)):
		if geno[z] == "0":
			geno[z] = str(record.REF[0])
		elif geno[z] == "./." or geno[z] == ".":
			geno[z] = "N"
		else:
			for i in range(n_alt+1):
				if str(geno[z]) == str(i):
					geno[z] = str(record.ALT[i-1].value)


	pairs = zip(reader.header.samples.names, geno)

	my_dict = dict(pairs)

	for name in reader.header.samples.names:
		OUT[name] += my_dict[name]

f = open(str(arguments.out[0])+".fa", "w")
for name_out, seq_out in OUT.items():
#		print(name_out)
#		print(seq_out)
		f.write(">" + name_out + "\n")
		if (arguments.rev_com):
			seq=Seq(seq_out)
			rev_seq=seq.reverse_complement()
			f.write(str(rev_seq)+"\n")
		else:
			f.write(seq_out+"\n")
f.close()

f = open(str(arguments.out[0])+"_p.fa", "w")

if (arguments.mt_code):
	codon_table=4
else:
	codon_table=1

print (codon_table)

for name_out, seq_out in OUT.items():
#		print(name_out)
#		print(seq_out)
		f.write(">" + name_out + "\n")
		seq=Seq(seq_out)
		if (arguments.rev_com):
			rev_seq=seq.reverse_complement()
			protein_seq = rev_seq.translate(table = codon_table)

			f.write(str(protein_seq)+"\n")
		else:
			protein_seq = seq.translate(table = codon_table)
			f.write(str(protein_seq)+"\n")
f.close()


