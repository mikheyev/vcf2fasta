from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from copy import deepcopy
import sys
import pdb
import bisect, re
from collections import defaultdict
"""
Takes a an original fasta file, a vcf file, a gff file and a sample number (starting from zero) and produces a new fasta and gff, where each coordinate is mapped to a new one.
Treats sites with more than one alt as Ns.
This version works for haploid samples only.
"""

GAPPEDAMBIG_DNA_ALPHABET = Alphabet.Gapped (IUPAC.ambiguous_dna, '-')
genome = defaultdict(dict)
for rec in SeqIO.parse("test.fa","fasta",alphabet=GAPPEDAMBIG_DNA_ALPHABET):
	genome["ref"][rec.id] = rec.seq.tomutable()

s = 2  #FIXME: this should be sysargv

	#                         __                                                          ___             __           
	#   _____________  ____ _/ /____     ____  ___ _      __   _________  ____  _________/ (_)___  ____ _/ /____  _____
	#  / ___/ ___/ _ \/ __ `/ __/ _ \   / __ \/ _ \ | /| / /  / ___/ __ \/ __ \/ ___/ __  / / __ \/ __ `/ __/ _ \/ ___/
	# / /__/ /  /  __/ /_/ / /_/  __/  / / / /  __/ |/ |/ /  / /__/ /_/ / /_/ / /  / /_/ / / / / / /_/ / /_/  __(__  ) 
	# \___/_/   \___/\__,_/\__/\___/  /_/ /_/\___/|__/|__/   \___/\____/\____/_/   \__,_/_/_/ /_/\__,_/\__/\___/____/  


#dictionaries of lists have the new coordinates
scaffolds_ivals = dict()
scaffolds_offsets = dict()
for line in open("test.vcf"):
	if line[0] == "#":	#skip over comments
		if line.find("#CHROM") == -1:
			continue
		samples = line.rstrip().split("\t")[9:]
		#make a copy of reference genome for every in vcf
		for sample in samples:
			for i in genome["ref"]:
				genome[sample][i] = deepcopy(genome["ref"][i])
		continue
	(chrom, pos, ID, ref, alt, qual, filt, info) = line.split("\t")[0:8]
	if chrom not in scaffolds_offsets:  # set the first offset to zero
		scaffolds_offsets[chrom] = []
		scaffolds_ivals[chrom] = []
		offset = 0
#	if pos == '6524':
#		pdb.set_trace()
	pos = int(pos)-1 + offset #vcfs are 1-indexed
	if genome["ref"][chrom][pos] != ref[0]:		
		pdb.set_trace()
	#extend reference if there is a longer insert
	maxalt = max(map(len,alt.split(",")))
	if maxalt > len(ref):
#		pdb.set_trace()
		genome["ref"][chrom][pos:pos+len(ref)] = ref+"-"*(maxalt - len(ref))
	else:
		maxalt = len(ref)
	#extend other alts to match the longest one
	snps = []
	for snp in alt.split(","):
		if len(snp) < maxalt:
			snp += (maxalt - len(snp))*"-"
		snps.append(snp)
	alt = ",".join(snps)
	#adjust genomes 
	variants = line.rstrip().split("\t")[9:]
	for idx,sample in enumerate(samples):
	# deal with the various possible cases of substitutions
		if variants[idx].split(":")[0] == "." :  #no data 
			genome[sample][chrom][pos:pos+len(ref)] = len(ref)*"N"+"N"*(maxalt - len(ref))
			continue
			#genome[sample][chrom][pos:pos+len(ref)] = "N"*len(ref)
		elif variants[idx].split(":")[0] == "0":
			genome[sample][chrom][pos:pos+len(ref)] = ref+"-"*(maxalt - len(ref))
		elif variants[idx].split(":")[0] != "0":
			currentAlt = alt.split(",")[int(variants[idx].split(":")[0])-1]
			genome[sample][chrom][pos:pos+len(ref)] = currentAlt+"-"*(maxalt - len(currentAlt))
	offset += (maxalt - len(ref))
	scaffolds_ivals[chrom].append(int(line.split("\t")[1])-1)  #convert oritinal position to VCF coordinates
	scaffolds_offsets[chrom].append(offset)

	#                                   __     _    ______________                             ___             __           
	#   _________  _____________  _____/ /_   | |  / / ____/ ____/  _________  ____  _________/ (_)___  ____ _/ /____  _____
	#  / ___/ __ \/ ___/ ___/ _ \/ ___/ __/   | | / / /   / /_     / ___/ __ \/ __ \/ ___/ __  / / __ \/ __ `/ __/ _ \/ ___/
	# / /__/ /_/ / /  / /  /  __/ /__/ /_     | |/ / /___/ __/    / /__/ /_/ / /_/ / /  / /_/ / / / / / /_/ / /_/  __(__  ) 
	# \___/\____/_/  /_/   \___/\___/\__/     |___/\____/_/       \___/\____/\____/_/   \__,_/_/_/ /_/\__,_/\__/\___/____/  
	#                                                                                                                       


#AlignIO.write( MultipleSeqAlignment([SeqRecord(seq = genome[i]['pbar_scf7180000350264'].toseq(),id=i) for i in genome]), "test.phy", "phylip")
#pdb.set_trace()
SeqIO.write(SeqRecord(seq = genome['ref']['pbar_scf7180000350264'].toseq(),id='pbar_scf7180000350264'),"pbar_scf7180000350264.fa","fasta")
outfile = open('pbar_scf7180000350264.gff',"w")
for line in open("test.gff"):
	if line[0] == "#" :	#skip over comments
		continue
	line = line.split("\t")
#	if line[3] == '13048': #'4810':
#		pdb.set_trace()
	offset = scaffolds_offsets[line[0]][bisect.bisect_left(scaffolds_ivals[line[0]],int(line[3]))]
#	pdb.set_trace()
	offset_start = int(line[3]) + offset
#	if int(line[4]) > scaffolds_ivals[line[0]][-1]:
#		offset_end = len(genome[line[0]])
#	else:
	offset = scaffolds_offsets[line[0]][bisect.bisect_left(scaffolds_ivals[line[0]],int(line[4]))]
	offset_end = int(line[4])+offset
	line[3] = str(offset_start)
	line[4] = str(offset_end)
	print >>outfile, "\t".join(line),
outfile.close()
	#                 _ __          ____           __       
	#  _      _______(_) /____     / __/___ ______/ /_____ _
	# | | /| / / ___/ / __/ _ \   / /_/ __ `/ ___/ __/ __ `/
	# | |/ |/ / /  / / /_/  __/  / __/ /_/ (__  ) /_/ /_/ / 
	# |__/|__/_/  /_/\__/\___/  /_/  \__,_/____/\__/\__,_/  
	#                                                       

# outfile = open("recoded.fa","w")
# for rec in genome:
# 	print >>outfile, ">"+rec
# 	print >>outfile, genome[rec]
# outfile.close()
