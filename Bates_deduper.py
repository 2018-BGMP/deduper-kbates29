import argparse
import numpy as np
import re


def get_arguments():
	parser = argparse.ArgumentParser(description= "Takes aligned reads and removes PCR duplicates, returns a file containing only unique reads. This program will only work with single end reads.")
	parser.add_argument( "-f", "--file",
	help="The full path of the file you wish to eliminate PCR duplicates from.",
	required=True, type=str)
	parser.add_argument("-p", "--paired",
	help = "Specifies if your input file is paired end reads. If passed an error message will print. Reads must be single end in order to work properly",
	action="store_true", required = False)
	parser.add_argument("-u", "--umi",
	help = " File which conatins the unique molecular identifiers used account for PCR duplicates",
	required= False, type = str)
	return parser.parse_args()

args = get_arguments()
file = args.file
umi = args.umi
if args.paired:
	print("Error: Program cannot handle paired end reads.")
	exit()



def start_fixer(start, cigar, flag):
	'''
	Checks the start postion for the forward and reverse strand.
	If forward strand, corrects for soft clipping.
	If reverse strand, corrects for soft clipping, insertions, deletions, and skipped regions.
	'''
	input_start = int(start)
	cig_string = cigar
	flag = int(bit)
	if (flag & 16) != 16:
		# forward strand
		soft = re.search("^([0-9]+[S])", cig_string)
		if soft != None:
			# soft clipping occured, adjust start leftwards
			sub = soft.group(1)
			sub = int(sub[:-1])
			correct_start = input_start - sub
			str_start = str(correct_start)
		else:
			# no soft clipping, start is correct
			str_start = str(input_start)
	else:
		# reverse strand
		# can ingnore insertions
		match = re.findall("[0-9]+[M]", cig_string)
		for m in match: # matches or mismatches
			m_fix = int(m[:-1])
			input_start = input_start + m_fix
		soft = re.search("([0-9]+[S])$", cig_string)
		if soft != None: # soft clipping
			sub = soft.group(1)
			sub = int(sub[:-1])
			input_start = input_start + sub
		deletions = re.search("([0-9]+[D])", cig_string)
		if deletions != None: # deletions
			delete = deletions.group(1)
			delete = int(delete[:-1])
			input_start = input_start + delete
		skip = re.search("([0-9]+[N])", cig_string)
		if skip != None: # skipped regions
			skipped = skip.group(1)
			skipped = int(skipped[:-1])
			input_start = input_start + skipped
		str_start = str(input_start)
	return(str_start)




###################
## Main Function ##
###################

outfile = re.search("([A-Za-z0-9]+\.sam)$", file).group(1) # creates outfile handle


if umi != None: # looks to see if umi option has been specified
	with open(umi, "r") as um:
		umi_dict = {}
		for line in um:
			umi_dict[line[:-1]] = {} # initialze dictionary where key is UMI and value is another dictionary which will contain tuple of key parameters to check
else:
	rand_dict = {} # initialze dictionary where key is Randomer and value is another dictionary where key is a tuple of parameters to check for PCR duplicate




with open(file, "r") as fh:
	with open(outfile + "_deduped", "w") as dd:
		for line in fh:
			if umi != None:
				# umi has been given
				if line.startswith("@"):
					dd.write(line)
				else:
					read = line.strip().split()
					read_array = np.array(read)
					qname = read[0] # qname for line
					umi = qname[-8:] # gets umi from last 8 chars of qname
					chrom = read[2] # gets chrom for read
					bit = read[1] # gets bitwise flag for read
					print(bit)
					start = read[3] # gets uncorrected start position for read
					cigar = read[5] # gets cigar string for read
					if umi in umi_dict:
						real_start = start_fixer(start, cigar, bit)
						stats = (bit, chrom, real_start) # creates tuple for key paramerters to look at for PCR duplicate
						if stats in umi_dict[umi]:
							# PCR duplicate can discard read
							continue
						else:
							# unique read add to dictionary and wrtie to _deduped
							umi_dict[umi][stats] = ""
							dd.write(line)
					else:
						continue # bad read can discard
			else:
				# randomer option specified
				if line.startswith("@"):
					dd.write(line)
				else:
					read = line.strip().split()
					read_array = np.array(read)
					qname = read[0]
					umi = qname[-8:]
					chrom = read[3]
					bit = read[2]
					start = read[4]
					cigar = read[6]
					if randomer in rand_dict:
						real_start = start_fixer(start, cigar, bit)
						stats = (bit, chrom, real_start)
						if stats in rand_dict[randomer]:
							# PCR duplicate can discard read
							continue
						else:
							# unique read add to dictionary and wrtie to _deduped
							rand_dict[randomer][stats] = ""
							dd.write(line)
					else:
						# new randomer write to file
						rand_dict[randomer] = {}
						real_start = start_fixer(start, cigar, bit)
						stats = (bit, chrom, real_start)
						rand_dict[randomer][stats] = ""
						dd.write(line)
