#!/usr/bin/env python

####Argeprase Setup####


#Importing argparse so that the program can be run from the command line, and re for identifying regular expressions for later use in the script. 
import re
import argparse


#Desiginating the required and optional arguments for this program. These include the sorted SAM file, a flag to include a known UMI file,
#and an option for paired end reads (which is not included in this iteration of the program). 

parser = argparse.ArgumentParser(description="Given a sorted SAM file (samtools sort, this program removes all of the PCR duplicates present and returns the file, in it's entirety, minus the duplicates.")
parser.add_argument("-f", "--file_path", help="Requires the unsorted SAM file that needs to be deduped.")
parser.add_argument("-u", "--UMI" ,help="This optional flag indicates that known UMIs (not randomers) are used in the SAM file and requires an UMI file to be passed to it.")
parser.add_argument("-p", "--paired_end", help="Indicates that the file is paired end", action="store_true")
args = parser.parse_args()







#Checking to see if a file is supplied for the known UMIs, if there is, create a list of known UMIs that will be later used in the program. 
#If a UMI file is not provided, print an error and exit the program.
if args.UMI:
	UmiList = []
	with open(args.UMI) as fh100:
	    for line in fh100:
		    line = line.strip()
		    UmiList.append(line)
else:
	print('Warning! Program is not designed to handle randomers. A file of known UMIs is required.')
	exit()

#Checking to see if the paired-end option was specified from the commmand line. If it is, print an error message and exit the program.
if args.paired_end:
	print("Warning! Program is not designed to handle paired end reads.")
	exit()
	
####Function Definitions####

	
def UmiChecker(UMI):
    '''Takes the UMI from each entry and checks to see if it is in the list of known UMIs provided. If the UMI is recognized, the function returns true, if not, it returns false.'''
    if UMI in UmiList:
        return(True)
    else:
        return(False)

def bit_checker(bit):
	'''Takes the bitwise flag and checks it for strandedness. Assumes read is mapped, otherwise returns None. Assume data are single-end reads.Returns “+” or “-”, depending on strand.'''
	strand = '+'
	if (bit & 4) == 4:
		print("Read is unmapped!")
		return None
	if (bit & 16) == 16:
		strand = '-'
	return strand
 		
def softclip_adjustment(starting_position,adjustment):
	'''Takes the original starting position of a soft clipped read, and adjusts it for accurate comparisons with other reads'''
	position = 0
	position = int(starting_position) - int(adjustment)
	return position

####Actual Script####
	
#Opening the sorted SAM file and a new file to write all the non-duplicate reads to desiginating them as fh and fh2 respectively.	
with open(args.file_path) as fh:
    with open(args.file_path + "_deduped", 'w') as fh2:
		#Establishing the dictionary that will be used to store all the reads for comparison and setting counters for PCR Duplicates, Non-Duplicates, and reads with incorrect UMIs.
        testdict = {}
        PCRDuplicates = 0
        NotPCRDuplicates = 0
        BadUmi = 0
		#Looping through each line and setting a variable for the first read.
        for line in fh:
            test = []
			#If a line starts with an @, it is immediatly written to the new file.
            if line.startswith('@'):
                fh2.write(line)
            #If not, each line is stripped of the new line character and split into its indvidiual columns. Following that, the UMI is deduced for each read 
			#using regular expressions, and variables for the strand, chromosome, and starting position are assigned based on their column values for the read
			#Finally, the test variable is set to equal the strand, chromosome, and starting position as a tuple.
            else:
                line2 = line.strip().split()
                umi = re.search('[A-Z]+$', line2[0])
                umi2 = umi.group(0)
                strand = bit_checker(int(line2[1]))
                chr = line2[2]
                startingposition = line2[3]
                test = [strand, chr, startingposition]
                #Using the UmiChecker function, each read is checked to see if the UMI is part of the provided list, if not, a tally is given to the BadUmi counter
                if UmiChecker(umi2) == False:
                    BadUmi += 1
                else:
					#Next the script checks to see if the current combination of UMI (set as the dictionary key), and strand, chromosome, and starting position (set as the dictionary value)
					#Is in in the dictionary already. If it is, tally one for the PCR duplicate counter.
                    if (umi2,test) in testdict.items():
                        PCRDuplicates += 1
                    else:
						#If the current UMI, strand, chromosome, starting position is not in the dictionary check to see if it is soft clipped, if not, add it to the dictionary and write it to the output file.
                        test2 = []
						#If it is soft clipped, adjust the starting position and then recheck the dictionary to see if it is in there. If it is, tally one for PCR duplicates
						#If it is not, add it to the dictionary and write it to the file..
                        if re.search('^\d+[S]', line2[5]):
                            adjust = line2[5].split('S')
                            adjust = adjust[0]
                            newstartingposition = softclip_adjustment(line2[3], adjust)
                            test2 = [strand, chr, str(newstartingposition)]
                            if (umi2, test2) in testdict.items():
                                PCRDuplicates += 1
                            else:
                                fh2.write(line)
                                NotPCRDuplicates += 1
                                testdict[umi2] = [strand, chr, line2[3]]
                        else:
                            fh2.write(line)
                            NotPCRDuplicates += 1
                            testdict[umi2] = [strand, chr, startingposition]
            
#Finally print some statistics regarding the amount of PCR Duplicates, and bad UMIs were present.
print("Number of PCR Duplicates " + str(PCRDuplicates))
print("Number of Incorrect UMIs " + str(BadUmi))