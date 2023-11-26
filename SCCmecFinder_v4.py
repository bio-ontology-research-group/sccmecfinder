# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 12:00:30 2016

@author: Hulya Kaya
"""

#!usr/bin/python

##############################################################
# IMPORT OF LIB AND START OF TIME
##############################################################
import time
import subprocess
import os, sys
import glob
from argparse import ArgumentParser
from sets import Set
t0 = time.time()


definition_SCCmec = {'SCCmec_type_I(1B)': set(['ccrA1', 'ccrB1', 'mecA', 'dmecR1', 'IS1272']), \
		     'SCCmec_type_II(2A)': set(['ccrA2', 'ccrB2', 'mecA', 'mecR1', 'mecI']), \
		     'SCCmec_type_III(3A)': set(['ccrA3', 'ccrB3', 'mecA', 'mecR1', 'mecI']), \
		     'SCCmec_type_IV(2B)': set(['ccrA2', 'ccrB2', 'mecA', 'dmecR1', 'IS1272']), \
		     'SCCmec_type_IV(2B&5)': set(['ccrA2', 'ccrB2', 'ccrC11','mecA', 'dmecR1', 'IS1272']), \
		     'SCCmec_type_V(5C2)': set(['ccrC11', 'mec-class-C2', 'mecA']), \
		     'SCCmec_type_V(5C2&5)': set(['ccrC11', 'ccrC12', 'mec-class-C2', 'mecA']), \
		     'SCCmec_type_VI(4B)': set(['ccrA4', 'ccrB4', 'mecA', 'dmecR1', 'IS1272']), \
		     'SCCmec_type_VII(5C1)': set(['ccrC11', 'mec-class-C1', 'mecA']), \
		     'SCCmec_type_VIII(4A)': set(['ccrA4', 'ccrB4','mecA', 'mecR1', 'mecI']), \
		     'SCCmec_type_IX(1C2)': set(['ccrA1', 'ccrB1' ,'mec-class-C2', 'mecA']), \
		     'SCCmec_type_X(7C1)': set(['ccrA1', 'ccrB6', 'mec-class-C1', 'mecA']), \
		     'SCCmec_type_XI(8E)': set(['ccrA1', 'ccrB3', 'mecC', 'mecR1', 'mecI']), \
		     'SCCmec_type_XII(9C2)': set(['ccrC21', 'mec-class-C2', 'mecA']), \
		     'SCCmec_type_XIII(9A)': set(['ccrC21', 'mecA', 'dmecR1', 'IS1272'])}

definition_SCCmec_classes = {'SCCmec_type_I(1B)': set(['ccr class 1', 'mec class B']), \
			'SCCmec_type_II(2A)': set(['ccr class 2', 'mec class A']), \
			'SCCmec_type_III(3A)': set(['ccr class 3', 'mec class A']), \
			'SCCmec_type_IV(2B)': set(['ccr class 2', 'mec class B']), \
			'SCCmec_type_IV(2B&5)': set(['ccr class 5', 'ccr class 2', 'mec class B']), \
			'SCCmec_type_V(5C2)': set(['ccr class 5', 'mec class C2']), \
			'SCCmec_type_V(5C2&5)': set(['ccr class 5&5', 'mec class C2']), \
			'SCCmec_type_VI(4B)': set(['ccr class 4', 'mec class B']), \
			'SCCmec_type_VII(5C1)': set(['ccr class 5', 'mec class C1']), \
			'SCCmec_type_VIII(4A)': set(['ccr class 4', 'mec class A']), \
			'SCCmec_type_IX(1C2)': set(['ccr class 1', 'mec class C2']), \
			'SCCmec_type_X(7C1)': set(['ccr class 7', 'mec class C1']), \
			'SCCmec_type_XI(8E)': set(['ccr class 8', 'mec class E']), \
			'SCCmec_type_XII(9C2)': set(['ccr class 9', 'mec class C2']), \
			'SCCmec_type_XIII(9A)': set(['ccr class 9', 'mec class A'])}

def performing_ccr_gene_complex_typing (ccrAB_genes, ccrC_genes):
	if all([gene in ccrAB_genes for gene in ["ccrA1", "ccrB1"]]) == True: classes.append("ccr class 1")
	if all([gene in ccrAB_genes for gene in ["ccrA2", "ccrB2"]]) == True: classes.append("ccr class 2")
	if all([gene in ccrC_genes for gene in ["ccrC11", "ccrC12"]]) == True: classes.append("ccr class 5&5")
	elif all([gene in ccrC_genes for gene in ["ccrC11"]]) == True: classes.append("ccr class 5")
	if all([gene in ccrAB_genes for gene in ["ccrA3", "ccrB3"]]) == True: classes.append("ccr class 3")
	if all([gene in ccrAB_genes for gene in ["ccrA4", "ccrB4"]]) == True: classes.append("ccr class 4")
	if all([gene in ccrAB_genes for gene in ["ccrA5", "ccrB3"]]) == True: classes.append("ccr class 6")
	if all([gene in ccrAB_genes for gene in ["ccrA1", "ccrB6"]]) == True: classes.append("ccr class 7")
	if all([gene in ccrAB_genes for gene in ["ccrA1", "ccrB3"]]) == True: classes.append("ccr class 8")
	if all([gene in ccrC_genes for gene in ["ccrC21"]]) == True: classes.append("ccr class 9")
	return classes

def performing_mec_gene_complex_typing (mec_genes, classes):
	if all ( [gene in mec_genes for gene in ["mecA", "mecR1", "mecI"]]) == True: classes.append("mec class A")
	if all ( [gene in mec_genes for gene in ["mecA", "dmecR1", "IS1272"]]) == True: classes.append("mec class B")
	if all ( [gene in mec_genes for gene in ["mecA", "mec-class-C1"]]) == True: classes.append("mec class C1")
	if all ( [gene in mec_genes for gene in ["mecA", "mec-class-C2"]]) == True: classes.append("mec class C2")    
	if all ( [gene in mec_genes for gene in ['mecC', 'mecR1', 'mecI']]) == True: classes.append("mec class E") 
	return classes

def performing_SCCmectyping (classes):
	if all([gene in classes for gene in ["ccr class 1", "mec class B"]]) == True: SCCmectyping.append("SCCmec_type_I(1B)")
	if all([gene in classes for gene in ["ccr class 2", "mec class A"]]) == True: SCCmectyping.append("SCCmec_type_II(2A)")
	if all([gene in classes for gene in ["ccr class 3", "mec class A"]]) == True: SCCmectyping.append("SCCmec_type_III(3A)")
	if all([gene in classes for gene in ["ccr class 2", "ccr class 5", "mec class B"]]) == True: SCCmectyping.append("SCCmec_type_IV(2B&5)")
	elif all([gene in classes for gene in ["ccr class 2", "mec class B"]]) == True: SCCmectyping.append("SCCmec_type_IV(2B)")
	if all([gene in classes for gene in ["ccr class 5", "mec class C2"]]) == True: SCCmectyping.append("SCCmec_type_V(5C2)")
	if all([gene in classes for gene in ["ccr class 5&5", "mec class C2"]]) == True: SCCmectyping.append("SCCmec_type_V(5C2&5)")    
	if all([gene in classes for gene in ["ccr class 4", "mec class B"]]) == True: SCCmectyping.append("SCCmec_type_VI(4B)")
	if all([gene in classes for gene in ["ccr class 5", "mec class C1"]]) == True: SCCmectyping.append("SCCmec_type_VII(5C1)")
	if all([gene in classes for gene in ["ccr class 4", "mec class A"]]) == True: SCCmectyping.append("SCCmec_type_VIII(4A)")
	if all([gene in classes for gene in ["ccr class 1", "mec class C2"]]) == True: SCCmectyping.append("SCCmec_type_IX(1C2)")
	if all([gene in classes for gene in ["ccr class 7", "mec class C1"]]) == True: SCCmectyping.append("SCCmec_type_X(7C1)")
	if all([gene in classes for gene in ["ccr class 8", "mec class E"]]) == True: SCCmectyping.append("SCCmec_type_XI(8E)")
	if all([gene in classes for gene in ["ccr class 9", "mec class C2"]]) == True: SCCmectyping.append("SCCmec_type_XII(9C2)")
	if all([gene in classes for gene in ["ccr class 9", "mec class A"]]) == True: SCCmectyping.append("SCCmec_type_XIII(9A)")
	return classes


##############################################################
#   COMMAND LINE OPTIONS (argparse)
##############################################################
parser = ArgumentParser(description='Prediction of SCCmec cassette in S. aureus')
parser.add_argument("-iDb", type=str, dest="fasta_file_db", default="", help="Fasta input file for MyDbFinder", required=True)
parser.add_argument("-iKm", type=str, dest="fasta_file_km", default="", help="Fasta input file for MyKmerFinder", required=True)
parser.add_argument("-k", type=str, dest="ID_threshold", default="", help="min. %ID threshold", required=True)
parser.add_argument("-l", type=str, dest="len_threshold", default="", help="min. % length threshold", required=True)
parser.add_argument("-o", type=str, dest="output_file", default="", help="Output file name", required=True)
parser.add_argument("-d", type=str, dest="output_dir", default="", help="Output directory name", required=True)
parser.add_argument("-db_dir", type=str, dest="database_dir", default="", help="Database directory name", required=True)
parser.add_argument("-sc_dir", type=str, dest="script_dir", default="", help="Script directory name", required=True)
parser.add_argument("-db_choice", type=str, dest="database_choice", default="", help="Database choice - reference or extended", required=True)
args = parser.parse_args()
print 'Finding the SCCmec cassette of:', args.fasta_file_db
	

##############################################################
#   INITIALIZATION
##############################################################
hit = {}
(classes, SCCmectyping, ccrAB_genes, ccrC_genes, mec_genes) = ([], [], [],[], [])
(total_gene_list, listofgenes_mec_set, missing, additional_complex, additional_gene, ccrC_genes_2, classes_set) = (set(), set(), set(), set(), set(), set(), set())
(MRSA, subtyping_myDb, subtyping_myKm) = ('',[],'')
perform_subtyping = 'No'

##############################################################
#   MAKING SPECIFIEC DIR IF NOT EXISTING
##############################################################
if not os.path.exists(args.output_dir):
	os.mkdir(args.output_dir)
mypath = args.output_dir + "/"

##############################################################
#   RUNNING MYDBFINDER AND MYKMERFINDER
#   both databases is harcoded
##############################################################

print 'Running the BLAST-based approach'
subprocess.check_output(["perl", args.script_dir +"/"+"CGE_MyDbFinder-1.1.pl", "-o", mypath, "-r", args.database_dir+"/"+"/single_genes_database_20171117.fasta", "-k", args.ID_threshold, "-l", args.len_threshold, "-i", args.fasta_file_db])

print 'Running the kmer-based approach'
if args.database_choice == 'reference':
	subprocess.check_output(["python2.7", args.script_dir +"/"+"findtemplate.py", "-i", args.fasta_file_km, "-t", args.database_dir+"/template_db/MyKmerFinder_reference_template", "-o", args.output_dir+"/"+"results_MyKmerFinder.txt"])
elif args.database_choice == 'extended':
	subprocess.check_output(["python2.7", args.script_dir +"/"+"findtemplate.py", "-i", args.fasta_file_km, "-t", args.database_dir+"/template_db/MyKmerFinder_extended_template", "-o", args.output_dir+"/"+"results_MyKmerFinder.txt"])
print "Done"

##############################################################
#   EXTRACTION OF DATA
##############################################################
#The BLAST-based approach
num = 1
try:
	fh = open(mypath+"results_tab_MyDbFinder.txt")
	first_line_MyDbFinder = fh.readline()
	line = fh.readline()
	num = 1
	while line:
		fields = line.split('\t')
		if len(fields) == 4:
			print 'Not a correct tab-separeted file.\nStopped execution of SCCmecFinder'
			exit()
		#Separating the genes based on gene complex - either ccr, mec gene complex or subtyping target
		if fields[0].split(':')[0][:3] == 'ccr':
	    		if (fields[0].split(':')[0][:4] == 'ccrA' or fields[0].split(':')[0][:4] == 'ccrB'):
				ccrAB_genes.append(fields[0].split(':')[0])
	    		if fields[0].split(':')[0][:4] == 'ccrC':
				ccrC_genes.append(str(fields[0].split("|")[0].split("-")[0])+str(num))
				num += 1
		elif (fields[0].split(':')[0][:3] == 'mec' or fields[0].split(':')[0][:2] == 'IS' or fields[0].split(':')[0][:3] == 'dme'):
			mec_genes.append(fields[0].split(':')[0])
			if fields[0].split(':')[0][:10] == 'mecALGA251': MRSA = fields[0].split(':')[0][:10]
			elif fields[0].split(':')[0][:4] == 'mecA': MRSA = fields[0].split(':')[0][:4]
			elif fields[0].split(':')[0][:9] == 'mec-class': MRSA = 'mecA'
		elif fields[0].split(':')[0][:3] == 'sub': subtyping_myDb.append(fields[0].split(':')[0].split("|")[0])
		line = fh.readline()
	fh.close()
except IOError, (errno, errtext):
	if errno==2:
		print 'Error! No such MyDbFinder result file\nStopped execution of SCCmecFinder'
		exit()	   
		

#The kmer-based approach
try:
	fh = open(mypath+"results_MyKmerFinder.txt")
	first_line_MyKmerFinder = fh.readline()
	line = fh.readline()
	while line:
		fields = line.split('\t')
		if len(fields) > 10:
			print 'Not a correct tab-separeted file.\nStopped execution of SCCmecFinder'
			exit()
		if float(fields[6].strip()) >= 50.00: #CUT-OFF VALUE - adjust if needed
			hit[fields[0]] = {'Score': fields[1], 'Expected': fields[2], 'z':fields[3], 'p_value': fields[4], 'query coverage [%]': fields[5], 'template coverage [%]': fields[6], 'depth': fields[7], 'Kmers in Template': fields[8], 'Description': fields[9]}
		line = fh.readline()
	fh.close()
except IOError, (errno, errtext):
	if errno == 2:
		print 'Error! MyKmerFinder did not produce a result file\nStopped execution of SCCmecFinder'

#Merging all genes together - needed to see if there is additional complex/genes

[total_gene_list.add(x) for x in ccrAB_genes]
[total_gene_list.add(x) for x in ccrC_genes]
[total_gene_list.add(x) for x in mec_genes]

##############################################################
#   PERFORMING TYPING BASED ON SINGLE GENES
##############################################################
performing_ccr_gene_complex_typing(ccrAB_genes, ccrC_genes)
performing_mec_gene_complex_typing (mec_genes, classes)
performing_SCCmectyping (classes)

[classes_set.add(x) for x in classes] #Converting classes from a list to a set

##############################################################
#   SCCMECTYPING BASED ON SINGLE GENES - the mec error - running BLAST-based approach again with a new database
##############################################################
redo = 'false'
if not SCCmectyping:
	for x in classes:
		if (x == 'ccr class 9' or x == 'ccr class 5&5' or x == 'ccr class 5' or x == 'ccr class 7' or x == 'ccr class 1'): redo = 'true'
	for x in classes:
		if x[:3] == 'mec': redo = 'false'
	if 'IS1272' in total_gene_list: redo = 'false'
	if 'mecI' in total_gene_list: redo = 'false'
if redo == 'true':
	for x in os.listdir(mypath):
		if not x == "results_MyKmerFinder.txt": os.rename(mypath+x, mypath+"old"+x)
	subprocess.check_output(["perl", args.script_dir +"/"+"CGE_MyDbFinder-1.1.pl", "-o", mypath, "-r", args.database_dir +"/"+"mec_database_20171117.fasta", "-k", args.ID_threshold, "-l", args.len_threshold, "-i", args.fasta_file_db])
	listofgenes_mec = []
	mec_genes_2 = []
	try:
		fh = open(mypath+'results_tab_MyDbFinder.txt', "r")
		line = fh.readline()
		while line:
			fields = line.split('\t')
			if len(fields) == 4:
				print 'Not a correct tab-separeted file.\nStopped execution of SCCmecFinder'
				exit()
			gene = fields[0].split(':')[0] + ', %ID ' + fields[1] + ', %cov ' + fields[2]
			if fields[0].split(':')[0][:3] == 'mec':
				listofgenes_mec.append(gene.split(",")[0])
				MRSA = 'mecA'
			line = fh.readline()
		fh.close()
	except IOError, (errno, errtext):
		if errno==2:
			print 'Error! No such MyDbFinder result file 2\nStopped execution of SCCmecFinder'
			exit()
	try:
		fh = open(mypath+"oldresults_tab_MyDbFinder.txt", "r")
		f = open(mypath+"results_tab2.txt", 'w') #temp. file		
		first_line_MyDbFinder = fh.readline()
		f.write(first_line_MyDbFinder.rstrip('\r\n') + "\n")
		line = fh.readline()
		while line:
			fields = line.split('\t')
			if len(fields) == 4:
				print 'Not a correct tab-separeted file.\nStopped execution of SCCmecFinder'
				exit()	 
			if not line[:4] == 'mecA': f.write(line.rstrip('\r\n') + "\n")
			if line[:4] == 'mecA':
				fh_mec = open(mypath+"results_tab_MyDbFinder.txt", "r")
				first = fh_mec.readline()
				line_mec = fh_mec.readline()
				f.write(line_mec.rstrip('\r\n') + "\n")
			line = fh.readline()
		fh.close()
	except IOError, (errno, errtext):
		if errno==2:
			print 'Error! No such MyDbFinder result file 3\nStopped execution of SCCmecFinder'
			exit()
	#Removing the unneeded files
	[os.remove(x) for x in glob.glob(mypath+"old*")]
	os.remove(mypath+"results_tab_MyDbFinder.txt")
	os.rename(mypath+"results_tab2.txt", mypath+"results_tab_MyDbFinder.txt")
	for x in mec_genes:
		if not x == 'mecA' or x == 'mecC': total_gene_list.discard(x)
		if x == 'mecA' or x == 'mecC': listofgenes_mec.append(x)			
	total_gene_list = total_gene_list | listofgenes_mec_set #Merging all genes
	performing_mec_gene_complex_typing (listofgenes_mec, classes) #Redoing the mec gene complex determination
	performing_SCCmectyping (classes) #Redoing the SCCmec typing


##############################################################
#   checking whether to perform subtyping also 
##############################################################
if len(SCCmectyping) == 1 and (SCCmectyping[0] == 'SCCmec_type_IV(2B)' or SCCmectyping[0] == 'SCCmec_type_V(5C2&5)' or SCCmectyping[0] == 'SCCmec_type_V(5C2)'): perform_subtyping = 'Yes'

##############################################################
#   PERFORMING SCCMECTYPING BASED ON WHOLE CASSETTE
#   checking whether to perform subtyping also 
##############################################################

hit_sorted = sorted(hit.items(), key=lambda x : x[1]['Score'], reverse=True) #hit_sorted is now a list and sorted based on the score
accession = ''
if len(hit_sorted) > 0:
	MyKmerFinder_hit = hit_sorted[0][0].split("|")[0] #Getting the best hit
	if not (MyKmerFinder_hit == 'SCCmec_type_IV(2B)' or MyKmerFinder_hit == 'SCCmec_type_V(5C2&5)' or MyKmerFinder_hit == 'SCCmec_type_V(5C2)'): perform_subtyping = 'No'
	else: subtyping_myKm = hit_sorted[0][0].split("|")[1] #Getting the subype of the best hit
	accession = hit_sorted[0][0].split("|")[-1] #Getting the accesionnumber
else:
	MyKmerFinder_hit = ''
	
##############################################################
#	COMPARING THE TYPINGS    
#	PRODUCING SCCMECFINDER RESULT FILE
##############################################################

f = open(mypath+args.output_file+".txt", 'w')  #Output file - result file

flag = 0
if MRSA:
	f.write('The input organism was prediced as a MRSA isolate\nThe ' + MRSA + ' gene was detected\n')
else:
	f.write('The input organism was prediced as a MSSA isolate\nThe mecA/mecC gene was not detected\n')

if SCCmectyping:
	if len(SCCmectyping) >= 2: #This will happen when the combination of gene complex can result into mulitple SCCmec element
		numberofSCC = str(len(SCCmectyping))
		f.write("Alert! Possible " + numberofSCC + " SCCmec cassettes were predicted based on detection of genes (see list). Look also the list of homology to whole cassettes below and the predicted gene complexes:\n")
		[f.write(x + "\n") for x in SCCmectyping]
		f.write("Following gene complexes based on prediction of genes was detected :\n")
		[f.write(x + "\n") for x in classes]	
	if len(SCCmectyping) == 1:
		MyDbFinder_hit = SCCmectyping[0] 
		if MyKmerFinder_hit == MyDbFinder_hit:
			if perform_subtyping == 'Yes':
				if len(subtyping_myDb) > 1:
					flag = 1
					f.write("Alert! Multiple subtype target genes detected.\n")
				if len(subtyping_myDb) == 1:
					MyDbFinder_hit = MyDbFinder_hit.split("(")[0].split("_")[0]+"_"+MyDbFinder_hit.split("(")[0].split("_")[1]+"_"+subtyping_myDb[0].split("-")[1]
					if MyDbFinder_hit != subtyping_myKm:
						flag = 1
						f.write("Alert! Contradicting subtype predictions.\n")
						MyKmerFinder_hit = subtyping_myKm
					if MyDbFinder_hit == subtyping_myKm: 
						MyKmerFinder_hit = subtyping_myKm
				if len(subtyping_myDb) == 0: 
					f.write("Alert! non-subtypeable\n")
			if flag == 0: f.write("One SCCmec element detected.\n")
			f.write("Prediction based on genes:\nPredicted SCCmec element: %s\n" % (MyDbFinder_hit))
			f.write("Prediction based on homology to whole cassette:\nPredicted whole cassette and %template coverage: " + MyKmerFinder_hit + " " + hit_sorted[0][1]['template coverage [%]'].strip() + "%\n")
### ADDTIONAL COMPPLEX/GENES
			#Finding additional complex and/or genes
			for x, y in definition_SCCmec_classes.items():
				if x == MyDbFinder_hit: additional_complex = classes_set.difference(y)
			for x, y in definition_SCCmec.items():
				if x == MyDbFinder_hit: additional_gene = total_gene_list.difference(y)
			if additional_complex:
				f.write("Be aware! Additional complex(es) was found.\n")
				[f.write(x + "\n") for x in additional_complex]
			elif additional_gene:
				f.write("Be aware! Additional gene(s) was found.\n")
				for x in additional_gene:
					if x[:5] == 'ccrC1': f.write(x[:5] + "\n")
					else: f.write(x + "\n")				
### CONTRADICTING TYPE-PREDICTIONS
		elif MyKmerFinder_hit != MyDbFinder_hit: 
			if MyKmerFinder_hit:
				f.write('Alert! ' +'Contradicting predictions\nPrediction based on genes:\nPredicted SCCmec element: %s\n' % (MyDbFinder_hit))
				f.write('Prediction based on homology to whole cassette:\nPredicted whole cassette and % template coverage: ' + MyKmerFinder_hit + " " + hit_sorted[0][1]['template coverage [%]'].strip() + "%\n")
				for x, y in definition_SCCmec_classes.items():
					if x == MyDbFinder_hit: additional_complex = classes_set.difference(y)
				if additional_complex:
					f.write("Be aware! Additional complex(es) was found.\n")
					[f.write(x + "\n") for x in additional_complex]
				elif additional_gene:
					f.write("Be aware! Additional gene(s) was found.\n")
					for x in additional_gene:
						if x[:5] == 'ccrC1': f.write(x[:5] + "\n")
						else: f.write(x + "\n")
### NO CASSETTE > 40%
			elif MyKmerFinder_hit == '':
				f.write("Alert! It was not possible to predict any whole SCCmec cassette with an template coverage of at least 40%\n")
				f.write("However it was possible to predicted a SCCmec cassette based on genes prediction:\n%s with following gene complexes: %s and %s\n" % (MyDbFinder_hit, classes[0], classes[1]))
elif not SCCmectyping:
	if (MRSA == 'MSSA' and total_gene_list and not MyKmerFinder_hit): f.write('Alert! The isolate migth have a pseudo element.\n')
	else: f.write("No SCCmec element was detected\n")
	f.write("Prediction based on genes:\nPredicted SCCmec element: none\nPrediction based on homology to whole cassette:\nPredicted whole cassette and %template coverage: none\n")
	if len(classes) != 0: #Checking if any classes was found. 
		f.write("\nHowever it was possible to identify one or more gene complex(es) based on prediction of genes:\n")
		[f.write(x + "\n") for x in classes]
	if MyKmerFinder_hit:
		f.write('\nAssuming that the cassette with the best homology, which is ' + hit_sorted[0][0] + ' with ' + hit_sorted[0][1]['template coverage [%]'].strip() + '% template coverage, is the best hit then following genes/complexes are missing:\n')
		for x, y in definition_SCCmec.items():
			if x == MyKmerFinder_hit: missing = y.difference(total_gene_list)	
		if missing:
			for x in missing:
				if x[:5] == 'ccrC1': f.write(x[:5] + "\n")
				else: f.write(x + "\n")
f.close()

##############################################################
#   PRINTING OUTPUT TO HTML
##############################################################
HTML_MyKmerFinder = open(mypath+"HTML_MyKmerFinder.txt", 'w')
HTML_MyKmerFinder.write("SCCmec elements\n")
if hit_sorted:
	HTML_MyKmerFinder.write(first_line_MyKmerFinder)
	if len(hit_sorted) >= 5: length = 5
	elif len(hit_sorted) < 5: length = len(hit_sorted)
        x = 0
        while x < length:
		HTML_MyKmerFinder.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                 (hit_sorted[x][0], hit_sorted[x][1]['Score'], hit_sorted[x][1]['Expected'], hit_sorted[x][1]['z'], hit_sorted[x][1]['p_value'], hit_sorted[x][1]['query coverage [%]'], hit_sorted[x][1]['template coverage [%]'] , hit_sorted[x][1]['depth'], hit_sorted[x][1]['Kmers in Template'], hit_sorted[x][1]['Description']))
	        x = x+ 1
elif not hit_sorted:
        	HTML_MyKmerFinder.write('no whole SCCmec cassette was found\n')
		
##############################################################
#   REMOVING FILES
##############################################################

[os.remove(x) for x in glob.glob(mypath+"MyKmerFinder_template.*")]
os.remove(mypath+"Hit_in_genome_seq.fsa")
os.remove(mypath+"results.txt")
os.remove(mypath+"Database_gene_seq.fsa")

print 'Prediction completed'
t1 = time.time()
print "Time used: %s sec" % (int(t1 - t0))
