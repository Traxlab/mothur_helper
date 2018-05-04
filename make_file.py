import csv
import os
import sys
#os.system("rm names.txt")
os.system("ls *.fastq >> names.txt")
with open("names.txt","r") as namefile:
	names = namefile.readlines()
	names = [x.strip() for x in names]
with open("fastq_info_batch","w") as info_file:
	for name in names:
		string = "fastq.info(fastq=" + name +")\n" 
		info_file.write(string)
sys.exit()
os.system("mothur fastq_info_batch")
sys.exit()
with open('stability.files', 'w') as csvfile:
	writer = csv.writer(csvfile,delimiter='	')
	for i in range(0,len(names)-1,2): 
		writer.writerow([names[i].split("_")[0],names[i],names[i+1]])
os.system("rm names.txt")
#os.system("mothur batchfile")
