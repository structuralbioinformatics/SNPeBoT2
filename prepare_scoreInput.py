import sys
import re
TF = sys.argv[1]
thread = sys.argv[2]
jobid = sys.argv[3]

pattern = r'([ATGC]+\.txt)$'

threadsprefix = re.sub(pattern, '', thread)
threadsprefix = threadsprefix.split("/")[-1]

fd = open(thread,'r')
executefile = fd.read()
fd.close()

#executefile = executefile.replace("\n","\\\\n")

#pdb_file = re.sub(r">dna[\s\S]+", ">dna\\\\n'+dna+';0\\\\n//\\\\n>dna_fixed\\\\n'+dna+';0\\\\n//\\\\n", executefile)
pdb_file = re.sub(
    r">dna[\s\S]+", 
    ">dna\\n+dna+;0\\n//\\n>dna_fixed\\n+dna+;0\\n//\\n", 
    executefile
)

newFile = f"""\nwith open('../../../Output_{jobid}/{TF}_seq_table_{jobid}.txt', 'r') as fd:
    for line in fd:
        dna = line.strip().upper()

        file_name = '{threadsprefix}'+dna+'.txt'

        newfile = '''{pdb_file}'''.replace('+dna+', dna)

        wd = open(file_name, 'w')
        wd.write(newfile)
        wd.close()
"""

#newFile = "\nwith open('../NNinput/class_and_seqs/"+TF+"_seq_table.txt','r') as fd:\n    for line in fd:\n        dna = line.strip().upper()\n\n        file_name = '"+threadsprefix+"'+dna+'.txt'\n\n        newfile = '"+pdb_file+"' \n\n        wd = open(file_name,'w')\n        wd.write(newfile)\n        wd.close()\n"

wd = open("create_threads.py","w")
wd.write(newFile)
wd.close()



