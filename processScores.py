import sys

TF=sys.argv[1]
prefix=sys.argv[2]
jobid=sys.argv[3]


classificationFile="Output_"+str(jobid)+"/"+TF+"_classification_"+str(jobid)+".csv"
ScoreLocation="Output_"+str(jobid)+"/"+"threads_"+str(jobid)+"/"+TF+"_scores/"+prefix
OutputFile="Output_"+str(jobid)+"/ModelInput_"+str(jobid)+".tsv"



print(classificationFile)
print(ScoreLocation)
print(OutputFile)
not_run = []
dictionary = {}
classDict = {}


with open(classificationFile,"r") as fd:
    for line in fd:
        # Each line in the classification file is of the format: chr5-74395388-A-G-RUNX3,GCTGACCACAA,GCTGACCGCAA,lost (ASB_ID,ref_dna,alt_dna,class) 
        line = line.split(",")
        ref = line[1]
        alt = line[2]
        #Following is a control for the cases where mutiple templated were tested
        #if len(ref.strip()) != threadLength[threader]:
        #    continue
        #ds = open(TF+"_DIMER_scores/"+ threader + ref.strip()  + ".txt","r")
        try:
            ds = open(ScoreLocation + ref.strip()  + ".txt","r")
        except:
            #print("no score file exists for "+threader + ref.strip()  + ".txt")
            not_run.append(alt.strip())
            continue
            #exit(1)
        scores = ds.read()
        try:
            refinput = scores.split("\n")[7].split("\t")[1:] 
            ds.close()
        except:
            print("failure on "+ ScoreLocation + ref.strip()  + ".txt")
            not_run.append(alt.strip())
            continue
        #ds = open(TF+"_DIMER_scores/"+ threader + alt.strip()  + ".txt","r")
        try:
            ds = open(ScoreLocation + alt.strip()  + ".txt","r")
        except:
            #print("no score file exists for "+threader + alt.strip()  + ".txt")
            not_run.append(alt.strip())#exit(1)
            continue
        scores = ds.read()
        try:
            altinput = scores.split("\n")[7].split("\t")[1:] 
            ds.close()
        except:
            print("failure on "+ ScoreLocation + ref.strip()  + ".txt")
            not_run.append(alt.strip())
            continue
        kmerInput = "\t".join(refinput)+"\t"+"\t".join(altinput)
        #the following is a replacement for multi templates (the commented that comes after)
        if line[0] in dictionary:
            dictionary[line[0]].append(kmerInput)
        else:
            dictionary[line[0]] = [kmerInput]
            classDict[line[0]] = line[3].strip()
        #if line[0]+"_"+str(len(alt.strip())) in dictionary:
        #    dictionary[line[0]+"_"+str(len(alt.strip()))].append(kmerInput)
        #else:
        #    dictionary[line[0]+"_"+str(len(alt.strip()))] = [kmerInput]
        #    classDict[line[0]+"_"+str(len(alt.strip()))] = line[3].strip()


#if len(not_run) > 0: 
    #print(not_run)
    #exit(1)
#wd = open("Homeobox_GVAT_Normalized_NNinput.tsv","a")
wd = open(OutputFile,"w")
line = 'Feature_0\tFeature_1\tFeature_2\tFeature_3\tFeature_4\tFeature_5\tFeature_6\tFeature_7\tFeature_8\tFeature_9\tFeature_10\tFeature_11\tFeature_12\tFeature_13\tFeature_14\tFeature_15\tFeature_16\tFeature_17\tFeature_18\tFeature_19\tFeature_20\tFeature_21\tFeature_22\tFeature_23\tFeature_24\tFeature_25\tFeature_26\tFeature_27\tFeature_28\tFeature_29\tFeature_30\tFeature_31\tFeature_32\tFeature_33\tFeature_34\tFeature_35\tFeature_36\tFeature_37\tFeature_38\tFeature_39\tFeature_40\tFeature_41\tFeature_42\tFeature_43\tFeature_44\tFeature_45\tFeature_46\tFeature_47\tFeature_48\tFeature_49\tFeature_50\tFeature_51\tFeature_52\tFeature_53\tFeature_54\tFeature_55\tFeature_56\tFeature_57\tFeature_58\tFeature_59\tFeature_60\tFeature_61\tFeature_62\tFeature_63\tFeature_64\tFeature_65\tFeature_66\tFeature_67\tFeature_68\tFeature_69\tFeature_70\tFeature_71\tFeature_72\tFeature_73\tFeature_74\tFeature_75\tFeature_76\tFeature_77\tFeature_78\tFeature_79\tFeature_80\tFeature_81\tFeature_82\tFeature_83\tFeature_84\tFeature_85\tFeature_86\tFeature_87\tFeature_88\tFeature_89\tFeature_90\tFeature_91\tFeature_92\tFeature_93\tFeature_94\tFeature_95\tFeature_96\tFeature_97\tFeature_98\tFeature_99\tFeature_100\tFeature_101\tFeature_102\tFeature_103\tFeature_104\tFeature_105\tFeature_106\tFeature_107\tFeature_108\tFeature_109\tFeature_110\tFeature_111\tLabel\tID\n'
wd.write(line)
for ID in dictionary:
    addition = "\t".join((dictionary[ID])) +"\t"+ classDict[ID] + "\t" + ID +"\n"
    wd.write(addition)
wd.close()


