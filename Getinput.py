import pickle 
import argparse

###
#Example usage:
#python Getinput.py --ASBInput Resources/E2F_from_loss_Test_FinalVersion.in --tf E2F4 --length 15 --jobID 1

###


parser = argparse.ArgumentParser()
parser.add_argument('--ASBInput', type=str, required=True, help='Input file with ASB data')
parser.add_argument('--tf', type=str, required=True, help='Transcription Factor name')
parser.add_argument('--length', type=int, required=True, help='Length of the sequence')
parser.add_argument('--jobID', type=str, required=True, help='Job ID for parallel processing and tracking')
args = parser.parse_args()
ASBInput = args.ASBInput
TF = args.tf
length = args.length
jobID = str(args.jobID)

def sequence_values(n):
    if n < 4:
        raise ValueError("Input must be >= 4") # Setting the lower boundary at 4 nucleotides.
    elif n <= 8:
        return [n - 1, 1]
    else:
        value1 = 7 + (n - 8) // 2
        value2 = 2 + (n - 9) // 2
        return [value1, value2]


def get_window(sequence,step):
            upstream = vals[0] - step
            downstream = vals[1] + step
            start = 25 -  upstream
            stop = 25 + downstream
            newseq = sequence[start:stop]
            return newseq.upper()


vals = sequence_values(length)

ids = {}

classesDict = {}
testIDs = set()

label = "unknown"
with open(ASBInput, "r") as fd:
    i = 1
    for raw_line in fd:
        raw_line = raw_line.strip()

        if not raw_line:
            continue  # skip empty lines

        fields = raw_line.split()  # handles tabs, spaces, mixed whitespace

        if len(fields) < 4:
            raise ValueError(
                f"Malformed input line {i}: expected ≥4 columns, got {len(fields)}\n{raw_line}"
            )

        id = fields[3]
        ids[id] = [fields[1], fields[2]]
        i += 1


#print(len(ids.keys()))

wc = open("Output_"+jobID+"/"+TF+"_classification_"+jobID+".csv","a")

positionsRun = set()
ASB = 0
dup=0


for x in ids:
      ASB+=1
      for i in range(0,8):
            if get_window(ids[x][0], i) in positionsRun:
                  dup+=1
            if get_window(ids[x][1], i) in positionsRun:
                  dup+=1
            first = get_window(ids[x][0], i)
            positionsRun.add(first)
            second = get_window(ids[x][1], i)
            positionsRun.add(second)
            addition = x + "," + first + "," + second + "," + label+ "\n"
            wc.write(addition)



#print(ASB)
#print(dup)
wc.close()
            
with open("Output_"+jobID+"/"+TF+"_seq_table_"+jobID+".txt","a") as wd:
      for x in positionsRun:
            add = x + "\n"
            wd.write(add)
#print(TF,len(positionsRun))



