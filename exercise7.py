#set working directory to repository with data
#add packages
import pandas 
import numpy
from plotnine import *
#1 create two plots that summarize info from fasta file 
# open fasta file
fasta=open("Lecture11.fasta","r")
#create lists for storing information about sequences
sequenceID=[]
sequenceLength=[]
percentGC=[]
meltingTemp=[]
#loop through each line of fasta file to process sequences
for Line in fasta:
    # remove newline character from file line
    Line=Line.strip()
    # if a sequence record
    if '>' in Line:
        # add the sequence ID (except the ">" character) to the sequenceID list
        sequenceID.append(Line[1:])
    # if a sequence line
    else:
        # get the number of characters in the sequence and convert to a float to avoid integer division
        seqLen=float(len(Line))
        # count the number of G's and C's
        nG=Line.count("G")
        nC=Line.count("C")
        # if the sequence is 14 or fewer bases calculate melting temperature
        if seqLen<=14:
            Tm=2*(nG+nC)+2*seqLen
        else:
            Tm=-9999
        # append values to the lists
        sequenceLength.append(seqLen)
        percentGC.append((nG+nC)/seqLen*100)
        meltingTemp.append(Tm)
# combine lists into dataframe
seqinfo = pandas.DataFrame(list(zip(sequenceID,sequenceLength,percentGC,meltingTemp)),columns=['sequenceID','sequenceLength','percentGC','meltingTemp'])
# close original file
fasta.close()
#FIRST GRAPH: HISTOGRAM OF SEQUENCE LENGTHS
seqlen=ggplot(seqinfo,aes(x="sequenceLength"))
seqlen+geom_histogram()+theme_classic
#SECOND GRAPH: HISTOGRAM OF GC CONTENT 
gccont=ggplot(seqinfo,aes(x="percentGC"))
gccont+geom_histogram()+theme_classic
#2 Create graph relating two variables
#add data
geyser=pandas.read_csv("old_faithful_erruptions.csv",header=0)
#check to make sure data loaded properly 
geyser.head()


#3

