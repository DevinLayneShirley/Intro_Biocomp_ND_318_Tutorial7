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
#check out summary table
seqinfo.head()

#FIRST GRAPH: HISTOGRAM OF SEQUENCE LENGTHS
seqlen=ggplot(seqinfo,aes('sequenceLength'))
seqlen+geom_histogram(binwidth=10,fill='cadetblue',color='black')+theme_classic()+ggtitle("Sequence Lengths")
#https://i.stack.imgur.com/fMx2j.png gives a bunch of colors in python,
#I chose to make the bars cadetblue with black outlines

#SECOND GRAPH: HISTOGRAM OF GC CONTENT 
gccont=ggplot(seqinfo,aes('percentGC'))
gccont+geom_histogram(binwidth=1,fill='mediumorchid',color='black')+theme_classic()+ggtitle("GC Content of the Sequences")

#2 Create graph relating two variables
#add data
geyser=pandas.read_csv("old_faithful_erruptions.txt",sep="\t",header=0)
#check to make sure data loaded properly 
geyser.head()
#create graph with treadline
ggraph=ggplot(geyser,aes('eruptions','waiting'))
ggraph+geom_point(color='firebrick')+xlab("Eruption Duration (min)")+ylab("Time Waited (min)")+stat_smooth(method="lm")#stat_smooth adds the trend line

#3 two figures that summarize data in data.txt
#load data
data=pandas.read_csv("data.txt",header=0)
#check to make sure data loaded properly
data.head()
#create barplot showing means of the four populations
means=ggplot(data)+theme_classic()+xlab("Populations")+ylab("Mean Number of Observations")
means+geom_bar(aes(x="factor(region)",y="observations",fill="region"),stat="summary",fun_y=numpy.mean)
#calculate means to check bar plot
data.groupby(['region'])['observations'].mean()
#this was correct, they are only small differences 
#THIS ONE WORKS

#create scatterplot of all observations #THIS IS ALL WRONG
scat=ggplot(data,aes(x="observations"))+theme_classic()
scat+geom_point(aes(color="region"))+ylab("Number of Observations")

scat=ggplot(data,aes("observations"))+theme_classic()
scat+geom_dotplot(aes(color="region"))+xlab
