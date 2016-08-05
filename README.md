# README #

##### 05.08.16
Ok made a new function that can read in paired end fastq, and align it to a ref
That function returns the SAM output which you could do other stuff to, BUT
I also made a function that converts the SAM output into a mutation Counter
Each unique combination of mutations have a count, e.g.:
(): 23000   # empty is WT
(SNP-34-A-T): 456
(INS-4-AATG, SNP-34-A-T): 203

I think it should be quite easy to go from this mutation counter to the old style summary
That function will enable the old plot types to work.

Ok made it. Now the functions in align.py outputs both new and old style.


##### 03.08.16
Trying to make it more general.
The current version only makes summaries of the mutation.
E.g. Per position
But whats needed in a current project (for Tadas) is also mutations per sequence
This requires a ...different... summing of the data
* I would still align (probably with bowtie2)
* I would then produce a dataframe containing some sort of per read-pair alignment summary
* From that summary it should be a short step to the current summary


##### Before time
I am trying to make a more general version of the software I used for the pparp project from which this code is forked.


### What is this repository for? ###

So basically it reads a bam/sam file and converts it into a dataframe of mutation frequencies along a DNA sequence