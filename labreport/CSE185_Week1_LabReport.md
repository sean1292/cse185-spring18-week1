# *E.coli* Sample Sequence Analysis for Determening Possible Antibiotic Resistance 
##### Sean Sadykoff

## Abstract
In no more than 100 words, briefly summarize what was done in the lab this week, what the findings were, and why they were important.

## Introduction

Seasonal influenza virus is estimated to account for 50000 yearly deaths in just the United States alone, and even though rapid advancements in medicine are made every day it remains a hot issue in medicine today (1). Flu virus’s genome consist of only approximately 13500 base pairs of negative-sense, single-stranded RNA that is divided into 8 segments containing 11 total genes (2). Due to high mutation rate, (1 mutation/genome/replication) after infecting the host, flu virus quickly mutates and recombines to potentially produce a population of quasispecies that can in turn produce a phenotype that is resistant to certain antibody or medication (2).  One of the main mechanisms for Influenzas’ infection of host cells involves a glycoprotein called hemagglutinin, which allows the virus to bind sialic acid that is present on the cell membranes in upper respiratory tract or red blood cells (2).  Using an antibody for that glycoprotein is common technique that helps prevent binding of the virus to the host cell and therefore stop the spread of the infection (2). However, previously mentioned fast mutation rate may cause a mutation to occur the hemagglutinin glycoprotein that would in turn make the antibody unable to bind and therefore make the treatment ineffective.   

Current technology allowed to significantly reduce the cost of the whole-genome sequencing as well as the time it takes to do so and by taking advantage of this technology doctors and scientists might be able to treat ill patient much faster and more efficient (3). Whenever doctor is faced with decision of properly determining appropriate medication to treat a given patient he/she may be able to take advantage of the whole-genome sequencing in order to determine the variations of flu virus that infected the patient. By being able to determine the types of mutations present in the population, doctors can in turn be able to better decide of what combination of antibiotics/anti-bodies to use in order to be able to efficiently battle the infection that they are faced with.   

In this lab we mapped paired-end sequencing data for E.coli  to a reference genome in order to determine the locations of possible genetic variants. By analyzing genetic variation present in our population, we may be able to explain certain antibiotic and antibody resistance as well as make predictions on possible treatment options. Using this method may allow to narrow down the choice of medication/combination of medications based on the specific characteristics present in the virus’s genetic code.


## Methods
This section should contain sufficient information so that other bioinformaticists could reproduce your results. You should briefly describe your raw data (what is it, what is the name of the reference) and describe what you did with it. You should write this in 2-3 paragraphs, not in a list. When you use a bioinformatics software program, do not write out the full command you typed, but do specify which program (ie ‘bwa-mem’ or ‘samtools tview’) you used and whether you used the default options. If you did not use the defaults, you should specify the exact settings you used. The first time you mention bioinformatics software or an online tool, you should cite it and specify which version of the tool you used. The correct citation for most software can be found by looking up its documentation on line (you don’t have to cite common tools like python or perl or the bash shell). If you write a custom script, (for example, our awk script from week 1), include that code in the labreport folder and reference it in your writeup.

## Results
This section should include the results of your data processing and data analysis, and may include tables with read lengths, pictures of quality distributions, or tables of gene names for examples. In the text, briefly restate how you got the results in full sentences, but in less detail than the methods, before you say what the results are (ie ‘reads were mapped to the reference and scanned to identify positions that likely contained mutations. We found….’). Refer to tables and figures by number, and include a brief descriptive title for each. Be sure to include any results specifically requested in the lab project tutorial. The results section should be as objective as possible, so please refrain from interpreting the meaning or significance here. It should be just the facts.

## Discussion
In 2-3 paragraphs, explain what you think the results mean, and why you are interpreting them this way. If you encountered any problems, or answered questions, discuss them and suggest ways to solve them with future experiments or analyses. Also include any information specifically requested in the tutorial.

## Citations
1. Palese, Peter. (2017). Influenza: A broadly protective antibody. Nature, 551(7680), pp.310-311.   
2. Gymrek, Melissa Ann. "CSE185_Spring2018_Week2_Monday" 10 Apr. 2018. CSE 185, UCSD. Microsoft PowerPoint presentation.
3. Peacock, Sharon. (2014). Health care: Bring microbial sequencing to hospitals, 509, pp.557–559.    
You can use any commonly used format you like, but be consistent. Lab reports will be submitted via turnitin to check for plagiarism, so be sure to cite other people’s ideas, and put everything in your own words (paraphrasing) if you aren’t using direct quotes.
