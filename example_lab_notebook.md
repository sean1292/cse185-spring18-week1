# Lab Notebook - Week 1
#### Name: Pikachu

## Date: 04/04/18 (Lab session 1)

### Part 2
Look at the first couple lines in my SAM file
```
head -n 5 trimpaired.sam
```

Output:
```
@SQ	SN:NC_000913.3	LN:4641652
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1039	CL:bwa mem NC_000913.3.fasta trimpair1.fastq trimpair2.fastq
SRR1363232.1	83	NC_000913.3	1087469	60	101M	=	1087258	-312	AGCCATGCGATCCAGGATGGCACGGGTTTTATCTGTTGAACCGTCATTTACGGCNNTAACTTCAATGTTCTCATAACGCTGTGCTAAAGCGGCGTGTATGG	9?DDDB?9::?DCCDCCAB@@@DB@>:@5(;>?C>CA=1=BD@AEHHHE@A@5.##GIIGF?<IIGJJIGGGGFDFJJIGEFJIEE:EGGHDHDBFDF@@@	NM:i:2	MD:Z:54A0A45	AS:i:97	XS:i:0
SRR1363232.1	163	NC_000913.3	1087258	60	101M	=	1087469	312	TTCGAATACGAGGATTACCGGTTACGGCACCCACACGCGGGTTGTACAACATCGGTTCCACAATATATGCCGCCGCATCGCGGTCTAATAACGCATCGCCA	@?@FDD?B=FAFA<BHIIBHGEGGIHHGG<GG=DGGEGGGF/=@;AD@CDBCDDBDBBDEACD@CCD>(:A@B@BDD<@BDBB>B7:>ACACBDD<99@B?	NM:i:0	MD:Z:101	AS:i:101	XS:i:0
SRR1363232.2	99	NC_000913.3	1081173	60	101M	=	1081413	341	CCCATCTTTCTACCCTGGAATAATCGTTTATATCCCTTGGCATTACCTCTCTTTGTTTACATTCCAACATCATTTTATAAACATTCCGCTTGTGTTTTTCT	CCCFFFFFHGHGHIJIJJIIJJJIIJIJIIJJGGHIGIEIJIAGGGGJJIJIGIGFB>DCHGJADHIIJGH@GDG@GHFGGGGBECCCFBDABC?CCD?CC	NM:i:0	MD:Z:101	AS:i:101	XS:i:0
```

Note: first two lines (starting with "@") are header lines. So there are only 3 reads shown here.

## Date: 04/06/18 (Lab session 2)

### Part 3
Run fastqc, which reports various quality control stats from the fastq files
```
fastqc -o . /storage/mgymrek/cse-185/public/week1/amp_res_1.fastq  /storage/mgymrek/cse-185/public/week1/amp_res_2.fastq 
```
Saved html files to Desktop week 1 folder.

### Part 4
sickle trims the ends of reads and outputs filtered fastq files to be used downstream
```
sickle pe \
    -f /storage/mgymrek/cse-185/public/week1/amp_res_1.fastq \
    -r /storage/mgymrek/cse-185/public/week1/amp_res_2.fastq \
    -t sanger \
    -o trimpair1.fastq \
    -p trimpair2.fastq \
    -s singletons.fastq
```

Rerun fastqc to see if the issues observed previously at the ends of the read are improved
```
fastqc -o . trimpair1.fastq trimpair2.fastq
```

Saved html files to Desktop week 1 folder.

#### Part 5
Running BWA MEM on our genome
```
bwa index NC_000913.3.fasta 
bwa mem NC_000913.3.fasta trimpair1.fastq trimpair2.fastq > trimpaired.sam
```

#### Part 6
Use samtools to get stats on my files and view the alignment
```
samtools flagstat trimpaired.sam # note in total includes supp alignments! accts for difference
samtools view -S -b trimpaired.sam > trimpaired.bam
samtools sort trimpaired.bam > sortedtrimpaired.bam
samtools index sortedtrimpaired.bam
samtools tview sortedtrimpaired.bam NC_000913.3.fasta
```


#### Part 7
Run Varscan to detect variants
```
samtools mpileup -f NC_000913.3.fasta sortedtrimpaired.bam > my.mpileup
java -jar VarScan.jar \
    mpileup2snp my.mpileup \
    --min-var-freq 0.7 \
    --variants \
    --output-vcf 1 > VarScan.vcf
awk '{if (NR>24) $1="Chromosome"; print}' VarScan.vcf > mymodVarScan.vcf
```
Saved mymodVarScan to Desktop week 1 folder.
