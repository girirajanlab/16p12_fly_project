#This pipeline generates a set of read counts for genes in the Drosophila transcriptome from raw FASTQ read files.
#Softwares used include: Trimmomatic v.0.36 for quality control, TopHat2 v.2.1.1 for alignment, and HTSeq v.0.6.1 for read count quantification. 
#The reads were aligned to the Drosophila Dmel6.08 reference genome and transcriptome.

#Set variables for sample name and Dmel reference genome files
runid="run_id" #Change to reflect sample name
btref="dmel6.08.fasta"
gffref="dmel6.08.embl"
gtfref="dmel6.08.embl.gtf"

#Use Trimmomatic to trim reads with poor quality, or those containing Illumina TruSeq adapters (use FastQC to determine appropriate parameters)
java -jar trimmomatic-0.36.jar SE $runid.fastq.gz $runid.trimmed.fq.gz LEADING:5 SLIDINGWINDOW:4:30 ILLUMINACLIP:truseq_19.fa:2:30:10

#Align reads to Dmel reference transcriptome
tophat -o $runid --transcriptome-index $gffref $btref $runid.trimmed.fq.gz

#Use Samtools to index output Tophat2 BAM file
cd $runid
samtools index accepted_hits.bam

#Use HTSeq-Count to quantify reads based on the Drosophila transcriptome
htseq-count -f bam -s reverse accepted_hits.bam $gtfref > $runid.htseq.txt