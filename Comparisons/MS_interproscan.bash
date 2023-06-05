#############################################################
#Trying to run Interproscan on DMRS:
#Note: In progress
#Done so far:
    #Installation
    #Testing on given data (works without -dp parameter, potentially try with dependencies??)
    #
#Using: blast_ref_intersect_SE_W_C_SE_L_H.out , for testing -> later update to all necessary files
#############################################################
cd scratch
mkdir interproscan
cd interproscan
#############################################################
#Install interproscan: https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html#index-hmm-models

#cedar1
tmux new-session -s Interproscan
tmux attach-session -t Interproscan
#From website:
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz.md5
#here for now
# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.62-94.0-64-bit.tar.gz.md5
# Must return *interproscan-5.62-94.0-64-bit.tar.gz: OK*

#Extract from tar
tar -pxvzf interproscan-5.62-94.0-*-bit.tar.gz

module load python
cd interproscan-5.62-94.0
python3 setup.py -f interproscan.properties
##############################################################
#Testing interproscan:
#Required java version:

#Try:
module load java/11.0.16_8
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
#Fails: not enough space.
#Try wil allocation:
salloc -c1 --time 1:00:00 --mem 120000m --account def-rieseber
module load java/11.0.16_8

#Try not within directory
cd..
interproscan-5.62-94.0/interproscan.sh -i interproscan-5.62-94.0/test_all_appl.fasta -f tsv -dp
#Permission errors 

#Try within folder
cd interproscan-5.62-94.0
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
#no memory error, many exceptions thrown, but time ran out 

#Try within folder again:
salloc -c1 --time 2:00:00 --mem 120000m --account def-rieseber
module load java/11.0.16_8
cd interproscan-5.62-94.0
./interproscan.sh -i test_all_appl.fasta -f tsv -dp

#Trying without -dp param
./interproscan.sh -i test_all_appl.fasta -f tsv 
#Works!!! - leave for now 

cd..
interproscan-5.62-94.0/interproscan.sh -i interproscan-5.62-94.0/test_all_appl.fasta -f tsv 
#Works without -dp as well

#Clean up directory:
rm *.tsv


##############################################################
#Trying with my own files:

#Set up desired fasta file:t
cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_output/blast_ref_intersect_Wild_W_C_Mat_Sen.out .
cp ~/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.transcript.fa  


git clone https://github.com/lh3/seqtk.git;
cd seqtk; make

module load StdEnv/2020
module load gcc/9.3.0
module laod r-bundle-bioconductor/3.16 

#.fasta and .fa - both file extensions for same format, but irterproscan took .fasta 
seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa blast_ref_intersect_Wild_W_C_Mat_Sen.out > intersect_map.fasta
#Empty file

######
#Trying without sebseq:
awk 'BEGIN{while((getline<"blast_ref_intersect_Wild_W_C_Mat_Sen.out")>0)l[">"$1]=1}/^>/{f=l[$1]}f' Dryas_octopetala_H1.protein.fa > intersect_map2.fasta
#empty file

#Trying, extract second column of blast output:
cut -f2 blast_ref_intersect_Wild_W_C_Mat_Sen.out > blast_out.txt

seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa blast_out.txt > intersect_map.fasta
#Works!!!

#Try run:
interproscan-5.62-94.0/interproscan.sh -i intersect_map.fasta -f tsv 
#Fails :not enough space

#Try with allocation:
salloc -c1 --time 2:00:00 --mem 120000m --account def-rieseber
module load java/11.0.16_8 #Need ve
interproscan-5.62-94.0/interproscan.sh -i intersect_map.fasta -f tsv 
# fails: try test again:
interproscan-5.62-94.0/interproscan.sh -i interproscan-5.62-94.0/test_all_appl.fasta -f tsv 
#Works: have all dependencies correct:
#Test with multiple protiens:
interproscan-5.62-94.0/interproscan.sh -i interproscan-5.62-94.0/test_proteins.fasta -f tsv
#Works: with protien files?

#Why exceptions thrown?
    #Gene sequence has *'s in it.
    #Possible solution: Filter those with * out, ask
salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber


##2)Filtering out any asterisk:
awk '/^>/ { if (seq != "" && seq !~ /\*/) { print header "\n" seq } header = $0; seq = "" } /^[^>]/ { seq = seq $0 } \
END { if (seq != "" && seq !~ /\*/) { print header "\n" seq } }' intersect_map.fasta > intersect_map2.fasta

##1) Remove duplicates: https://www.biostars.org/p/165001/
#awk '/^>/{f=!d[$1];d[$1]=1}f' intersect_map.fasta > intersect_map2.fasta
interproscan-5.62-94.0/interproscan.sh -i intersect_map2.fasta -f tsv 

#'ageListenerInvoker.run(DefaultMessageListenerContainer.java:1073)
#        at java.base/java.lang.Thread.run(Thread.java:829)
#2023-06-02 14:29:08,432 [amqEmbeddedWorkerJmsContainer-4] [uk.ac.ebi.interpro.scan.jms.worker.LocalJobQueueListener:222] ERROR - StepExecution with errors - stepName: stepPrositePatternRunBinary
#2023-06-02 14:29:08,595 [main] [uk.ac.ebi.interpro.scan.jms.master.StandaloneBlackBoxMaster:190] WARN - StepInstance 16 is being re-run following a failure.
#
#Failed with:
#InterProScan analysis failed. Exception thrown by StandaloneBlackBoxMaster. Check the log file for details

#Try with more cpus
salloc -c4 --time 7:00:00 --mem 120000m --account def-rieseber
module load java/11.0.16_8
#Try again:
interproscan-5.62-94.0/interproscan.sh -i intersect_map2.fasta -f tsv 