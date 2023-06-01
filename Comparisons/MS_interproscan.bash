#############################################################
#Trying to run Interproscan on DMRS:
#Note: In progress
#https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html#obtaining-the-core-interproscan-software
#Using: blast_ref_intersect_SE_W_C_SE_L_H.out , for testing -> later update to all necessary files
#############################################################

cd scratch
mkdir interproscan
cd interproscan

#############################################################
#Install interproscan: 

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
#Set up other files:


#cp blast_ref_intersect_SE_W_C_SE_L_H.out .