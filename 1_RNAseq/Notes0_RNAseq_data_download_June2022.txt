########################################
# Dryas differential gene expression data
# Data download
# First download failed but second worked
######################################
# June 3 2022
# data download from Genome Quebec

    # In the project, select the appropriate 'Read sets' tab.
    # Click the main check box to select all lines of the table or specific read sets by checking the checkbox on the line.
    # Then click the 'Download Read Files' button.
    # On the popup window, select the second option 'Download files from selected reads'.
    # Some new options will show up. Select the 'MD5' one.
    # Make sure the desired file formats are selected by filling the checkboxes. Click then 'Download'.
    # You will get a readSets.md5 file containing the MD5 checksum values for the selected readsets.

# https://ces.genomequebec.com/nanuqMPS/project/ProjectPage/projectId/21617

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21617&tech=NovaSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

# downloading

FINISHED --2022-06-04 23:49:43--
Total wall clock time: 1d 2h 11m 39s
Downloaded: 100 files, 178G in 1d 2h 6m 4s (1.94 MB/s)


#################################
# look at raw data

cp -v NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11* ..
gunzip NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11*
gzip: NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R2.fastq.gz: invalid compressed data--crc error
gzip: NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R2.fastq.gz: invalid compressed data--length error

# maybe error in this file??
# what about other R2 files??
# also in some other files

#################################
# try data download again

tmux new-session -s RNA_download
tmux attach-session -t RNA_download

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21617&tech=NovaSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

# downloading more from some...

ls -1 | wc -l
# should have 4*70 files = 280
287 files

#-----------------------------------
# check gunzip on all and any that fail delete and re download

gunzip -kv *.fastq.gz

rm NS.1889.002.NEBNext_dual_i7_104---NEBNext_dual_i5_104.T_W12_R1.fastq.gz
rm NS.1889.002.NEBNext_dual_i7_107---NEBNext_dual_i5_107.T_W2_R2.fastq.gz
rm NS.1889.002.NEBNext_dual_i7_111---NEBNext_dual_i5_111.A_C2_R1.fastq.gz
rm NS.1889.002.NEBNext_dual_i7_112---NEBNext_dual_i5_112.A_C6_R1.fastq.gz

# since many files are corrupt

#-------------------------------
# just restart the download

# remove all files

tmux attach-session -t RNA_download
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

rm *

read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21617&tech=NovaSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

FINISHED --2022-09-09 03:53:56--
Total wall clock time: 4h 51m 47s
Downloaded: 284 files, 669G in 4h 49m 38s (39.4 MB/s)

#-----------------------
ls -1 | wc -l
# should be 284

#-----------------------------
# check md5 sums 

# on local machine go to directory with downloaded file
cd /home/Owner/MyDocuments/Hackery\ Backup/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/RNAseq/Data
scp -v readSets_*.md5 celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

# back on Cedar or your server
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

md5sum -c readSets_fastq.md5

# should look like this when checking
NS.1889.003.NEBNext_dual_i7_138---NEBNext_dual_i5_138.L_W10_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_138---NEBNext_dual_i5_138.L_W10_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_157---NEBNext_dual_i5_157.N_C3_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_157---NEBNext_dual_i5_157.N_C3_B_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_161---NEBNext_dual_i5_161.N_W4_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_161---NEBNext_dual_i5_161.N_W4_B_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_D8---NEBNext_dual_i5_D8.N_W11_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_D8---NEBNext_dual_i5_D8.N_W11_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_160---NEBNext_dual_i5_160.N_W3_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_160---NEBNext_dual_i5_160.N_W3_B_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_149---NEBNext_dual_i5_149.A_C24_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_149---NEBNext_dual_i5_149.A_C24_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_147---NEBNext_dual_i5_147.A_C22_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_147---NEBNext_dual_i5_147.A_C22_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_153---NEBNext_dual_i5_153.A_W24_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_153---NEBNext_dual_i5_153.A_W24_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_151---NEBNext_dual_i5_151.A_W22_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_151---NEBNext_dual_i5_151.A_W22_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_146---NEBNext_dual_i5_146.A_C21_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_146---NEBNext_dual_i5_146.A_C21_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_144---NEBNext_dual_i5_144.L_W5_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_144---NEBNext_dual_i5_144.L_W5_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_148---NEBNext_dual_i5_148.A_C23_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_148---NEBNext_dual_i5_148.A_C23_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_150---NEBNext_dual_i5_150.A_W21_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_150---NEBNext_dual_i5_150.A_W21_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_152---NEBNext_dual_i5_152.A_W23_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_152---NEBNext_dual_i5_152.A_W23_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_162---NEBNext_dual_i5_162.N_W8_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_162---NEBNext_dual_i5_162.N_W8_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_145---NEBNext_dual_i5_145.L_W6_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_145---NEBNext_dual_i5_145.L_W6_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_143---NEBNext_dual_i5_143.L_W2_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_143---NEBNext_dual_i5_143.L_W2_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_134---NEBNext_dual_i5_134.L_C2_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_134---NEBNext_dual_i5_134.L_C2_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_142---NEBNext_dual_i5_142.L_W17_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_142---NEBNext_dual_i5_142.L_W17_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_135---NEBNext_dual_i5_135.L_C3_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_135---NEBNext_dual_i5_135.L_C3_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_136---NEBNext_dual_i5_136.L_C5_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_136---NEBNext_dual_i5_136.L_C5_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_140---NEBNext_dual_i5_140.L_W14_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_140---NEBNext_dual_i5_140.L_W14_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_132---NEBNext_dual_i5_132.L_C15_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_132---NEBNext_dual_i5_132.L_C15_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_B8---NEBNext_dual_i5_B8.T_C5_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_B8---NEBNext_dual_i5_B8.T_C5_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_137---NEBNext_dual_i5_137.L_C8_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_137---NEBNext_dual_i5_137.L_C8_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_155---NEBNext_dual_i5_155.N_C6_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_155---NEBNext_dual_i5_155.N_C6_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_C8---NEBNext_dual_i5_C8.T_C7_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_C8---NEBNext_dual_i5_C8.T_C7_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_141---NEBNext_dual_i5_141.L_W16_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_141---NEBNext_dual_i5_141.L_W16_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_139---NEBNext_dual_i5_139.L_W13_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_139---NEBNext_dual_i5_139.L_W13_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_133---NEBNext_dual_i5_133.L_C16_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_133---NEBNext_dual_i5_133.L_C16_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_158---NEBNext_dual_i5_158.N_C9_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_158---NEBNext_dual_i5_158.N_C9_B_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_154---NEBNext_dual_i5_154.A_W6_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_154---NEBNext_dual_i5_154.A_W6_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_159---NEBNext_dual_i5_159.N_W10_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_159---NEBNext_dual_i5_159.N_W10_B_R2.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_156---NEBNext_dual_i5_156.N_C2_B_R1.fastq.gz: OK
NS.1889.003.NEBNext_dual_i7_156---NEBNext_dual_i5_156.N_C2_B_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_108---NEBNext_dual_i5_108.T_W3_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_108---NEBNext_dual_i5_108.T_W3_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_110---NEBNext_dual_i5_110.T_W9_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_110---NEBNext_dual_i5_110.T_W9_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_129---NEBNext_dual_i5_129.SE2W_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_129---NEBNext_dual_i5_129.SE2W_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_101---NEBNext_dual_i5_101.T_C16_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_101---NEBNext_dual_i5_101.T_C16_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_99---NEBNext_dual_i5_99.T_C14_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_99---NEBNext_dual_i5_99.T_C14_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_126---NEBNext_dual_i5_126.SE2C_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_126---NEBNext_dual_i5_126.SE2C_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_115---NEBNext_dual_i5_115.A_C7_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_115---NEBNext_dual_i5_115.A_C7_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_128---NEBNext_dual_i5_128.SE1W_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_128---NEBNext_dual_i5_128.SE1W_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_98---NEBNext_dual_i5_98.T_C12_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_98---NEBNext_dual_i5_98.T_C12_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_103---NEBNext_dual_i5_103.T_W1_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_103---NEBNext_dual_i5_103.T_W1_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_105---NEBNext_dual_i5_105.T_W13_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_105---NEBNext_dual_i5_105.T_W13_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_124---NEBNext_dual_i5_124.N_W7_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_124---NEBNext_dual_i5_124.N_W7_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_127---NEBNext_dual_i5_127.SE3C_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_127---NEBNext_dual_i5_127.SE3C_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_107---NEBNext_dual_i5_107.T_W2_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_107---NEBNext_dual_i5_107.T_W2_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_100---NEBNext_dual_i5_100.T_C15_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_100---NEBNext_dual_i5_100.T_C15_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_122---NEBNext_dual_i5_122.N_W2_R1.fastq.gz: OK                                                                            [0/1970]
NS.1889.002.NEBNext_dual_i7_122---NEBNext_dual_i5_122.N_W2_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_97---NEBNext_dual_i5_97.T_C1_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_97---NEBNext_dual_i5_97.T_C1_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_131---NEBNext_dual_i5_131.L_C14_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_131---NEBNext_dual_i5_131.L_C14_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_125---NEBNext_dual_i5_125.SE1C_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_125---NEBNext_dual_i5_125.SE1C_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_130---NEBNext_dual_i5_130.SE3W_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_130---NEBNext_dual_i5_130.SE3W_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_102---NEBNext_dual_i5_102.T_C20_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_102---NEBNext_dual_i5_102.T_C20_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_123---NEBNext_dual_i5_123.N_W6_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_123---NEBNext_dual_i5_123.N_W6_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_116---NEBNext_dual_i5_116.A_W11_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_116---NEBNext_dual_i5_116.A_W11_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_104---NEBNext_dual_i5_104.T_W12_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_104---NEBNext_dual_i5_104.T_W12_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_120---NEBNext_dual_i5_120.N_C5_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_120---NEBNext_dual_i5_120.N_C5_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_118---NEBNext_dual_i5_118.N_C1_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_118---NEBNext_dual_i5_118.N_C1_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_114---NEBNext_dual_i5_114.A_C11_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_114---NEBNext_dual_i5_114.A_C11_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_112---NEBNext_dual_i5_112.A_C6_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_112---NEBNext_dual_i5_112.A_C6_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_119---NEBNext_dual_i5_119.N_C10_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_119---NEBNext_dual_i5_119.N_C10_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_121---NEBNext_dual_i5_121.N_C7_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_121---NEBNext_dual_i5_121.N_C7_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_117---NEBNext_dual_i5_117.A_W7_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_117---NEBNext_dual_i5_117.A_W7_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_113---NEBNext_dual_i5_113.A_W1_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_113---NEBNext_dual_i5_113.A_W1_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_106---NEBNext_dual_i5_106.T_W15_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_106---NEBNext_dual_i5_106.T_W15_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_109---NEBNext_dual_i5_109.T_W5_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_109---NEBNext_dual_i5_109.T_W5_R2.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_111---NEBNext_dual_i5_111.A_C2_R1.fastq.gz: OK
NS.1889.002.NEBNext_dual_i7_111---NEBNext_dual_i5_111.A_C2_R2.fastq.gz: OK

#-------------------------------------
# tar input folder and save in Dryas folder
tar -zcvf RNAseq_data_Sept2022.tar.gz /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/

#####################################
# try downloading to here as well - fresh download

tmux new-session -s RNA
tmux attach-session -t RNA

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/

read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21617&tech=NovaSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

ls -1 | wc -l
# should have 4*70 files = 280

FINISHED --2022-09-08 21:43:15--
Total wall clock time: 6h 31m 9s
Downloaded: 284 files, 677G in 6h 28m 46s (29.7 MB/s)

#----------------------------
# try unzipping all files
gunzip -kv *.fastq.gz

#------------------------------
# check files

cp -v /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data/readSets_*.md5 /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/

md5sum -c readSets_fastq.md5

# all okay

#---------------------------------
# move non .gz files to different folder
mkdir fastq_files
mv *.fastq fastq_files/
mv NS*_I2.fastq index_files/

########################
