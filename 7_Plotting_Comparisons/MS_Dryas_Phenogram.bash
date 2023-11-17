#########################################################
#November 2023
#Steps:
    #Get 
    #Convert to Gene: start-end
    #Repeat with 
    #Get genome file for phenogram
#########################################################


wget https://ritchielab.org/files/RL_software/ruby_install.sh
wget https://ritchielab.org/files/RL_software/pheno_gram.rb
sh ruby_install.sh
ruby pheno_gram.rb
#######################################################
#Acquire Dryas_Genome.txt:

module load bioawk
bioawk -c fastx '{ print $name, length($seq) }' < Dryas_octopetala_H1.supercontigs.fa

#Hand Edit to be tsv file, with header:
#######################################################
grep -E "SE_W_C|Wild_W_C|Parent_W_C|POS|RNA" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_W_C_SE_Wild_P_RNA_Phenogram.txt   
grep -E "SE_L_H|Wild_Lat_L_H|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_L_H_SE_Wild_Phenogram.txt   
grep -E "Sval|Alas|Swed|Nunavut|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_Loc_W_C_Phenogram.txt   
grep -E "Mat_Sen|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_Phenology_Phenogram.txt  
#######################################################
ruby pheno_gram.rb -i Dryas_W_C_SE_Wild_P_RNA_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Warming Phenogram" -o Dryas_W_C_SE_Wild_P_RNA_Phenogram -f jpg
ruby pheno_gram.rb -i Dryas_Loc_W_C_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Location Phenogram" -o Dryas_Loc_W_C_Phenogram -f jpg
ruby pheno_gram.rb -i Dryas_Phenology_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Phenology Phenogram" -o Dryas_Phenology_Phenogram -f jpg
#Too many, figure out way to subset
#ruby pheno_gram.rb -i Dryas_L_H_SE_Wild_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Latitude Phenogram" -o Dryas_L_H_SE_Wild_Phenogram -f jpg
#######################################################\
scp -v msandler@cedar.computecanada.ca:/home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Dryas_Phenogram/*.jpg .


