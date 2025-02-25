##############################
# Notes to make a circos plot
# Dec 2024
##############################
#Step 1: Install Circos and Prepare Tools
#First, install Circos. It’s a software for creating circular visualizations, and it’s available from Circos website.

cd /home/celphin/scratch/Dryas/circos/data

#############################
#Step 2: Organize Data

# 1. bedGraph Files
# Each bedGraph file will correspond to a separate track in the Circos plot.
# You might need to convert the bedGraph file into a more suitable form, like:
# chromosome start_position end_position value

# CpG methylation
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/RagTag_C_SwedC.bedGraph .

# CHH methylation
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/CHH/RagTag/RagTag_C_SwedC_CHH.bedGraph .

# difference in methylation in Sweden
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/RagTag_W_SVAL.bedGraph .
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/CHH/RagTag/RagTag_W_SwedC_CHH.bedGraph .

# DMRs
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/BedGraphs_From_Metilene/RagTag_bedGraph/RagTag_Sweden_W_C.bedGraph .
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/CHH/RagTag/RagTag_Sweden_CHH_W_C_CHH.bedGraph .

# FST - note diff ref
cp /home/celphin/scratch/Dryas/BS-Snper/pop_gen/FST/Sweden_W_C.weir.fst .

# Centromeres
DoctH0-1_RagTag	26480000	27330000
DoctH0-2_RagTag	8135000	8810000
DoctH0-3_RagTag	16830000	16920000
DoctH0-4_RagTag	4690000	5340000
DoctH0-5_RagTag	1695000	2635000
DoctH0-6_RagTag	11560000	11600000
DoctH0-7_RagTag	16840000	17670000
DoctH0-8_RagTag	16785000	17405000
DoctH0-9_RagTag	2695000	2995000

# Reference 
cp /home/celphin/scratch/Dryas/CrossMap/OldRagtagDryOcto.fasta .
 

#-------------------------
#2. Gene Annotations
#If your gene annotations are in GFF or GTF format, you can convert them into a BED format, which Circos can use. You can use tools like gff2bed for this:

cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/RagTag_Dryas_genes.gff3 .
module load StdEnv/2023 bedops/2.4.41 
gff2bed <  RagTag_Dryas_genes.gff3  >  RagTag_Dryas_genes.bed

#Make sure your annotations are in the correct format:
#chromosome start_position end_position gene_id

# extract columns
grep "CDS" RagTag_Dryas_genes.bed > RagTag_Dryas_genes1.bed
awk -F'\t' '{print $1, $2, $3, $10}' RagTag_Dryas_genes1.bed | sed 's/ID=\([^;]*\);.*/\1/'  > RagTag_Dryas_genes2.bed

#--------------------------
#3. TE Annotations
#TE annotations should also be in a BED format (or converted if necessary). 
#Similar to gene annotations, these will contain the TE locations.

cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/RagTag_Dryas_TEs.gff3 .
gff2bed < RagTag_Dryas_TEs.gff3 >  RagTag_Dryas_TEs.bed

#Example:
#chromosome start_position end_position TE_type

awk -F'\t' '{print $1, $2, $3, $8}'  RagTag_Dryas_TEs.bed >  RagTag_Dryas_TEs1.bed


#----------------------
# Karyotypes
seqkit stats OldRagtagDryOcto.fasta

grep -A 1 ">" OldRagtagDryOcto.fasta | grep -v ">" | awk '{print length($0)}'

nano karyotype.Dryas.txt
chr	- DoctH0-1_RagTag	1	0	28726255	
chr	- DoctH0-2_RagTag	2	0	29810504	chr2
chr	- DoctH0-3_RagTag	3	0	29331009	chr3
chr	- DoctH0-4_RagTag	4	0	27845443	chr4
chr	- DoctH0-5_RagTag	5	0	21333714	chr5
chr	- DoctH0-6_RagTag	6	0	22319059	chr6
chr	- DoctH0-7_RagTag	7	0	26945421	chr7
chr	- DoctH0-8_RagTag	8	0	20495270	chr8
chr	- DoctH0-9_RagTag	9	0	19685959	chr9


#################################
#Step 3: Create a Configuration File for Circos

#Circos uses a configuration file that specifies how the data should be plotted. 
# Below is a basic structure of a Circos configuration file (circos.conf):
nano circos.conf

# # Circos configuration file

# Chromosome name, size and color definition
karyotype = data/karyotype.Dryas.txt

<ideogram>

<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 20p
fill      = yes

</ideogram>

<image>
# Included from Circos distribution.
#<<include etc/image.conf>>
</image>

# Plotting the methylation bedGraph file as a track
<plots>
  <plot>
    type = histogram
    file = data/RagTag_C_SwedC.bedGraph
    r0   = 0.75r
    r1   = 1r
    color = red
	width=2p
  </plot>
</plots>

# Gene annotations track
<plots>
  <plot>
    type = line
    file = data/RagTag_Dryas_genes2.bed
    r0   = 0.6r
    r1   = 0.7r
    color = blue
  </plot>
</plots>

# Transposable element annotations track
<plots>
  <plot>
    type = line
    file = data/RagTag_Dryas_TEs1.bed
    r0   = 0.55r
    r1   = 0.6r
    color = green
  </plot>
</plots>



# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
#<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
#<<include etc/housekeeping.conf>>

###############################
#Step 4: Build the Circos Plot

cd /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/circos/0.69-9/
mkdir /home/celphin/scratch/Dryas/circos/etc/
cp ./etc/* /home/celphin/scratch/Dryas/circos/etc/
cd  etc/
mkdir /home/celphin/scratch/Dryas/circos/etc/tracks/
cp ./tracks/* /home/celphin/scratch/Dryas/circos/etc/tracks/
cd /home/celphin/scratch/Dryas/circos/


#Once your configuration file (circos.conf) and data files (e.g., bedgraph_data.txt, gene_annotations.bed, 
#te_annotations.bed) are ready, you can run Circos to generate the plot:

module load StdEnv/2023 circos/0.69-9
circos -conf circos.conf

  # Tried to obtain (x,y) coordinate from angle [-89.9945114346755] and radius
  # [_undef_], but one of these was not defined.

  # If you are having trouble debugging this error, first read the best practices
  # tutorial for helpful tips that address many common problems

      # http://www.circos.ca/documentation/tutorials/reference/best_practices

