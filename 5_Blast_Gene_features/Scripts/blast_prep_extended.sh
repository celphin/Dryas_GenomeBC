#########################################
#blast_prep_extended.sh
    #Maps bedgraph onto fasta file
    #Ensures no zero 
    #Parameter: blastname, string
    #creates blastname.fasta file
#########################################
# https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
# bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>
#########################################
blastname=$1
directory="blast_bedgraphs"
ref_genome=Dryas_octopetala_H1.supercontigs.fa

cat "${directory}/${blastname}.bedGraph" | awk '{print $1 "\t" $2 "\t" $3}' > "${blastname}.bed"

# add 500bp on each side of the DMR
awk '{$3+=1000}1' "${blastname}.bed" > "${blastname}_extended.bed"
awk '{$2-=1000}1' "${blastname}_extended.bed" > "${blastname}_extended2.bed"

#changing all negative values to zero:
awk '{if($2<0) sub($2,0)}1' "${blastname}_extended2.bed" > "${blastname}_extended3.bed"

#formatting
sed -e 's/ /\t/g' "${blastname}_extended3.bed" > "${blastname}_extended3_tab.bed"

#get fasta
bedtools getfasta -fi $ref_genome \
-bed "${blastname}_extended3_tab.bed" \
> "${blastname}.fasta"

#clear intermediate files
rm *.bed
############################################################