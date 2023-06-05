#sh rename_files.sh to run
#input: newnames.csv - csv file of format oldname,newname
#File for methylseq of mature flowers - in Seedling_Phenology folder
#############################################
while IFS=, read -r oldname newname
do
    mv $oldname $newname
done < newnames.csv