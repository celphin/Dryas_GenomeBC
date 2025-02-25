###########################################################
# Sum Bedgegraph files from methylseq for viewing in IGV
# Nov 2024
#########################################

tmux new-session -s bedgraph
tmux attach-session -t bedgraph

cd ~/scratch/Dryas/bedgraph_old_ref

cp -r /home/celphin/projects/rpp-rieseber/celphin/Dryas/Methylation_calling/May2021_Parents/bismark_methylation_calls/bedGraph* ~/scratch/Dryas/bedgraph_old_ref

gunzip *deduplicated.bedGraph.gz

module load StdEnv/2023 bedtools/2.31.0

# test
bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph | \
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' > W_SVAL.bedGraph

bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph > W_SVAL_total.bedGraph

head W_SVAL_total.bedGraph
# Do1_00107       8       9       0       0       0       0       0       0
# Do1_00107       9       10      0       0       0       0       0       0
# Do1_00107       137     138     12.5    0       16.6666666666667        0       0       0
# Do1_00107       138     139     28.5714285714286        0       7.69230769230769        12.5 00
# Do1_00107       247     248     0       0       46.1538461538462        0       0       14.2857142857143
# Do1_00107       248     249     0       0       8.69565217391304        0       0       0
# Do1_00107       572     573     30      50      40      0       25      38.0952380952381
# Do1_00107       573     574     41.6666666666667        50      27.2727272727273        57.1428571428571      62.5    54.5454545454545
# Do1_00107       580     581     77.7777777777778        0       75      33.3333333333333     027.2727272727273
# Do1_00107       581     582     82.6086956521739        25      36.3636363636364        16.6666666666667      50      81.8181818181818

head W_SVAL.bedGraph
# Do1_00107       8       9       0
# Do1_00107       9       10      0
# Do1_00107       137     138     4.86111
# Do1_00107       138     139     8.12729
# Do1_00107       247     248     10.0733
# Do1_00107       248     249     1.44928
# Do1_00107       572     573     30.5159
# Do1_00107       573     574     48.8546
# Do1_00107       580     581     35.564
# Do1_00107       581     582     48.7429

# Working!

#############################
# Now for all combinations

module load StdEnv/2023 bedtools/2.31.0

bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph > W_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_SVAL_total.bedGraph > W_SVAL.bedGraph

bedtools unionbedg -i C*SVAL*.deduplicated.bedGraph > C_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_SVAL_total.bedGraph > C_SVAL.bedGraph

#---------------
bedtools unionbedg -i W*LAT*.deduplicated.bedGraph > W_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_Swed_total.bedGraph > W_SwedC.bedGraph

bedtools unionbedg -i C*LAT*.deduplicated.bedGraph > C_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_Swed_total.bedGraph > C_SwedC.bedGraph

#---------------
bedtools unionbedg -i W*CASS*.deduplicated.bedGraph > W_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_CASS_total.bedGraph > W_CASS.bedGraph

bedtools unionbedg -i C*CASS*.deduplicated.bedGraph  > C_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_CASS_total.bedGraph > C_CASS.bedGraph

#---------------
bedtools unionbedg -i W*WILL*.deduplicated.bedGraph > W_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_WILL_total.bedGraph > W_WILL.bedGraph

bedtools unionbedg -i C*WILL*.deduplicated.bedGraph > C_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_WILL_total.bedGraph > C_WILL.bedGraph

#---------------
bedtools unionbedg -i W*FERT*.deduplicated.bedGraph > W_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_FERT_total.bedGraph > W_FERT.bedGraph

bedtools unionbedg -i C*FERT*.deduplicated.bedGraph > C_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_FERT_total.bedGraph > C_FERT.bedGraph

#---------------
bedtools unionbedg -i W*DRY*.deduplicated.bedGraph > W_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_DRY_total.bedGraph > W_DRY.bedGraph

bedtools unionbedg -i C*DRY*.deduplicated.bedGraph > C_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_DRY_total.bedGraph > C_DRY.bedGraph

#---------------
bedtools unionbedg -i W*MEAD*.deduplicated.bedGraph > W_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_MEAD_total.bedGraph > W_MEAD.bedGraph

bedtools unionbedg -i C*MEAD*.deduplicated.bedGraph > C_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_MEAD_total.bedGraph > C_MEAD.bedGraph

#---------------
bedtools unionbedg -i W*ALAS*.deduplicated.bedGraph  > W_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_ALAS_total.bedGraph > W_ALAS.bedGraph

bedtools unionbedg -i C*ALAS*.deduplicated.bedGraph > C_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_ALAS_total.bedGraph > C_ALAS.bedGraph

####################################
# Total bedGraph intersections

tmux new-session -s bedgraph
tmux attach-session -t bedgraph

mkdir /home/celphin/scratch/Dryas/DMR_bedGraph_files
cd /home/celphin/scratch/Dryas/DMR_bedGraph_files

# CHH methylation and DMRs
cp /home/celphin/scratch/Dryas/CHG_CHH/*.bedGraph .

# DMRs
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/BedGraphs_From_Metilene/*.bedGraph .
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/Bedgraphs_Intersected_Subtracted/*.bedGraph .
cp /home/celphin/scratch/Dryas/CrossMap/file_for_converting/Bedgraphs_Intersected_Subtracted/*.bedgraph .


#-------------------
# test
module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/DMR_bedGraph_files

#Bedtools intersect all but Svalbard
bedtools intersect -u -a X.bedGraph -b Y.bedGraph  > X_Y_intersect.bedgraph


#-------------------
# Run on all

nano intersect_bedGraph.sh

#!/bin/bash

# Set the directory containing the bedGraph files
bedGraph_dir="/home/celphin/scratch/Dryas/DMR_bedGraph_files"  
output_dir="/home/celphin/scratch/Dryas/DMR_bedGraph_files/intersections" 
matrix_file="intersection_matrix.txt"  

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Create an array of all bedGraph files in the directory
bedGraph_files=($bedGraph_dir/*.bedGraph)

# Initialize matrix as an empty string
matrix=""

# First, create the header row with file names (to be the column names in the matrix)
header_row=""

for ((i=0; i<${#bedGraph_files[@]}; i++)); do
  header_row+=$(basename "${bedGraph_files[$i]}")$'\t'
done

# Remove the last tab from the header row
header_row=${header_row%$'\t'}

# Add the header to the matrix
matrix+="$header_row"$'\n'

# Loop over all pairs of bedGraph files
for ((i=0; i<${#bedGraph_files[@]}; i++)); do
  # Create a row for the i-th file (start with the file name)
  row="$(basename "${bedGraph_files[$i]}")"

  for ((j=0; j<${#bedGraph_files[@]}; j++)); do
    # Get the file names
    file1="${bedGraph_files[$i]}"
    file2="${bedGraph_files[$j]}"

    # If it's the diagonal (when i == j), we want to count the total number of lines in the file
    if [ "$i" -eq "$j" ]; then
      intersection_count=$(wc -l < "$file1")
    else
      # Otherwise, run bedtools intersect to get the intersection and the count
      intersection_output="$output_dir/$(basename "$file1")_$(basename "$file2")_intersection.bed"
      intersection_count=$(bedtools intersect -a "$file1" -b "$file2" -u > "$intersection_output"; wc -l < "$intersection_output")
    fi

    # Append the intersection count to the current row, separated by a tab
    row+="$intersection_count"$'\t'
  done

  # Remove the last tab from the row and add it to the matrix
  row=${row%$'\t'}
  matrix+="$row"$'\n'
done

# Save the matrix to the file
echo -e "$matrix" > "$matrix_file"

# Print a message indicating the file was saved
echo "Matrix saved to $matrix_file"


#-------------------------------
# extended bedGraph files

cd /home/celphin/scratch/Dryas/DMR_bedGraph_files
mkdir extensions

# Set the directory containing your .bedGraph files
input_directory="/home/celphin/scratch/Dryas/DMR_bedGraph_files"
output_directory="/home/celphin/scratch/Dryas/DMR_bedGraph_files/extensions"

# Loop through all .bedGraph files in the input directory
for file in "$input_directory"/*.bedGraph; do
    # Get the base name of the file (without path)
    filename=$(basename "$file")
    
    # Apply awk to modify the second and third columns, ensuring tab-delimited output
    awk -F'\t' 'BEGIN {OFS="\t"} { $2=$2-500; $3=$3+500; print }' "$file" > "$output_directory/extended_$filename"
    
    echo "Processed $file and saved to $output_directory/extended_$filename"
done


# Loop through all .bedGraph files in the input directory
# remove negatives

for file in "$input_directory"/*.bedGraph; do
    # Get the base name of the file (without path)
    filename=$(basename "$file")
    
    # Apply awk to filter out lines with negative second column
    awk -F'\t' 'BEGIN {OFS="\t"} { if ($2 >= 0) print }' "$file" > "$output_directory/noneg_$filename"
    
    echo "Processed $file and saved to $output_directory/$filename"
done
#----------------------
# Set the directory containing the bedGraph files
bedGraph_dir="/home/celphin/scratch/Dryas/DMR_bedGraph_files/extensions"  
output_dir="/home/celphin/scratch/Dryas/DMR_bedGraph_files/extensions/intersections" 
matrix_file="extended_intersection_matrix.txt"  

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Create an array of all bedGraph files in the directory
bedGraph_files=($bedGraph_dir/noneg_*.bedGraph)

# Initialize matrix as an empty string
matrix=""

# First, create the header row with file names (to be the column names in the matrix)
header_row=""

for ((i=0; i<${#bedGraph_files[@]}; i++)); do
  header_row+=$(basename "${bedGraph_files[$i]}")$'\t'
done

# Remove the last tab from the header row
header_row=${header_row%$'\t'}

# Add the header to the matrix
matrix+="$header_row"$'\n'

# Loop over all pairs of bedGraph files
for ((i=0; i<${#bedGraph_files[@]}; i++)); do
  # Create a row for the i-th file (start with the file name)
  row="$(basename "${bedGraph_files[$i]}")"

  for ((j=0; j<${#bedGraph_files[@]}; j++)); do
    # Get the file names
    file1="${bedGraph_files[$i]}"
    file2="${bedGraph_files[$j]}"

    # If it's the diagonal (when i == j), we want to count the total number of lines in the file
    if [ "$i" -eq "$j" ]; then
      intersection_count=$(wc -l < "$file1")
    else
      # Otherwise, run bedtools intersect to get the intersection and the count
      intersection_output="$output_dir/$(basename "$file1")_$(basename "$file2")_intersection.bed"
      intersection_count=$(bedtools intersect -a "$file1" -b "$file2" -u > "$intersection_output"; wc -l < "$intersection_output")
    fi

    # Append the intersection count to the current row, separated by a tab
    row+="$intersection_count"$'\t'
  done

  # Remove the last tab from the row and add it to the matrix
  row=${row%$'\t'}
  matrix+="$row"$'\n'
done

# Save the matrix to the file
echo -e "$matrix" > "$matrix_file"

# Print a message indicating the file was saved
echo "Matrix saved to $matrix_file"
