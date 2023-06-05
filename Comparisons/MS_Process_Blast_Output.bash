##########################################
#Processing Blast output files - to have clean gene names
#Goal - remove any duplicates, leave only fragments with gene names
#Note: This file is mainly scratch work for now
##########################################

cd scratch
mkdir process_blast_output
cd process_blast_output

tmux new-session -s Process_Blast
tmux attach-session -t Process_Blast