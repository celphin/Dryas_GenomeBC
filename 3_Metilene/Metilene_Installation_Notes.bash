#######################################################
#Notes on required installations for Metilene folder
    #Installing ggplot to compute canada
    #Installing metilene to desired directory
######################################################

#Installing ggplot: https://docs.alliancecan.ca/wiki/
module load nixpkgs/16.09 
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python

#Install R ggplot2 : https://docs.alliancecan.ca/wiki/R
mkdir /home/username/R/x86_64-pc-linux-gnu-library/3.6/
export R_LIBS_USER=/home/username/R/x86_64-pc-linux-gnu-library/3.6/
R -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")'

#####################################################
#Installing metilene
#In desired directory (specified in script later)
wget http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-8.tar.gz
tar -xvzf metilene_v02-8.tar.gz
cd metilene_v02-8
make
#####################################################