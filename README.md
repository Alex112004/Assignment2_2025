Forma de utilizar:
sudo apt update
micromamba activate
git clone https://github.com/hawaiidatascience/metaflowmics.git

cd metaflowmics/metaflowmics/Pipeline-16S/
nextflow run main.nf -c nextflow.config -profile docker -resume
