**Forma de utilizar:**


Atualizar a vm:

//t sudo apt update

Ativar o micromamba/conda:

/t micromamba activate / conda activate

Instalar o :

/t git clone https://github.com/hawaiidatascience/metaflowmics.git

Abir a pasta da pipeline do 16S:

/t cd metaflowmics/metaflowmics/Pipeline-16S/

Correr nextflow com as novvas configurações e utilizar o comando -profile docker para correr dentro de uma imagem docker com os softwares necessários                                                                       

/t nextflow run main.nf -c nextflow.config -profile docker -resume
