**Forma de utilizar:**


Atualizar a vm:

    sudo apt update

Ativar o micromamba/conda:

    micromamba activate / conda activate

Instalar o :

    git clone https://github.com/hawaiidatascience/metaflowmics.git

Abir a pasta da pipeline do 16S:

    cd metaflowmics/metaflowmics/Pipeline-16S/

Correr nextflow com as novvas configurações e utilizar o comando -profile docker para correr dentro de uma imagem docker com os softwares necessários:

    nextflow run main.nf -c nextflow.config -profile docker -resume

**Colaboradores:**
    Alexandre Soares
    Sara Sampaio
    Maria Eduarda Oiveira
    Biatrix Venâncio
