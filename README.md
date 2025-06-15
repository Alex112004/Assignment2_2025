**Versões dos softwares utilizados:**
    
    FastQC v0.12.1    
    multiqc, version 1.29
    nextflow version 25.04.3.5949
    R version 4.4.3 (2025-02-28) -- "Trophy Case"

**Codigos a correr:**

Atualizar a vm:

    sudo apt update
    
Ativar o micromamba/conda:

    micromamba activate / conda activate

Instalar o fastqc:

    micromamba install  fastq / conda install  fastqc

Instalar o multiqc:

    micromamba install multiqc / conda install multiqc

Instalar o docker:

    sudo apt install docker.io docker-compose
    sudo usermod -aG docker $USER
    
Código para correr o fastqc:
    
    fastqc *.gz
    
Código para correr o multiqc:
    
    multiqc ./

instalar o Silva:

    mkdir ~/databases
    cd ~/databases
    wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz
    tar xfv silva.seed_v138_1.tgz
    rm silva.seed_v138_1.tgz

Clunar o repositorio do metaflowmics:

    git clone https://github.com/hawaiidatascience/metaflowmics.git

Abir a pasta da pipeline do 16S:

    cd metaflowmics/metaflowmics/Pipeline-16S/

Correr nextflow com as novvas configurações e utilizar o comando -profile docker para correr dentro de uma imagem docker com os softwares necessários:

    nextflow run main.nf -c nextflow.config -profile docker -resume

Arquivo diversidade alfa:

    cd ~/drought_results/postprocessing cut -f 2 alpha-diversity.97.summary | cut -d "_" -f 1 |sed 's/group/type/' | paste alpha-diversity.97.summary - > alpha-div.tsv

**Colaboradores:**
    
- Alexandre Soares
    
- Sara Sampaio
    
- Maria Eduarda Oiveira
    
- Biatrix Venâncio
