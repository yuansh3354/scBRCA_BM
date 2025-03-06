GENOME_CELLRANGER="~/reference/genome/hg38/10x/refdata-gex-GRCh38-2020-A"
ID="S1"
fq="S1"
SAMPLE="S1"
cellranger count --id "$ID"_cellranger710 \
                 --transcriptome=$GENOME_CELLRANGER \
                 --fastqs=$fq \
                 --localcores=32 \
                 --localmem=512 \
                 --sample=$SAMPLE \
                 --nosecondary
