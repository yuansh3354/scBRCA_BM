conda activate cpdb 
cd CellphoneDB

file_count=CellphoneDB/counts.txt
file_anno=CellphoneDB/meta.txt
outdir=CellphoneDB/Output
if [ ! -d ${outdir} ]; then
mkdir ${outdir}
fi
cellphonedb method statistical_analysis ${file_anno} ${file_count} --counts-data hgnc_symbol --output-path ${outdir} --threshold 0.01 --threads 64 

