conda activate Numbat
pileup_and_phase="~/R/x86_64-pc-linux-gnu-library/4.2/numbat/bin/pileup_and_phase.R"
gamp="genetic_map_hg38_withX.txt"
snpvcf='genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf'
hg38_1000g="1000G_hg38"
ncores=16
myPath=~
OUTdir=$myPath/result/OUT_Numbat
mkdir -p $OUTdir
cd $OUTdir

if [ -f $myPath/Numbat.list ]; then
  rm $myPath/Numbat.list
  ls -d $myPath/Brain/*_cellranger710 >> $myPath/Numbat.list
  ls -d $myPath/CSF/*_cellranger710 >> $myPath/Numbat.list
else
  ls -d $myPath/Brain/*_cellranger710 >> $myPath/Numbat.list
  ls -d $myPath/CSF/*_cellranger710 >> $myPath/Numbat.list
fi

cat $myPath/Numbat.list | while read id
do
SAMPLEID=$(echo "$id" | cut -d'/' -f7)
echo $SAMPLEID
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
barcodes=$cellranger_outDir/filtered_feature_bc_matrix/barcodes.tsv.gz
mkdir $SAMPLEID
Rscript $pileup_and_phase \
    --label ${SAMPLEID} \
    --samples ${SAMPLEID} \
    --bams $BAM \
    --barcodes $barcodes \
    --outdir $OUTdir/$SAMPLEID \
    --gmap $gamp \
    --snpvcf $snpvcf \
    --paneldir $hg38_1000g \
    --ncores $ncores 
done

if [ -f $myPath/Numbat.ctc.list ]; then
  rm $myPath/Numbat.ctc.list
  ls -d $myPath/Blood/*_cellranger710 >> $myPath/Numbat.ctc.list
else
  ls -d $myPath/Blood/*_cellranger710 >> $myPath/Numbat.ctc.list
fi

cat $myPath/Numbat.ctc.list | while read id
do
SAMPLEID=$(echo "$id" | cut -d'/' -f7)
echo $SAMPLEID
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
barcodes=$cellranger_outDir/filtered_feature_bc_matrix/barcodes.tsv.gz
mkdir $SAMPLEID
Rscript $pileup_and_phase \
    --label ${SAMPLEID} \
    --samples ${SAMPLEID} \
    --bams $BAM \
    --barcodes $barcodes \
    --outdir $OUTdir/$SAMPLEID \
    --gmap $gamp \
    --snpvcf $snpvcf \
    --paneldir $hg38_1000g \
    --ncores $ncores 
done

if [ -f $myPath/Numbat.PL.list ]; then
  rm $myPath/Numbat.PL.list
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710 >> $myPath/Numbat.PL.list
else
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710 >> $myPath/Numbat.PL.list
fi

cat $myPath/Numbat.PL.list | while read id
do
SAMPLEID=$(echo "$id" | cut -d'/' -f7)
echo $SAMPLEID
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
barcodes=$cellranger_outDir/filtered_feature_bc_matrix/barcodes.tsv.gz
mkdir $SAMPLEID
Rscript $pileup_and_phase \
    --label ${SAMPLEID} \
    --samples ${SAMPLEID} \
    --bams $BAM \
    --barcodes $barcodes \
    --outdir $OUTdir/$SAMPLEID \
    --gmap $gamp \
    --snpvcf $snpvcf \
    --paneldir $hg38_1000g \
    --ncores $ncores 
done 
    
conda activate SComatic
OUTdir=$myPath/result/OUT_SComatic
cd $myPath
SCOMATIC=SComatic
REF=reference/genome/hg38/hg38.fa
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

mkdir -p $OUTdir
cd $OUTdir

if [ -f $myPath/SComatic.list ]; then
  rm $myPath/SComatic.list
  ls -d $myPath/Brain/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/CSF/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/Blood/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710/ >> $myPath/SComatic.list
else
  ls -d $myPath/Brain/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/CSF/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/Blood/*_cellranger710/ >> $myPath/SComatic.list
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710/ >> $myPath/SComatic.list
fi

cat $myPath/SComatic.list | while read id
do
SampleName=$(echo "$id" | cut -d'/' -f7)
echo $SampleName
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
OUT=$OUTdir/$SampleName
mkdir -p $OUT
META=OUT_Numbat/$SampleName/${SampleName}.cell_barcode_annotations.tsv
echo "bash script/SComatic.sh -S" $SCOMATIC "-s" $SampleName "-b" $BAM "-m" $META "-r" $REF "-e" $editing "-p" $PON "-o" $OUT >> mySComatic.sh
done

cat mySComatic.sh | parallel

snpEff=snpEff/snpEff.jar
for id in *_cellranger710
do
echo $id

tsv=${id}/Step4_VariantCalling/${id}.calling.step2.tsv
filted_tsv=${id}/Step4_VariantCalling/cellSNP.filter.vcf
awk -F'\t' '$12 > 9 && $11 > 19' $tsv > $filted_tsv


out=${id}/Step4_VariantCalling/${id}.calling.step2.snpEff.tsv
finalresult=${id}/Step4_VariantCalling/snpEff.anntag.vcf
java -Xmx512G -jar $snpEff hg38 $filted_tsv > $out
perl -alne 'next if $_ =~ /^#/;$F[7] =~ /(ANN=\S+)/;print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$1"' $out > $finalresult
sed -i '1s/^/CHROM\tStart\tEnd\tREF\tALT\tFILTER\tCell_types\tANN\n/' $finalresult
echo `wc -l $finalresult`
done



conda activate gatk
snpEff=snpEff/snpEff.jar
for id in *_cellranger710
do
echo $id
tsv=${id}/pileup/${id}/cellSNP.base.vcf
filted_tsv=${id}/pileup/${id}/cellSNP.filter.vcf
sed 's/AD=//g; s/DP=//g; s/OTH=//g' $tsv > $filted_tsv
sed 's/;/\t/g' $filted_tsv > ${id}/pileup/${id}/cellSNP.filter-1.vcf
awk -F'\t' '$8 > 9 && $9 > 19' ${id}/pileup/${id}/cellSNP.filter-1.vcf > $filted_tsv

out=${id}/pileup/${id}/snpEff.tsv
finalresult=${id}/pileup/${id}/snpEff.anntag.vcf
java -Xmx512G -jar $snpEff hg38 $filted_tsv > $out
perl -alne 'next if $_ =~ /^#/;$F[7] =~ /(ANN=\S+)/;print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$1"' $out > $finalresult
sed -i '1s/^/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tANN\n/' $finalresult
echo `wc -l $finalresult`
done


