conda activate Numbat

myPath=~
OUTdir=$myPath/result/OUT_Numbat
cd $OUTdir
if [ -f $myPath/VCF_ANNO.list ]; then
  rm $myPath/VCF_ANNO.list
  ls $OUTdir >> $myPath/VCF_ANNO.list
else
  ls $OUTdir >> $myPath/VCF_ANNO.list
fi

ANNOVAR=~/annovar

cat $myPath/VCF_ANNO.list | while read id
do
	echo $id
	vcf=$id/pileup/$id/cellSNP.base.vcf
	avinput=$id/pileup/$id/${id}.avinput
	
	convert2annovar.pl -format vcf4old $vcf > $avinput
	table_annovar.pl $avinput $ANNOVAR/humandb/ -buildver hg38 -out ${id}_anno -protocol clinvar_20221231,refGene,cytoBand,genomicSuperDups,ljb26_all,cosmic70,esp6500siv2_all,1000g2015aug_all -operation f,g,r,r,f,f,f,f -nastring . -csvout -polish
done


for i in *hg38_multianno.csv
do
	sample=`echo $i|awk -F '.' '{print $1}'`
	echo $sample
	cut -f '1-10' $i|sed '1d'|sed "s/$/,${sample}/">> MY_ANNOVAR.csv
done

sed -i '1s/^/Chr,Start,End,Ref,Alt,CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,cytoBand,genomicSuperDups,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,VEST3_score,CADD_raw,CADD_phred,GERP++_RS,phyloP46way_placental,phyloP100way_vertebrate,SiPhy_29way_logOdds,cosmic70,esp6500siv2_all,1000g2015aug_all,Tumor_Sample_Barcode\n/' MY_ANNOVAR.csv


