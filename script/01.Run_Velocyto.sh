conda activate Velocyto
myPath=~
OUTdir=$myPath/result/OUT_Velocyto
mkdir -p $OUTdir
cd $OUTdir

rmsk_gtf=hg38_repeat_rmsk.gtf 
cellranger_gtf=refdata-gex-GRCh38-2020-A/genes/genes.gtf 


if [ -f $myPath/Velocyto.list ]; then
  rm $myPath/Velocyto.list
  ls -d $myPath/Brain/*_cellranger710/ >> $myPath/Velocyto.list
  ls -d $myPath/CSF/*_cellranger710/ >> $myPath/Velocyto.list
else
  ls -d $myPath/Brain/*_cellranger710/ >> $myPath/Velocyto.list
  ls -d $myPath/CSF/*_cellranger710/ >> $myPath/Velocyto.list
fi

cat $myPath/Velocyto.list | while read id
do
SampleName=$(echo "$id" | cut -d'/' -f7)
echo $SampleName
  time velocyto run10x -@ 64 -m $rmsk_gtf $id $cellranger_gtf 
done

if [ -f $myPath/Velocyto.ctc.list ]; then
  rm $myPath/Velocyto.ctc.list
  ls -d $myPath/Blood/*_cellranger710/ >> $myPath/Velocyto.ctc.list
else
  ls -d $myPath/Blood/*_cellranger710/ >> $myPath/Velocyto.ctc.list
fi

cat $myPath/Velocyto.ctc.list | while read id
do
SampleName=$(echo "$id" | cut -d'/' -f7)
echo $SampleName
  time velocyto run10x -@ 64 -m $rmsk_gtf $id $cellranger_gtf 
done

if [ -f $myPath/Velocyto.PL.list ]; then
  rm $myPath/Velocyto.PL.list
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710/ >> $myPath/Velocyto.PL.list
else
  ls -d $myPath/Primary_and_lymphnode/*_cellranger710/ >> $myPath/Velocyto.PL.list
fi

cat $myPath/Velocyto.PL.list | while read id
do
SampleName=$(echo "$id" | cut -d'/' -f7)
echo $SampleName
  time velocyto run10x -@ 64 -m $rmsk_gtf $id $cellranger_gtf 
done
