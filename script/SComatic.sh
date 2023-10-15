
while getopts "S:s:b:m:r:e:p:o:" opt; do
  case ${opt} in
    S) SCOMATIC=$OPTARG ;;
    s) sample=$OPTARG ;;
    b) BAM=$OPTARG ;;
    m) META=$OPTARG ;;
    r) REF=$OPTARG ;;
    e) editing=$OPTARG ;;
    p) PON=$OPTARG ;;
    o) output_dir=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2 ; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
  esac
done

output_dir1=$output_dir/Step1_BamCellTypes
output_dir2=$output_dir/Step2_BaseCellCounts
output_dir3=$output_dir/Step3_BaseCellCountsMerged
output_dir4=$output_dir/Step4_VariantCalling
output_dir5=$output_dir/Step5_CellTypeCallableSites
output_dir6=$output_dir/Step6_UniqueCellCallableSites
output_dir7=$output_dir/Step7_SingleCellAlleles
STEP4_1=$output_dir4/${sample}.calling.step1.tsv
STEP4_2_pass=${output_dir4}/${sample}.calling.step2.pass.tsv

mkdir -p $output_dir $output_dir1 $output_dir2 $output_dir3 $output_dir4 $output_dir5 $output_dir6 $output_dir7

if [ ! -d "temp_code" ]; then
  mkdir temp_code
fi

echo "Step1: Splitting alignment file in cell type specific bams `date`"
python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $BAM \
        --meta $META \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1
echo "Step 1 ends at `date`"
echo "Step2: Collecting base count information `date`"
for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 32

  rm -rf $temp
done
echo "Step2 ends at `date`"
echo "Step3: Merging base count matrices `date`"
python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv

echo "Step2 ends at `date`"

echo "Step4: Detection of somatic mutations `date`"
python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF
python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON

bedtools intersect -header -a ${output_dir4}/${sample}.calling.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${sample}.calling.step2.pass.tsv
echo "Step4 ends at `date`"