#10:45 2025/5/28 
##$$$$$$$$$$$$$$$ Example of fastq_list $$$$$$$$$$$$$$$$$$$$$$$$$#
# /data5/Wuky/4CHIATAC10G/3YOfinal/yHSC-5 yHSC-5
# Command：nohup bash CHIATAC.sh config_CPUatac_mm10_l96.sh fastq_list11.txt > chiatac11.log 2>&1 &
# tested by wky, 2025/5/28

#!/bin/bash

# 2025.03.17 ljh@julab修改自tjongh@JGM GT
# 主要修改：
# 1. 由于julab没有集群，将有关sulrm的部分都修改掉了
# 2. 为了适配julab环境，做了一些容器路径的修改
# 2021.02.05 原文件注释
# updated with rm blacklist  2021.04.19 (avail only for hg38 )
# hichipper is NOT yet used for peak calling
# Backbone: ChIA-PET Utilities (by Chee Hong Wong https://github.com/cheehongsg/CPU/wiki), macs2, samtools, bedtools, kent, juicer, and powerful BASH :-)
# The Jackson Laboratory for Genomic Medicines
# Genome Ref is not included in the container so this is a lite version and the run must have the ref copy step.

# fastq_list has 2 columns:
# prefix with full path of fastq files & RUN_ID

#$$$$$$$$$$$$$$$ Example of fastq_list $$$$$$$$$$$$$$$$$$$$$$$$$#
# /projects/wei-lab/proj-in-situ-chia-pet/silvia/WT2_GT18-06298_AGGCAGAA-AGAGGATA_S2 wt2
# /projects/wei-lab/proj-in-situ-chia-pet/silvia/MUT1_GT18-06300_TAAGGCGA-ATAGAGAG_S4 mut1

# tested by linjiahao, 2025/03/17

# 设置输入参数
configfile=$1 #configuration general to a set of libraries
fastq_list=$2 #fastq_file with RUN id list

#####
## 检查输入参数
## 如果参数数量不大于1，则输出错误信息并退出脚本
die() {
  echo >&2 "$@"
  exit 1
}
[ $# -gt 1 ] || die "config file and fastq_list in this order required!"

## 获取当前工作目录，并将其赋值给变量 maindir
maindir=$PWD

## @date 2025.03.17 ljh@julab
## 因为bedtools使用的是singularity环境
## 所以需要先加载singularity环境
# conda activate singularity

## 引入配置文件
source $configfile

## 因为singularity run 默认工作目录是home，所以需要手动把当前目录绑定到容器的home目录下
## 不然运行命令时，容器会找不到文件
## 同时绑定参考基因组的目录，否则bedtools工具无法找到参考基因组文件
scom="singularity run -B $PWD:$PWD -B $refdir:$refdir $addpath/cpu0.0.1a-r2.sumner.sif"

## 打包软件
pigz="$scom pigz"
## samtools
samtools="$scom samtools"
## bedtools
bedtools="$scom bedtools"
# bedtools="bedtools"
## ChiaSig
chiasigprog="$scom ChiaSig"
## juice容器，并绑定当前目录到容器的当前目录下

####这里做了小修改，去掉sudo指令@date 2025.05.28 wky
#juicertool="sudo -S $singularityPath run -B $PWD:$PWD  $addpath/juicer_1.22.01.sif"
juicertool="$singularityPath run -B $PWD:$PWD  $addpath/juicer_1.22.01.sif"

## 定义变量 (除了参考基因组目录，其他变量都在容器cpu0.0.1a-r2.sumner.sif内)
## 参考基因组目录
fasta="$refdir/$genome/bwa/$genome.fa"
## 注释文件（暂时不清楚是什么注释文件）
annot=/opt/gt/annotation/${genome}.Genes.itxAnnotation.bed
## cpu0.0.1a-r2.sumner.sif容器中bedtools的位置，应该位于启动服务器时绑定的目录下
blackbed=${refdir}${genome}/bwa/$genome.blacklist.bed

## 注释perl脚本
annotPerl=/opt/gt/scripts/itxAnnotation.pl
## 提取摘要的脚本路径
extractsummary=/opt/gt/scripts/extract_summary_cpustat.sh
## 提取矩阵的脚本路径（暂时不清楚）
matrixprog=/opt/gt/scripts/clusters2matrix.$genome.py
## 暂时不清楚
Rciftif=/opt/gt/scripts/printCifTif.R

## 定义生成结果后缀#defined names
pairlabel='singlelinker.paired' ## 配对数据（有互作，ATAC开放）
singlabel='singlelinker.single' ## 无配对数据（无互作，ATAC开放）
nonelabel='none'                ## 不匹配linker的数据（无互作，ATAC不开放）

## 自定义创建结果目录的函数
function odir {
  outputDir=$1
  if [ ! -d $outputDir ]; then
    echo "Create $outputDir"
    mkdir $outputDir
  fi
}

# 定义一个函数，用于按顺序提交作业（nohup）
function pbsjob {

  # 函数接受两个参数：fasta文件路径，不包括后缀，和RUN作业名称
  fqfile=$1
  RUN=$2

  ## 结果文件名称定义
  read1=${fqfile}${R1suf} ## fastq文件名称定义
  read2=${fqfile}${R2suf}
  pbs="$RUN.insituChiapet.sh"                            ## 作业文件名称定义
  peakfile="${RUN}_peaks.narrowPeak"                     ## peak文件名称定义
  cis_prefix="${RUN}.e${extbp}.clusters.cis.chiasig"     ## cis文件名称定义
  trans_prefix="${RUN}.e${extbp}.clusters.trans.chiasig" ## trans文件名称定义
  chiasigoutput="${cis_prefix}.BE2.sigf.interactions"    ## chiasig输出文件名称定义
  shortname=${RUN}.chiapet
  ## 创建结果目录
  odir $RUN
  cd $RUN

  cat <<ENDHERE >$pbs

#!/bin/bash

echo "Workdir: \$PWD"
LOGFILE=${shortname}.log
export genome=$genome

##  复制util.sh脚本到当前目录下
# $scom cp /opt/gt/scripts/util.sh $PWD/util.sh
# sed 's/\/fastscratch\/ref_chiapet\//\/data\/home\/ruanlab\/huangjiaxiang\/pipeline\/ChIATAC\//' util.sh > util2.sh
# source util2.sh



echo "Reference genome:"
/bin/ls ${fasta}*
# /bin/rm util.sh util2.sh


## 执行linker检测并生成不同类别的fastq文件
echo "Linker detection on: ${RUN} " 2>\${LOGFILE}

## 检测linker
$cpuprog stag -A ${linker} -W -T 18 -t ${NTHREAD} -O ${RUN} $read1 $read2  2>>\${LOGFILE}
echo "--- linker detection completed ----" >>\${LOGFILE}

## 统计linker检测结果
$cpuprog stat -s -p -T 18 -t ${NTHREAD} ${RUN}.cpu 2>>\${LOGFILE} 1>${RUN}.stat
echo "--- statistics done  ----" >>\${LOGFILE}

echo "--- pigziiping   ----" >>\${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.singlelinker.paired.fastq 2>>\${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.singlelinker.single.fastq 2>>\${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.none.fastq 2>>\${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.conflict.fastq 2>>\${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.tied.fastq 2>>\${LOGFILE}
echo "--- pigz completed ----" >>\${LOGFILE}


#第二步骤Pair双端匹配
#-- perform hybrid bwa men and bwa all mapping, de-duplication, span computation, and tag clustering --#
# mapping
echo  START  ${RUN} cpu memaln .. | tee -a \${LOGFILE}
echo  Mapping paired tags .. | tee -a \${LOGFILE}
$cpuprog memaln -T $minMapScore -t ${NTHREAD} $fasta ${RUN}.$pairlabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$pairlabel.sam
$pigz -p ${NTHREAD} ${RUN}.$pairlabel.sam | tee -a \${LOGFILE}
echo  ENDED pair mapping | tee -a \${LOGFILE}

#pairing
echo  STARTED ${RUN} cpu pair .. | tee -a \${LOGFILE}
echo  Pairing paired tags .. | tee -a \${LOGFILE}
$cpuprog pair -s $selfbp -S -t ${NTHREAD} ${RUN}.$pairlabel.sam.gz 1>${RUN}.$pairlabel.stat.xls 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a \${LOGFILE}

# span
echo  STARTED ${RUN} cpu span .. | tee -a \${LOGFILE}
echo  Computing span of paired tags .. | tee -a \${LOGFILE}
$cpuprog span -s $selfbp -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 2>> \${LOGFILE} 1>${RUN}.$pairlabel.UU.span.xls
echo  ENDED ${RUN} span pair .. | tee -a \${LOGFILE}

# deduplication
echo  STARTED ${RUN} cpu dedup .. | tee -a \${LOGFILE}
echo  De-duplicating paired tags UU .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 1>${RUN}.$pairlabel.UU.dedup.lc 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu dedup .. | tee -a \${LOGFILE}

# deduplicated span
echo  STARTED ${RUN} cpu dedup span.. | tee -a \${LOGFILE}
echo  Computing span of paired tags UU nr .. | tee -a \${LOGFILE}
$cpuprog span -s $selfbp -t ${NTHREAD} ${RUN}.$pairlabel.UU.nr.bam 2>> \${LOGFILE} 1>${RUN}.$pairlabel.UU.nr.span.xls
echo  ENDED ${RUN} cpu dedup span.. | tee -a \${LOGFILE}

# cluster tags
echo  STARTED ${RUN} cpu clustering.. | tee -a \${LOGFILE}
$cpuprog cluster -s $selfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${RUN}.e$extbp -j -x -v 1 -g ${RUN}.$pairlabel.UU.nr.bam 2>&1 | tee -a \${LOGFILE}
echo  ENDED ${RUN} cpu clustering.. | tee -a \${LOGFILE}


#没有tag的匹配------------------------------------------------------------------------------
#None tag

# mapping
echo  STARTED ${RUN}.$nonelabel cpu memaln .. | tee -a \${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$nonelabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$nonelabel.sam
$pigz -p ${NTHREAD} ${RUN}.$nonelabel.sam 2>&1  | tee -a \${LOGFILE}
echo  ENDED ${RUN} cpu memaln .. | tee -a \${LOGFILE}

# pairing
echo Pairing $nonelabel tags .. | tee -a \${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$nonelabel.sam.gz 1>${RUN}.$nonelabel.stat.xls 2>> \${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a \${LOGFILE}

# span
echo  STARTED ${RUN}.$nonelabel cpu span .. | tee -a \${LOGFILE}
$cpuprog span -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 2>> \${LOGFILE} 1>${RUN}.$nonelabel.UU.span.xls
echo  ENDED ${RUN}.$nonelabel span pair .. | tee -a \${LOGFILE}

# deduplication
echo  STARTED ${RUN}.$nonelabel cpu dedup .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 1>${RUN}.$nonelabel.UU.dedup.lc 2>> \${LOGFILE}


#单端匹配------------------------------------------------------------------------------

#1tag

# mapping
echo STARTED ${RUN}.$singlabel cpu memaln .. | tee -a \${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$singlabel.fastq.gz 2>> \${LOGFILE} 1>${RUN}.$singlabel.sam
$pigz -p ${NTHREAD} ${RUN}.$singlabel.sam 2>&1  | tee -a \${LOGFILE}

# pairing to get bam
echo Pairing $singlabel tags .. | tee -a \${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$singlabel.sam.gz 1>${RUN}.$singlabel.stat.xls 2>> \${LOGFILE}

# deduplication
echo STARTED ${RUN}.$singlabel cpu dedup .. | tee -a \${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$singlabel.UxxU.bam 1>${RUN}.$singlabel.UxxU.dedup.lc 2>> \${LOGFILE}
echo  ENDED ${RUN}.$singlabel cpu dedup .. | tee -a \${LOGFILE}


#第三步骤生成高质量的bam------------------------------------------------------------------------------
#Filter out non primary reads
## 对于ATAC和CHIA-pet分析，主要关注其中的主要比对（高质量比对）
## 所以这里需要过滤掉非主要比对（低质量）
## 2048是bam文件中过的一个标志位，代表非主要比对（non primary reads）
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$pairlabel.UU.nr.bam  | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$pairlabel.F2048.bam
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$singlabel.UxxU.nr.bam | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$singlabel.F2048.bam
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$nonelabel.UU.nr.bam  | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$nonelabel.F2048.bam

echo -e "Converting file formats..\n" >> \${LOGFILE}
$samtools sort -@ $NTHREAD -o ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$pairlabel.F2048.bam
$samtools sort -@ $NTHREAD -o ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$singlabel.F2048.bam    
$samtools sort -@ $NTHREAD -o ${RUN}.$nonelabel.nr.sorted.bam ${RUN}.$nonelabel.F2048.bam  

#接下来做下游的peakcalling，forBASIC still have blacklist.  Removed secondary alignment with -F2048; for clean view we have NDP ${RUN}_treat_pileup.clip.sorted.bdg
$samtools merge -@ $NTHREAD ${RUN}.forBASIC.bam ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$nonelabel.nr.sorted.bam

echo 'Convert bam to bed for macs2 callpeak and clean from blacklist... '
$bedtools bamtobed -i ${RUN}.$pairlabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$pairlabel.nr.sorted.bed
$bedtools bamtobed -i ${RUN}.$singlabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$singlabel.nr.sorted.bed
$bedtools bamtobed -i ${RUN}.$nonelabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$nonelabel.nr.sorted.bed
ls $RUN.*.nr.sorted.bed

$scom macs2 callpeak --keep-dup all --nomodel --shift $shiftsize  --extsize $peakext  -B --SPMR -t ${RUN}.$pairlabel.nr.sorted.bed ${RUN}.$singlabel.nr.sorted.bed ${RUN}.$nonelabel.nr.sorted.bed  -f BED -g $organism -n $RUN  --qvalue $macsq >> \${LOGFILE}

echo  Generating coverage density.. | tee -a \${LOGFILE}
#$将结果变为可视化的bw文件，samtools view -H ${RUN}.$pairlabel.UU.nr.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.genome.length
$samtools view -H ${RUN}.forBASIC.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.genome.length

$scom bedClip ${RUN}_treat_pileup.bdg ${RUN}.genome.length ${RUN}_treat_pileup.clip.bdg
LC_COLLATE=C sort -k1,1 -k2,2n ${RUN}_treat_pileup.clip.bdg > ${RUN}_treat_pileup.clip.sorted.bdg

$scom bedGraphToBigWig ${RUN}_treat_pileup.clip.sorted.bdg ${RUN}.genome.length ${RUN}.treat_pileup.NDP.bw
echo ENDED ${RUN} coverage density generated. | tee -a \${LOGFILE}

# Make bedgraph
$bedtools genomecov -ibam ${RUN}.forBASIC.bam -bg > ${RUN}.forBASIC.bedgraph

# Sort bedgraph
$scom bedSort ${RUN}.forBASIC.bedgraph ${RUN}.forBASIC.sorted.bedgraph

# Make bigwig
$scom bedGraphToBigWig ${RUN}.forBASIC.sorted.bedgraph  ${RUN}.genome.length ${RUN}.coverage.forBASIC.bw
echo -e  "Done converting file formats, ${RUN}.coverage.forBASIC.bw generated.\n" >> \${LOGFILE}


#最后一步输出可视化的文件.hic等等------------------------------------------------------------------------------
#Significant interaction
#get ipet > 1
echo  Creating BE2 file and clean blacklist using $blackbed  >> \${LOGFILE}
zcat ${cis_prefix}.gz | awk '{ if ( \$7 > 1 ) print }' | $bedtools pairtobed -type neither -a stdin -b $blackbed > ${cis_prefix}.BE2

#run chiasig
$chiasigprog -c $PET -t $NTHREAD -p -m $selfbp ${cis_prefix}.BE2

#Annotation
echo "Annotating..."
$scom perl $annotPerl  $chiasigoutput $annot $peakfile $RUN | tee -a \${LOGFILE}

#dump out itx
echo  Creating itx files  | tee -a \${LOGFILE}
$scom Rscript /opt/gt//scripts/itx_from_annotated.R $RUN | tee -a \${LOGFILE}


#Get the summary
$scom $extractsummary --conf $maindir/$configfile -r $RUN


#generate matrix
echo "Generating matrix"
cat ${cis_prefix}.gz ${trans_prefix}.gz > $RUN.pooled.clusters.gz
$scom python $matrixprog ${RUN}.pooled.clusters.gz $RUN $binsize


#TIF CIF
echo "Calculating normalized matrix and print TIF & CIF"
$scom Rscript $Rciftif $RUN $binsize

## Juicebox可视化，生成hic文件
# Heatmap for Juicebox
#sed 's/chr//' ${RUN}.genome.length | sort -k1,1 -V | grep -v '^M' >  ${RUN}.nochr.genome.length #remove chrM
sed 's/chr//' ${RUN}.genome.length | sed 's/^M/ZZZZ/' | sort -k1,1 -V | sed 's/^ZZZZ/M/' >  ${RUN}.nochr.genome.length
echo ${sudopw} | $juicertool pre --threads $NTHREAD -r 2500000,1000000,500000,250000,100000,50000,25000,10000 -k KR,GW_KR ${RUN}.e${extbp}.juice.gz $RUN.hic  ${RUN}.nochr.genome.length  2>> \${LOGFILE}

$pigz -p ${NTHREAD} $RUN.*.bedgraph
$pigz -p ${NTHREAD} $RUN.*.bdg #keep the clip.sorted.bdg

echo "CPU pipeline run $RUN  completed!"
ENDHERE
  chmod +x $pbs
  ./$pbs
  cd $maindir
}

## 按顺序提交作业
while read myarg; do
  pbsjob $myarg
done <$fastq_list
