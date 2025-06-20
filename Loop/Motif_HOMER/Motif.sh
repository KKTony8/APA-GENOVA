21:05 2025/5/26
###1. 下载和安装HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl --install

#检查是否下载安装成功
findMotifsGenome.pl -h

###2.loop的bedpe文件需要转换为bed文件，并且去除重叠区域
#首先把空格分隔转换成TAB分隔文件(即^I分隔)
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6}' old.loop.bedpe > old.loop.tab.bedpe
awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6}' young.loop.bedpe > young.loop.tab.bedpe

#提取bedpe文件
cut -f1-3 old.loop.tab.bedpe > old.anchor1.bed
cut -f4-6 old.loop.tab.bedpe > old.anchor2.bed
cut -f1-3 young.loop.tab.bedpe > young.anchor1.bed
cut -f4-6 young.loop.tab.bedpe > young.anchor2.bed

#合并bed文件
cat old.anchor1.bed old.anchor2.bed | sort -k1,1 -k2,2n | bedtools merge > old.anchors.bed
cat young.anchor1.bed young.anchor2.bed | sort -k1,1 -k2,2n | bedtools merge > young.anchors.bed

###3.Old的差异loop文件做motif分析
findMotifsGenome.pl \
/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/old.anchors.bed \
/data5/Wuky/CHIFQ/ref/mm10/bwa/mm10.fa \
/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/old_motifs_output \
-size given

###4.Young的差异loop做motif分析
findMotifsGenome.pl \
/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/young.anchors.bed \
/data5/Wuky/CHIFQ/ref/mm10/bwa/mm10.fa \
/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/young_motifs_output \
-size given
