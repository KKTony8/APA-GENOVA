#### 第一步： 生成Accept和Reject diffTAD的两个文件
## 21:11 2025/6/9 wky
hicDifferentialTAD \
  --targetMatrix /data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/oHSC-8.mcool::resolutions/25000 \
  --controlMatrix /data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/yHSC-A.mcool::resolutions/25000 \
  --tadDomains /data5/Wuky/CHIFQ/software/YOResult/hic60/TADnew2/oHSC_TADs_domains.bed \
  --outFileNamePrefix /data5/Wuky/CHIFQ/software/YOResult/hic60/TADnew2/TADdiff/oHSC_vs_yHSC_diffTAD \
  --mode all \
  --modeReject one \
  --pValue 0.05 \
  --threads 8


#### 第二步： 先建立好track.init文件，然后用pyGenomeTracks来绘制轨道图，任意选取染色体位置绘制
## 15:31 2025/6/9 wky

pyGenomeTracks   --tracks tracks.ini   --region 1:4250000-20000000   --outFileName yHSC_vs_oHSC_diffTAD_region.pdf   --width 100   --dpi 300

