dir1="G:/【结题报告】KJ-DB-PM-BJ161062-03-2018-P0/upload/4_de/4_2_singleDE_table"
setwd(dir1)
list.files()
file1=read.csv(file = "T1_C1.anno.csv",stringsAsFactors = F,sep = '\t')
colnames(file1)[which(grepl(pattern = "ENSG.*",x = file1[10,],perl = TRUE,ignore.case = T))]='Gene'
colnames(file1)[which(grepl(pattern = "^Log2FoldChange$",x = colnames(file1),ignore.case = T,perl = T))]=paste("T1C1","Log2FoldChange",sep = "_")
file1_use=file1[,c("Gene","T1_count","T1_normalize","C1_count","C1_normalize","T1C1_Log2FoldChange","GeneName")]
file1_use[,"GeneGeneName"]=paste(file1_use$Gene,file1_use$GeneName,sep = "_")
file1_use=file1_use[apply(X = file1_use,MARGIN = 1,FUN = function(x){sum(as.numeric(x[c(2,4)])<10)<2}),]
#file1_use=file1_use[apply(X = file1_use,MARGIN = 1,FUN = function(x){as.numeric(x[6])>0.5849625}),]


dir2="G:/【结题报告】KJ-DB-PM-BJ161062-04-2018-P0/upload/4_de/4_2_singleDE_table"
setwd(dir2)
list.files()
file2=read.csv(file = "T2_C2.anno.csv",stringsAsFactors = F)
colnames(file2)[which(grepl(pattern = "ENSG.*",x = file2[10,],perl = TRUE,ignore.case = T))]='Gene'
colnames(file2)[which(grepl(pattern = "^Log2FoldChange$",x = colnames(file2),ignore.case = T,perl = T))]=paste("T2C2","Log2FoldChange",sep = "_")
file2_use=file2[,c("Gene","T2_count","T2_normalize","C2_count","C2_normalize","T2C2_Log2FoldChange","GeneName")]
file2_use[,"GeneGeneName"]=paste(file2_use$Gene,file2_use$GeneName,sep = "_")
file2_use=file2_use[apply(X = file2_use,MARGIN = 1,FUN = function(x){sum(as.numeric(x[c(2,4)])<10)<2}),]
#file2_use=file2_use[apply(X = file2_use,MARGIN = 1,FUN = function(x){as.numeric(x[6])>0.5849625}),]



files_use=merge(x = file1_use,y=file2_use,by.x = 'GeneGeneName',by.y='GeneGeneName',all = T)
rownames(files_use)=files_use$GeneGeneName
data_use=files_use[which(apply(X = files_use,MARGIN = 1,FUN = function(x){max(abs(as.numeric(x[c(7,14)])))>=0.5849625})),c("C1_normalize","T1_normalize","C2_normalize","T2_normalize")]
###NA
for (rowIndex in c(1:nrow(data_use))) {
  for(colIndex in c(1:ncol(data_use))){
    if(is.na(data_use[rowIndex,colIndex])){
      data_use[rowIndex,colIndex]=0
    }
  }
}



# library(NGCHM)
colnames(data_use)=stringi::stri_extract(str = colnames(data_use),regex = ".*(?<=_)")
mat=as.matrix(data_use)
mat=t(scale(t(mat),center = T,scale = T))
hm <- chmNew ('my-first-ngchm_YX',mat)
chmExportToPDF(chm = hm1, filename = "temp.pdf",shaidyMapGen = "F:/min-labs_paper/work/WHX/ShaidyMapGen.jar",overwrite = T)
chmExportToFile(hm1,paste('YANGXIANG','.ngchm',sep = "_"),shaidyMapGen = "F:/min-labs_paper/work/YANGXIANG/ShaidyMapGen.jar",overwrite = T)

BREAKS_UP=seq(median(mat),quantile(unlist(mat),0.90),length.out = 50)
(COL_UP=colorRampPalette(c("black","red"))(length(BREAKS_UP)-1))
#(COL_UP=colorRampPalette(c("white","red"))(length(BREAKS_UP)))
(BREAKS_DOWN=seq(quantile(mat,0.10),median(mat),length.out = 50))
(COL_DOWN=colorRampPalette(c("green","black"))(length(BREAKS_DOWN)))

(BREAKS_DOWN_UP=unique(c(BREAKS_DOWN,BREAKS_UP)))
(COL_DOWN_UP=unique(c(COL_DOWN,COL_UP)))
message(paste("COL",length(COL_DOWN_UP),"BREAKS",length(BREAKS_DOWN_UP),sep = "_"))
#####默认的方法只是改改颜色
# cmap1 <- chmNewColorMap (BREAKS_DOWN_UP,COL_DOWN_UP)
# layer1 <- chmNewDataLayer ('layer.name', mat, cmap1)
# hm <- chmNew (paste('WHXT123vsC1',sep = "_"),layer1)
# chmExportToFile(hm,paste('WHX','.ngchm',sep = "_"),shaidyMapGen = "F:/min-labs_paper/work/YANGXIANG/ShaidyMapGen.jar",overwrite = T)
# chmExportToPDF(hm,filename = paste('YANGXIANG','.pdf',sep = "_"),shaidyMapGen = "F:/min-labs_paper/work/YANGXIANG/ShaidyMapGen.jar",overwrite = T)
# 

library(gplots)
tiff(filename = "temp.tiff",height = 1800,width = 800,res = 90*3,type = "windows")
svg(filename = "temp.svg",width = 10,height = 10)
win.graph();gplots::heatmap.2(x = mat,scale ='none',breaks = BREAKS_DOWN_UP,col = COL_DOWN_UP
                              ,trace="none"
                              ,hclustfun = function(x){hclust(x,method = 'ward.D2')}
                              #,distfun = function(x){dist(x,method = "manhattan")}
                              ,density.info = "none"
                              ,key = FALSE
                  ,Colv = FALSE
                  ,labRow = FALSE
                  ,labCol = FALSE
                  )

dev.off()

####画热图使用pvalue sign(log2foldchage)
file1_use_pvalue=file1[,c("Gene","GeneName","padj","T1C1_Log2FoldChange")]
file1_use_pvalue[,"T1C1_signPvalue"]=-log10(file1_use_pvalue$padj)*sign(file1_use_pvalue$T1C1_Log2FoldChange)
file1_use_pvalue[,"GeneGeneName"]=paste(file1_use_pvalue$Gene,file1_use_pvalue$GeneName,sep = "_")

file2_use_pvalue=file2[,c("Gene","GeneName","padj","T2C2_Log2FoldChange")]
file2_use_pvalue[,"T2C2_signPvalue"]=-log10(file2_use_pvalue$padj)*sign(file2_use_pvalue$T2C2_Log2FoldChange)
file2_use_pvalue[,"GeneGeneName"]=paste(file2_use_pvalue$Gene,file2_use_pvalue$GeneName,sep = "_")
files_use_merge=merge(x = file1_use_pvalue,y=file2_use_pvalue,by.x = 'GeneGeneName',by.y='GeneGeneName')

data_use_merge=files_use_merge[,c("GeneGeneName","T1C1_signPvalue","T2C2_signPvalue")]
rownames(data_use_merge)=data_use_merge$GeneGeneName
data_use_merge=data_use_merge[,-1]

data_use_merge=data_use_merge[order(data_use_merge$T2C2_signPvalue,data_use_merge$T1C1_signPvalue),]
hm <- chmNew ('YX_heatmap_pvalue',as.matrix(data_use_merge),rowOrder = NULL)
chmExportToPDF(chm = hm, filename = "YX_heatmap_signpvalue.pdf",shaidyMapGen = "F:/min-labs_paper/work/WHX/ShaidyMapGen.jar",overwrite = T)
chmExportToFile(hm,paste('YX_heatmap_signvalue','.ngchm',sep = "_"),shaidyMapGen = "F:/min-labs_paper/work/YANGXIANG/ShaidyMapGen.jar",overwrite = T)
