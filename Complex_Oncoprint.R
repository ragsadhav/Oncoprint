#Include the required libraries
library(ComplexHeatmap)
library(data.table)
library(tibble)

#Read files
mat = fread(file = "Complex oncoprint draw.csv",
              sep = ",",
              header = T, stringsAsFactors = F,
              data.table = F)
mat = as.matrix(column_to_rownames(df = mat, var = "V1"))

mat_ps = fread(file = "Pathway matrix.csv",
              sep = ",",
              header = T, stringsAsFactors = F,
              data.table = F)
mat_ps = as.matrix(column_to_rownames(df = mat_ps, var = "V1"))

mat_ds = fread(file = "Oncogene matrix.csv",
              sep = ",",
              header = T, stringsAsFactors = F,
               data.table = F)
mat_ds = as.matrix(column_to_rownames(df = mat_ds, var = "V1"))

#0/0 location matrix
r <- which(apply(mat,2,rev) == "0/0",arr.ind = TRUE)
mat[mat == "0/0"] = ""

#Oncoprint (first heatmap)
alter_fun = 
  function(x, y, w, h, v) {
    n = sum(v)
    h = h*0.99
    grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.2, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.97, 1/n*h, 
                    gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  }


col = c("missense variant" = "darkgreen","frameshift variant" = "blue","splice acceptor variant" = "purple","splice donor variant" = "orange","inframe deletion" = "red","stop gained" = "black")
n_mutated_patients = apply(mat, 1, function(x) sum(x != ""))
n_mutated_genes = apply(mat, 2, function(x) sum(x != ""))
ht1 = oncoPrint(mat, name = "cases",alter_fun = alter_fun, get_type = function(x) strsplit(x, ";")[[1]],col = col,column_order = NULL,show_column_names = F,column_title = " ",
              row_order = NULL,row_names_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(title = "Alternations", at = c("missense variant","frameshift variant","inframe deletion","splice acceptor variant","splice donor variant","stop gained"),labels = c("missense variant","frameshift variant","inframe deletion","splice acceptor variant","splice donor variant","stop gained")),
              top_annotation = HeatmapAnnotation(barplot = anno_barplot(n_mutated_genes,border =FALSE,gp = gpar(fill = "grey", col = NA)),annotation_height = unit(0.7, "cm")))

#Heatmap for pathways (Second heatmap)
ha_cn2 = HeatmapAnnotation(cn = anno_text(colnames(mat_ps), rot = -60, just = "left",offset = unit(1, "npc") - unit(1, "mm"), gp = gpar(fontsize = 8,fontface = "bold")), annotation_height = unit(6, "cm"))
ht2 = Heatmap(mat_ps, col = c("0" = "white", "1" = "purple"), 
        rect_gp = gpar(col = "grey"), show_row_names = FALSE, cluster_columns = TRUE,
        show_column_dend = FALSE,  show_column_names = FALSE,bottom_annotation = ha_cn2,
        show_heatmap_legend = FALSE, width = unit(5.5, "cm"), column_title = " ")

#Heatmap for annotation (Third heatmap)
ht3 = Heatmap(mat_ds, col = c("0" = "white", "1" = "deeppink","2" = "dodgerblue4", "3" = "limegreen"), 
        rect_gp = gpar(col = "black"), show_row_names = FALSE, cluster_columns = FALSE,
        show_column_dend = FALSE,  show_column_names = FALSE,heatmap_legend_param = list(title = "Annotation",at = c("1","2","3"),labels = c("Oncogene","Tumor suppressor","Cancer gene census")),
        show_heatmap_legend = T, width = unit(0.5, "cm"), column_title = " ")


#Draw the heatmap in a file
png(filename="plot.png", width = 10000, height = 10000, res = 900)

#Display the heatmap in plots section
draw(ht1+ht2+ht3, padding = unit(c(1, 1, 1, 1), "cm"))

#Add the gridlines as the groups
decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "305_MMG120")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("VMMG", x, unit(1, "npc") + unit(1, "mm"),hjust=2,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})

decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "1050P_MMG41")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("PMMG", x, unit(1, "npc") + unit(1, "mm"),hjust=1.5,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})

decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "476_MMG35")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("476MMG", x, unit(1, "npc") + unit(1, "mm"),hjust=1,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})

decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "476_T37")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("476PT", x, unit(1, "npc") + unit(1, "mm"),hjust=2,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})

decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "153_TDH32")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("153PT", x, unit(1, "npc") + unit(1, "mm"),hjust=2,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})
decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "153_LMTDH31")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 3))
  grid.text("153LMT", x, unit(1, "npc") + unit(1, "mm"),hjust=2,vjust = -5,gp=gpar(fontsize=6.5,fontface = "bold"))
})
for (a in 1:nrow(r)){
  rw = r[a,2]
  cl = r[a,1]
  decorate_heatmap_body("cases", {
    grid.rect((rw-0.5)/ncol(mat), (cl-0.5)/nrow(mat),1/ncol(mat),1/nrow(mat), default.units = "npc",gp = gpar(fill = "grey89",col = "white"))
  })
}



dev.off()





