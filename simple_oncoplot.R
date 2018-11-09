library(ComplexHeatmap)
library(data.table)
library(tibble)

mat <-  fread(file = "153_ONCOPRINT_DRAW.csv",sep = ",",header = T, stringsAsFactors = F,data.table = F)
mat = as.matrix(column_to_rownames(df = mat, var = "V1"))
r <- which(apply(mat,2,rev) == "0/0",arr.ind = TRUE)
mat[mat == "0/0"] = ""

alter_fun = 
  function(x, y, w, h, v) {
  n = sum(v)
  h = h*0.99
  grid.rect(x, y, w-unit(0.2, "mm"), h-unit(0.2, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.97, 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  }


col = c("missense variant" = "darkgreen","disruptive inframe deletion" = "blue", "frameshift variant" = "orange","splice acceptor variant" = "red","splice donor variant" = "purple","inframe deletion" = "pink", "inframe insertion" = "yellow", "stop gained" = "black")
#png(filename="plot.png", width = 20000, height = 10000, res = 1000)
h1 = oncoPrint(mat, name = "mat",alter_fun = alter_fun, get_type = function(x) strsplit(x, ";")[[1]],col = col,column_order = NULL,show_column_names = T,column_title = " ",
          row_order = NULL,row_names_gp = gpar(fontsize = 12),
          heatmap_legend_param = list(title = "Alternations", at = c("missense variant","frameshift variant", "disruptive inframe deletion","splice acceptor variant","splice donor variant","inframe deletion","inframe insertion", "stop gained"),labels = c("missense variant","frameshift variant","disruptive inframe deletion","splice acceptor variant","splice donor variant","inframe deletion","inframe insertion","stop gained")))

draw(h1)
decorate_heatmap_body("mat", {
  i = which(colnames(mat) == "153PT18")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(fill="white", col = "white",lwd = 5))
  })

for (a in 1:nrow(r)){
  rw = r[a,2]
  cl = r[a,1]
  decorate_heatmap_body("mat", {
    grid.rect((rw-0.5)/ncol(mat), (cl-0.5)/nrow(mat),1/ncol(mat),1/nrow(mat), default.units = "npc",gp = gpar(fill = "grey89",col = "white"))
  })
}


#dev.off()




