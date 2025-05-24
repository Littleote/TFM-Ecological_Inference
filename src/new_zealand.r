library(ei.Datasets)
NZ <- ei_NZ_2020
counts <- mapply(function(t) {nrow(t)}, as.list(NZ[,3])[[1]])
biggest <- rev(order(counts))
name <- as.list(NZ[biggest,2])[[1]]
print(cbind(name, counts))