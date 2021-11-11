library(ape)
library(dplyr)
library(reshape2)
library(ggplot2)

#tree = rtree(10)
#saveRDS(tree, file = "tree_ex.rds")
tree = readRDS(file = "tree_ex.rds")
plot(tree)

r = 10
c = 100

genes_mat = cbind(round(matrix(runif(r*c, max = 0.6), r, c)), round(matrix(rep(1,r*c), r, c)), round(rbind(matrix(rep(1,r*c/2),r/2,c),matrix(rep(0,r*c/2),r/2,c))))

# remove occurrences based on completion
low_comp <- function(vec, comp){
  for (i in 1:length(vec)){
    prob_comp = rbinom(1,1,comp)
    if (prob_comp == 0){vec[i] = 0}}
  return(vec)}

# apply the function
rownames(genes_mat) = c('t7','t2','t8','t10','t1','t4','t9','t6','t3','t5')
completions = runif(10, min = 0.6)
names(completions) = rownames(genes_mat)
for (i in 1:r){genes_mat[i,] = low_comp(genes_mat[i,], completions[i])}

phylo_dists = melt(cophenetic(tree))
#phylo_dists = phylo_dists[!duplicated(t(apply(phylo_dists[c("Var1", "Var2")], 1, sort))), ]
phylo_dists = phylo_dists[phylo_dists$Var1 != phylo_dists$Var2,]

out_test = data.frame()
for (gene in 1:ncol(genes_mat)){
  data = data.frame(occurrence=genes_mat[,gene],completion=completions,genome=rownames(genes_mat))
  occ = sum(data$occurrence)
  pres_prob = prod(data$completion[data$occurrence == 1])
  abs_prob = prod(data$completion[data$occurrence == 0])
  mpd = mean(phylo_dists$value[(phylo_dists$Var1 %in% data$genome[data$occurrence == 1]) & (phylo_dists$Var2 %in% data$genome[data$occurrence == 1])])
  out_test = rbind(out_test, data.frame(gene=gene,occ=occ,pres_prob=pres_prob,abs_prob=abs_prob,mpd=mpd))}

out_test$category = c(rep('random',c), rep('core',c), rep('phylo',c))
ggplot(out_test, aes(x=occ,y=mpd,color=category,size=abs_prob)) + geom_point(alpha=0.5)


comp_data = expand.grid(g1 = names(completions), g2 = names(completions))
comp_data = comp_data[comp_data$g1 != comp_data$g2,]
comp_data$joint_proba = vapply(1:nrow(comp_data), function(i) completions[names(completions) == comp_data$g1[i]] * completions[names(completions) == comp_data$g2[i]], FUN.VALUE=numeric(1))
comp_data$PD = vapply(1:nrow(comp_data), function(i)  phylo_dists$value[(phylo_dists$Var1 == comp_data$g1[i]) & (phylo_dists$Var2 == comp_data$g2[i])], FUN.VALUE = numeric(1))

for (gene in 1:ncol(genes_mat)){
  for (i in 1:nrow(comp_data)){
    s1 = as.character(comp_data$g1[i])
    s2 = as.character(comp_data$g2[i])
    comp_data[i,paste0('gene_',as.character(gene))] = floor((genes_mat[s1,gene] + genes_mat[s2,gene])/2)
  }}
to_rm = which(colMeans(comp_data[5:ncol(comp_data)])==0) + 4
comp_data = comp_data[-to_rm]

  