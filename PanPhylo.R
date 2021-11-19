library(ape)
library(dplyr)
library(reshape2)
library(ggplot2)
library(bnlearn)
library(ade)

#################### CREATE DATA #################### 
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
to_rm = which(colMeans(genes_mat)==0)
genes_mat = genes_mat[,colSums(genes_mat) > 0]

phylo_dists = melt(cophenetic(tree))
#phylo_dists = phylo_dists[!duplicated(t(apply(phylo_dists[c("Var1", "Var2")], 1, sort))), ]
phylo_dists = phylo_dists[phylo_dists$Var1 != phylo_dists$Var2,]


#################### CREATE NETWORK ####################
library(phytools)
tips = tree$tip.label
tree_mat = as.data.frame(tree$edge)
tree_mat$V2 = vapply(tree_mat$V2, function(x) ifelse(x <= length(tree$tip.label), tree$tip.label[x], paste0('node_',as.character(x))), FUN.VALUE = character(1))
tree_mat$V1 = paste0('node_',as.character(tree_mat$V1))
tree_mat$Length = tree$edge.length

tree_depths = node.depth(tree)
tree_mat$joint_proba_missing = vapply(1:nrow(tree_mat), function(x) ifelse(tree_mat$V2[x] %in% names(completions), 1 - completions[names(completions) == tree_mat$V2[x]], -999), FUN.VALUE = numeric(1))
tree_mat$joint_proba_missing[tree_mat$joint_proba_missing == -999] = NA

names(tree_depths) = c(tree$tip.label,paste0('node_',as.character(10+(1:tree$Nnode))))
for (x in seq(2,max(tree_depths))){
  nodes = names(tree_depths)[tree_depths == x]
  for (node in nodes){
    tree_mat$joint_proba_missing[tree_mat$V2 == node] = prod(tree_mat$joint_proba_missing[tree_mat$V1 == node])}}
root = unique(tree_mat$V1[!(tree_mat$V1 %in% tree_mat$V2)])
tree_mat = rbind(tree_mat, data.frame(V1 = 'root', V2 = root, Length = 0, joint_proba_missing = prod(tree_mat$joint_proba_missing[tree_mat$V1 == root])))



tree_net = empty.graph(nodes = c(tree$tip.label, paste0('node_',as.character((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)))))
for (i in 1:nrow(tree_mat)){tree_net = set.arc(tree_net, from = tree_mat$V1[i], to = tree_mat$V2[i])}

full_df = data.frame()
for (gene in 1:ncol(genes_mat)){
gene_data = genes_mat[,gene]
gene_test = matrix(ncol = length(tree_net$nodes))

values_leaves = as.data.frame(t(gene_data))
values_nodes = list()
for (node in 1:tree$Nnode){
  node_n = as.character(length(tree$tip.label)+node)
  node_name = paste0('node_',node_n)
  desc_nodes = getDescendants(tree, node_n)
  desc_leaves = desc_nodes[desc_nodes <= length(tree$tip.label)]
  desc_leaves = vapply(desc_leaves, function(i) tree$tip.label[i], FUN.VALUE = character(1))
  prop_gene = mean(gene_data[names(gene_data) %in% desc_leaves])
  values_nodes[node_name] = prop_gene}
values = cbind(values_leaves, as.data.frame(t(values_nodes)))
full_df = rbind(full_df, values)}
full_df[] <- lapply(full_df, as.numeric)
dfit = bn.fit(tree_net, full_df[1,])
dnet = bn.net(dfit)


















colnames(tree_mat) = c('from','to')
tree_mat = tree_mat[!(tree_mat$to %in% 1:length(tree$tip.label)),]
nodes_dists = dist.nodes(tree)
nodes_i = (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
nodes_dists = nodes_dists[nodes_i,nodes_i]

n_tips = length(tree$tip.label)
nod_beg = n_tips+1
nod_end = length(tree$tip.label)+tree$Nnode
nodes_depths = node.depth(tree)[nod_beg:nod_end]


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


  