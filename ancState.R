library(ape)
library(phytools)
#read in the tree
mytree = read.tree("out_treemix0.treeout")
#check the tree
mytree
#get tree objects
str(mytree)
class(mytree) #phylo
#find the edge lengths or branch length of the tree
mytree$edge.length
#assign each object of the tree to some vector
branches <- mytree$edge
species <- mytree$tip.label
populations <- mytree$tip.label       
brlength <- mytree$edge.length
nodes <- mytree$Nnode #27


#designate character state (Medicago(1) or native(0) for each tip label (population)). Note the outgroup are given 0 for native plant 
char.tree <- c(0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1)
length(char.tree) #28 for 28 populations
names(char.tree) <- species
#since the branch lengths are lesser than 0, these lengths need to be modified to run the ancestral state reconstruction as I was getting an error. To do this, the recommended is using the compute.brlen() funtion in library(ape)
##now this is the new tree with modified branch lengths
newtree <- compute.brlen(mytree)

######using Ace ##############
ERreconstruction <- ace(char.tree, newtree, type="discrete", model="ER")
SYMreconstruction <- ace(char.tree, newtree, type="discrete", model="SYM")
ARDreconstruction <- ace(char.tree, newtree, type="discrete", model="ARD")

##transition rates
ERreconstruction$rates #2.715391
SYMreconstruction$rates #2.715391
ARDreconstruction$rates #2.649871 2.425755

#chi square to decide the model
1-pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 1) #0.900594 chose ARD model

######### using make.simmap ###################
#just trying things
mytree.anc <- make.simmap(newtree, x = char.tree, model="ARD") #still chosing this model
mytree.er  <- make.simmap(newtree, x = char.tree, model="ER")
mytree.sym  <- make.simmap(newtree, x = char.tree, model="SYM")
#the real deal
mytree.anc <- make.simmap(newtree, x = char.tree, model="ARD",nsim=2000,Q="mcmc",burnin=10000,samplefreq=50,pi=c(1.0,0.0))
mytree.mapsum <- describe.simmap(mytree.anc, plot=FALSE)
mytree.mapsum$count  
mytree.mapsum$count[,2]
mean(mytree.mapsum$count[,2])
quantile(mytree.mapsum$count[,2],probs=c(0.5,0.025,0.95))  
#50%  2.5% 95% 
#    4     1    11 
quantile(mytree.mapsum$count[,2],probs=c(0.5,0.05,0.95)) 
# 50%  5% 95%  
# 4   2  10 
quantile(mytree.mapsum$count[,3],probs=c(0.5,0.05,0.95))
# 50%  5% 95%  
# 7   4  11 
quantile(mytree.mapsum$count[,3],probs=c(0.5,0.025,0.95))
# 50% 2.5%  95% 
#   7    3   12 
##getting posterior probabilities for transitions from native to alfalfa
mean(mytree.mapsum$count[,2] >= 2) #0.9545

plot(mytree.mapsum)


##plot the tree
par(mar = c(3, 3, 3, 3))
plotTree(mytree.anc[[1]], ftype = "i")
nodelabels(node = as.numeric(rownames(mytree.mapsum$ace)), pie = mytree.mapsum$ace, piecol = cols, cex = 0.5)
tiplabels(pie=to.matrix(char.tree, sort(unique(char.tree))), piecol=cols, cex=0.2)
add.simmap.legend(colors = cols, x = 0, y = 24)
dev.copy(pdf, "ancState_tree.pdf")
dev.off()




