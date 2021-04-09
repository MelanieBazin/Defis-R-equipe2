###
#  Melanie_Script clustering level 2
#  D'après scrpit de : Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
###

######
# Variable à choisir
######

#fichier de donnée
expMatrix = read.table("Mito_Genes.txt", header = T, row.names = 1) 
#nombre de cluster que l'on veux faire
nb_cluster = 4 
#algorithe à utiliser pour faire les clusters
method = "kmeans"               # "kmeans" -ou- "HCL" 
#méthode de calcule des distances
distance = "Euclidean"          #"Euclidean" -ou- "Correlation"
#type de graphique que l'on veux en sortie
graph = "heatmap" #c("heatmap","profils") -ou- "heatmap" -ou- "profils" 


##########

# This function is useful to draw gene expression profiles
plotGenes <- function(expData, title = "", yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(0, ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}


##########

# Choisir le mode de calcule des distances
if (distance == "Euclidean"){
  matDist = dist(expMatrix)
}else if (distance == "Correlation"){
  matDist = as.dist(1 - cor(t(expMatrix)))
}  

# Choisir le type d'algorithme utilisé pour faire les clusters
if (method  == "kmeans"){
  res = kmeans(matDist, nb_cluster)
  vecCluster = res$cluster
}else if(method  == "HCL"){
  res = hclust(matDist)
  vecCluster = cutree(res, nb_cluster)
}


pdf(paste0("Lvl2_profil_",nb_cluster,"clusters_",method,"_",distance,".pdf"))
  for (i in 1:nb_cluster){
    geneCluster = names(which(vecCluster == i))
    cluster = expMatrix[geneCluster,]
    
    if (is.element(TRUE,graph_type == "profils")) {
      plotGenes(cluster, 
                title = paste("Cluster", selected_cluster, "\n",
                              "Distance :",distance,"-",
                              "Algorithme :", method, "avec ", nb_cluster, "clusters")
                , yMax = max(expMatrix)) 
    }
    if (is.element(TRUE,graph_type == "heatmap")){
      heatmap(as.matrix(cluster),
              Colv = NA, Rowv = NA,
              main  = paste("\n","\n","Cluster", selected_cluster,"\n","Distance :",distance,"-",
                            "Algorithme :", method,  nb_cluster, "clusters"))
    }
  }
dev.off()

