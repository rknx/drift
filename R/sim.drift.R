#' Simulate Genetic Drift
#' 
#' This function allows you to simulate the change in allelic frequencies due to random drift in a Wright-Fisher model
#' @param n.pop Number of effective population
#' @param n Ploidy: Haploid=1, Diploid=2...
#' @param freq initial frequency of allele (0-1)
#' @param n.gen Number of generations to simulate Default:100
#' @param n.gene Number of genes/loci tested
#' @param n.time Number of simulations to perform
#' @param s.plot Save the image of plot? (TRUE/FALSE) Default:FALSE
#' @param s,table Save the table of numerical valuest? (TRUE/FALSE) Default:FALSE
#' @return Frequency Plots. If save is true, image and tabel text file are saved in working dierctory.
#' @usage sim.drift(n.pop,n,freq,n.gen,n.gene,n.time,s.plot,s.table)
#' @name sim.drift
#' @export
#' @examples
#' sim.drift(30,2,0.25,100,10,3,TRUE,FALSE)
#' sim.drift(30,2,0.25,100,10,3)

#Libraries
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2")}
library(ggplot2)                        # For plotting graph
library(reshape2)                        # For refrshing the plot for new simulations

# Monte-carlo function
sim.drift <- function (n.pop, n=1, freq, n.gen=100, n.gene, n.time, s.plot=FALSE, s.table=FALSE) {
  n.pop <- n.pop * n									# Adjust allele count for diploid
  timestamp <- floor(as.numeric(Sys.time()))			# Unique timestmamp for the files
  
  #Simulation
  for (k in 1:n.time) {									# For each simulation
    X = array(0, dim=c(n.gen,n.gene))					# Make array for data storage
    X[1,] = rep(n.pop*freq,n.gene)						# Start with initial allele count
    for(j in 1:n.gene){									# For each simulation
      for(i in 2:n.gen){								# For each generation
        X[i,j] = rbinom(1,n.pop,prob=X[i-1,j]/n.pop)	# Algorithm for frequency change by random drift
      }  
    }
    
    # Change table format and plot it
    X <- data.frame(X/n.pop)							# Normalize the dataframe
    colnames(X) = paste0("Line",1:n.gene)
    print(ggplot(melt(X), aes(x = rep(c(1:n.gen), n.gene), y = value, colour = variable)) + geom_line() + labs(title = "Simulations of Genetic Drift", x="Generation", y="Allele Frequency") + ylim(0,1) + theme(legend.position = "none"))
														# Plotting
    # Saving the files
    if (s.table == TRUE) {write.table(round(X,digits=3), file=paste0("Genetic Drift Simulation ",timestamp,"-",k,".txt"), sep="\t")}   # Optional, save CSV
    if (s.plot == TRUE) {ggsave(paste0("Genetic Drift Simulation ",timestamp,"-",k,".png"), type="cairo-png")}        # Optional, save graph
  }
  print(paste0("The timestamp is: ",timestamp))
}