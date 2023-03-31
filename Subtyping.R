library(Rtsne)
library(CancerSubtypes)
library(umap)
library(otrimle)
library(dplyr)
library(survival)
require(survRM2)
library("survminer")
require("survival")
set.seed(235)
#read the gene expression profiles of patients
expression <- read.table("LUNG_Gene_Expression.txt")
#read the clinical data of patients
Survdata = read.table("LUNG_Survival.txt",sep="\t", header = TRUE)


summarise(expression)
expression[1:3,1:3]
dim(expression)
dim(Survdata)

######## missing value imputation ###########
index=which(is.na(expression))
expression <- data.imputation(expression,fun="median")
expression = data.normalization(expression,type="feature_Mean",log2=FALSE)
dim(expression)
#transpose the data
lung_expr <- t(as.matrix(expression))
dim(lung_expr)

##########################################


# Initialize hyper-parameters

LUNGparameters<- tibble(min_dist = numeric(),
                                # n_neighbors = numeric(),
                                n_components = numeric(),
                                p_value = character(),
                                #gst = double(),
                                RLEDMIN = numeric())




##### Hyper- parameter tuning ######

for (min_dist in seq(0.001,1.0,0.001)){
  #for (n_neighbors in seq(5,99))
  #{
  for (n_components in seq(2,50,2))
  {
    set.seed(2523)
    umap.lung <- umap(lung_expr,  
                      min_dist = min_dist, n_components = n_components)
    
    
    
    
    
    #perform the GMM based clustering
    
    otri=otrimle(umap.lung$layout, 4, initial = NULL, logicd = NULL, npr.max =6/100, erc = 50,
                 iter.max = 100, tol = 1e-06, ncores = NULL, monitor = TRUE)
    
    
    
    ## record cluster labels
    C<-otri$cluster
    C
    
    clinical <- cbind(Survdata,C)
    colnames(clinical)<-c('PatientID', 'Survival', 'Death', 'C')
    
    ### remove noise cluster information from both expression and clinical data ###
    
    noise<-length(which(C==0))
    lung_expr_new <- lung_expr
    clinical_new <- clinical
    if(noise > 0 )
    {
      idx = which(clinical_new$C==0)
      clinical<-clinical_new[-idx,]
       lung_expr_new <- lung_expr_new[-idx,]
    }
    
    #check if any cluster has single observation, then simply continue through the loop
    
    C1 <- length(which(C==1))
    C2 <- length((which(C==2)))
    C3 <- length(which(C==3))
    C4 <- length((which(C==4)))
    
    if(C1 <2 | C2 <2 | C3 <2 |C4 <4 ){
      next
    }
    dim(lung_expr_new)
    dim(clinical_new)
    
    
    #fit the Survival data into the function
    clinical_new<-as.data.frame(clinical)
    fitQ <- survfit(Surv(Survival,Death) ~ as.factor(C), data = clinical_new)
    
    #calculate Pvalue
    p <- surv_pvalue(fitQ)
    print("pvalue is")
    print(p[2])
    pval <- c(p[2])


    cv<-otri$size
    
    #size of each cluster
    cv
 
    
    
    #quantify the separability of survival curves in terms of various metrics a  the details of the metrics can be found in the papers

    # this function is adapted from https://github.com/angy89
    survival.curves.separation <- function(data, cluster, tau){
      
      ans <- list()
      tmp           <- table(cluster)
      ClusterLabel  <- as.numeric(names(tmp))
      ClusterSize   <- as.numeric(tmp)
      K             <- length(ClusterLabel)
      n             <- sum(ClusterSize)
      
      ## Set a time grid
      time_grid  <- min(data$Survival) : min(max(data$Survival), tau)
      time_scale <- 1/diff(range(time_grid))
      
      ## Create estimated survival curve matrix: [time x clusters]
      H <- matrix(NA, nrow=length(time_grid), ncol=K)
      colnames(H) <- paste('Cluster_', ClusterLabel, sep='')
      
      ## Estimate the KM curve clusterwise on a common support
      for(k in 1:K){
        ## Compute Kaplan-Meier estimator on the kth cluster
        km    <- survfit(Surv(Survival, Death)~1, data = data[cluster==ClusterLabel[k] , ] )
        ## Construct the KM estimator function
        KMfun <- stepfun(x=km$time[-1], y=km$surv)
        H[,k] <- KMfun(time_grid)
      }
      
      ## construct matrix of pairwise L1 distances
      D <- matrix(0, ncol=K, nrow=K)
      for (i in 1:K){
        for(j in 1:K)
          if(i!=j){
            D[i,j] <- D[j,i] <- sum( abs( H[ , i]  -  H[ , j] ))
          }
      }
      ## Some scaling is given so that these numbers are somewhow interpretable
      ## for the same number of clusters independently of the time interval
      D <- D * time_scale
      
      
      ## Metric 1: min pairwise L1 distance
      iut <- which(upper.tri(D, diag=FALSE))
      ans$L1min <- min(D[iut])
      
      ## Metric 2: compute the summed L1 distance of each
      ## cluster to the nearest one
      diag(D)   <- NA
      ans$L1sum <- sum( D[iut] ) / length(iut)
      
      
      return(ans)
    }
    
    
    
    # this function is adapted from https://github.com/angy89
    
    rmst.separation <- function(data, cluster, tau) {
      ans <- list()
      tmp           <- table(cluster)
      ClusterLabel  <- as.numeric(names(tmp))
      ClusterSize   <- as.numeric(tmp)
      K             <- length(ClusterLabel)
      n             <- sum(ClusterSize)
      
      
      ## Compute the minimum of the largest observed time in each of the two groups
      max.time <- rep(0, K)
      for (k in 1:K) {
        max.time[k] <- max(data$Survival[cluster == ClusterLabel[k]])
      }
      TAU  <- min(max.time, tau)
      
      
      ## Names
      ##    * RMST = restricted mean survival time
      ##    * LER  = life expectancy ratio
      ##
      ## LER: Life Expectancy Ratio  Matrix
      ##      LER[i,j] = LER[j,i] = max{RMST[i] / RMST[j], RMST[j] / RMST[i]}
      ##      note that here we don't have a baseline group so we define the ratio
      ##      always using in the denominator the group that have smaller RMST
      ##
      ## LED: Life Expectancy Difference
      ##    LED[i,j] = LED[j,i] = abs(RMST[i] - RMST[j])
      ##    note that here we don't have a baseline group so we define tha abs difference
      ##
      LER <- LED <-  matrix(0, ncol = K, nrow = K)
      for (i in 1:K) {
        for (j in 1:K)
          if (i != j) {
            ## First select data from  the two groups
            idx <- { cluster == ClusterLabel[i] | cluster == ClusterLabel[j]  }
            x   <- data[idx,]
            ##  Create a 0-1 vector, with gr==1 if cluster==ClusterLabel[i]
            gr0  <- cluster[idx]
            gr   <- ifelse(gr0 == ClusterLabel[i], 1, 0)
            u    <- rmst2(time = x$Survival, status = x$Death, arm = gr, tau = TAU)
            
            rmst_i <- u$RMST.arm1$rmst[1]
            rmst_j <- u$RMST.arm0$rmst[1]
            
            LER[i,j]  <- LER[j, i] <- max(rmst_i / rmst_j, rmst_j / rmst_i  )
            LED[i, j] <- LED[j, i] <- abs(rmst_i - rmst_j)
          }
      }
      
      ## index of the upper triangle
      iut <- which(upper.tri(LER, diag = FALSE))
      
      ## metric: min of pairwise LER discrepancy
      ans$LERmin <- min(LER[iut])
      
      ## metric: scaled summed pairwise LER discrepancy
      ans$LERsum <- sum(LER[iut]) / length(iut)
      
      ## metric: min of pairwise LED discrepancy
      ans$LEDmin <- min(LED[iut])
      
      ## metric: scaled summed pairwise LED discrepancy
      ans$LEDsum <- sum(LED[iut]) / length(iut)
      
      
      return(ans)
    }
    
    # this function is adapted from https://github.com/angy89
    survPval = function(SurvDat,CLg,nYears=5){
      
      fit <- survfit(Surv(Survival, Death) ~ CLg,data = SurvDat,subset = Survival < (365 * nYears))
      suv <- survminer::ggsurvplot(fit, risk.table = TRUE, risk.table.height = 0.5,
                                   xlim = c(0,5000), break.time.by = 500, pval = TRUE)
      pVal <- survminer:::surv_pvalue(fit,method = "survdiff",data=SurvDat)
      return(list(fit=fit,suv = suv,pVal=pVal))
      
    }
    
    
    ## Computes the discrepancy between survival curves in temrs of RMST
    ## (restricted mean survival time).
    
    ##
    ## Outputs: a list with multiple separation measures
    ##
    ##
  
    
    ###
    
    
    #########################################################
    ########################################################
    ########################################################
  
    survRes = list()
    nYears = 11
    tau = 1825
    
    dim(clinical_new)
    clinical<-clinical[order(clinical$C),]
    l2d1 <- dim(clinical_new)[1]
    ClustGroup = clinical_new[1:l2d1,4:4]
    SurvDat = clinical_new[1:l2d1,2:3]
    CLg = ClustGroup
    
    
    LNormDist = c(unlist(survival.curves.separation(SurvDat, CLg,tau)),
                  unlist(rmst.separation(SurvDat, CLg,tau)))
    
   
    str(LNormDist)
    print('the value is)))))))))))))))))))))))))))))))))))))))))))))))))')
    print(LNormDist[5])
    print(p[2])
    
    ledmin <- as.integer(LNormDist[5])
    p <- as.double(p[2])
    
    LUNGparameters <- LUNGparameters %>% add_row(min_dist =min_dist,
                                                                 # n_neighbors = n_neighbors,
                                                                 n_components = n_components,
                                                                 p_value = (p = as.character(p)),
                                                                 #gst = x,
                                                                 RLEDMIN = ledmin)



    
    }}
    
    
    
    
    
    
    
    
    
write.table(LUNGparameters, file = "LUNGparameters.txt", sep = "\t",
                             row.names = FALSE)
    
View(LUNGparameters)
    
    
LUNG.hyperparameters <-  as_tibble(as.matrix(read.delim("LUNGparameters.txt")))


# keep only those rledmin values with significant  p value and high rledmin
#compared to benchmark papers

rem <- which(LUNG.hyperparameters$RLEDMIN <6 )
LUNG.hyperparameters <- LUNG.hyperparameters[-rem,]
rem <- which(LUNG.hyperparameters$p_value >0.06)
LUNG.hyperparameters <- LUNG.hyperparameters[-rem,]

#plot the hyperparameters and RLED_min

LUNG.plot <- ggplot(LUNG.hyperparameters, aes(x = n_components, y= RLEDMIN))
LUNG.plot  + geom_point(col = "cyan")+ geom_line( col = "red") + theme_light() + scale_x_continuous(breaks = seq(10, 60, by = 1)) +  labs(x = "n_components", y ='RledMin' )

  #3D interactive plot
library(plotly)
plot_ly(x=LUNG.hyperparameters$min_dist, y=LUNG.hyperparameters$n_components, z=LUNG.hyperparameters$RLEDMIN, type="scatter3d", mode="markers" )
    
#From LUNG.hyperparameters pull up the min_dist and n_components correspoding to which RLED_min is maximum
#rerun the UMAP and Otrimle but with the same seed and draw the Kaplan-Meier plots
  
  
  set.seed(2523)
  umap.lung <- umap(lung_expr,  
                    min_dist = _, n_components = _) # substitute _ with values that maximize RLED_min 
  
  
  otri=otrimle(umap.lung$layout, 4, initial = NULL, logicd = NULL, npr.max =6/100, erc = 50,
               iter.max = 100, tol = 1e-06, ncores = NULL, monitor = TRUE)
  
  
  C<-otri$cluster
  C
  
  clinical <- cbind(Survdata,C)
  colnames(clinical)<-c('PatientID', 'Survival', 'Death', 'C')
  
  
  noise<-length(which(C==0))
  lung_expr_new <- lung_expr
  clinical_new <- clinical
  if(noise > 0 )
  {
    idx = which(clinical_new$C==0)
    clinical<-clinical_new[-idx,]
    lung_expr_new <- lung_expr_new[-idx,]
  }
  
  
  #fit the survivla data into the function
  clinical_new<-as.data.frame(clinical)
  fitQ <- survfit(Surv(Survival,Death) ~ as.factor(C), data = clinical_new)
  
  #calculate Pvalue
  p <- surv_pvalue(fitQ)
  print("pvalue is")
  print(p[2])
  
  
  
  
  
  
  
  #drawing survival curves
  ggsurvplot(fitQ,
             pval = TRUE,
             pval.coord = c(0, 0.02),
             palette = c("red","blue","green", "black"),
             conf.int = FALSE,
             font.x = c(30, "plain", "black"),
             font.family="Roman",
             font.y = c(30, "plain", "black"),
             pval.size=12,
             censor.size=5,
             font.title=c(30,"Roman", "black"),
             font.tickslab=c(32, "plain","black"),
             font.legend=c(32, "plain","black"),
             #legend.title = "Fig.A",
             break.time.by=  1000,
             xlab = "Time (Days)",
             legend.title = "Subtype",

             legend.labs = c("C=1", "C=2", "C=3", "C=4" ),
             legend=c(0.8,0.75))


  
     
     
    #Cox Analysis
    fitcox <- coxph(Surv(as.numeric(Survival),as.numeric(Death)) ~ as.factor(C), data = clinical)
    
    summary(fitcox)
    ftest <- cox.zph(fitcox,terms=FALSE)
    ggcoxzph(ftest)
    
    dev.off()







require(org.Hs.eg.db)
library(limma)
require(clusterProfiler)


source("gene_mapping.R")
source("pathway_analysis.R")

C <- clinical$C
red_dat<-lung_expr
dim(red_dat)
red_dat <- cbind(red_dat, C)
dim(red_dat)
length(C)
C <- c(C)
pathway_matrix <- red_dat
#pathway_matrix <- pathway_matrix[-106,]


dim(pathway_matrix)
pathway_matrix[1:3,1:3]

pathway_df<-as.data.frame(pathway_matrix)
str(pathway_df)
pathway_df[1:3,1:3]
pathway_df <- pathway_df[order(pathway_df$C),]
#pathway_df[400:410,1:3]
length(pathway_df$C)
pathway_df[1:3,1:14]
clusterR<-pathway_df$C

clusterR
FR = pathway_analysis(cluster=clusterR,pathway_df)
dev.off()
FR$dp
