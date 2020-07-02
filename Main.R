#### Initialization ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(tidyverse)
library(magrittr)
library(GeneNet)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggforce)
library(openxlsx)
library(readxl)
library(tictoc)

# import helper functions
source("HelperFunctions.R")

#### Load Data ----

tic()

# download unprocessed glycomics data from figshare
file <- "GlycanData_figshare.xls"
load.web.file(
  url="https://ndownloader.figshare.com/files/23528105",
  md5sum = "adf8cde4bea1ad2ca6045b02cef3a2c5",
  outfile = file
)

# define data file name
file <- "GlycanData_figshare.xls"
# extract sheet names
sheets <- excel_sheets(path = file)
# define cohort names (exclude legend)
cohort_names <- sheets[sheets != "Legend"]

# define number of bootstraps
# Note: depending on this number, the ranking of the normalizations in the barplot below might change slightly, 
# but the results will be qualitatively equivalent
nboot <- 100 # warning: will take substantial amount of time for nboot=1000

# loop over all cohorts
results <- lapply(cohort_names %>% {names(.)=.;.}, function(cohort){
  printf("%s...\n", cohort)
  platform <- "LC-ESI-MS"
  if(cohort=="CRC") platform <- "UPLC"
  if(cohort=="LLS") platform <- "MALDI"
  
  # load data
  data_raw <- read_excel(file, sheet = cohort, col_names = T) %>% as.data.frame()
  # convert data to numeric 
  data_raw <- apply(data_raw,2,as.numeric) %>% as.data.frame()
  # impute zeros with NAs
  data_raw[data_raw == 0] <- NA
  # remove NAs
  data_raw <- data_raw %>% drop_na()
  
  # get age
  age <- data_raw$Age 
  # get glycan data
  data <- data_raw[2:dim(data_raw)[2]]
  
  # get adjacency matrix
  adja <- get_adja(platform = platform)
  
  ##### Normalize Data ----
  
  # for LC-ESI-MS data all normalizations are also performed per IgG subclass
  dt <- list()
  
  # Quotient Normalization
  dt$quot <- quotNorm(data)$X
  if (platform %in% c("LC-ESI-MS")) {
    dt$quotSub <- quotNorm(data)$X
    dt$quotSub[,1:20] <- quotNorm(data[,1:20])$X
    dt$quotSub[,21:40] <- quotNorm(data[,21:40])$X
    dt$quotSub[,41:50] <- quotNorm(data[,41:50])$X
    
    dt$quotSubLog <- log(dt$quotSub)
  }
  dt$quotLog <- log(dt$quot)
  
  # Total Area
  dt$TA <- tanorm(data)
  if (platform %in% c("LC-ESI-MS")) {
    dt$TASub <- tanorm(data)
    dt$TASub[,1:20] <- tanorm(data[,1:20])
    dt$TASub[,21:40] <- tanorm(data[,21:40])
    dt$TASub[,41:50] <- tanorm(data[,41:50])
    
    dt$TASubLog <- log(dt$TASub)
  }
  dt$TALog <- log(dt$TA)
  
  # Quantile Normalization
  dt$quant <- quantile_normalisation(data, 2) #column wise
  if (platform %in% c("LC-ESI-MS")) {
    dt$quantSub <- quantile_normalisation(data, 2)
    dt$quantSub[,1:20] <- quantile_normalisation(data[,1:20],2)
    dt$quantSub[,21:40] <- quantile_normalisation(data[,21:40],2)
    dt$quantSub[,41:50] <- quantile_normalisation(data[,41:50],2)
    
    dt$quantSubLog <- log(dt$quantSub)
  }
  dt$quantLog <- log(dt$quant)
  
  # Raw data
  dt$raw <- data
  dt$rawLog <- log(dt$raw)
  
  # Median Normalization
  dt$med <- mednorm(data)
  
  # Rank normalization
  dt$rank <- ranknorm(data, dir = 2)
  dt$rankLog <- log(dt$rank)
  
  # TAQuotient Normalization
  dt$TAquot <- quotNorm(tanorm(data))$X
  if (platform %in% c("LC-ESI-MS")) {
    dt$TAquotSub <- quotNorm(tanorm(data))$X
    dt$TAquotSub[,1:20] <- quotNorm(tanorm(data[,1:20]))$X
    dt$TAquotSub[,21:40] <- quotNorm(tanorm(data[,21:40]))$X
    dt$TAquotSub[,41:50] <- quotNorm(tanorm(data[,41:50]))$X
    
    dt$TAquotSubLog <- log(dt$TAquotSub)
  }
  dt$TAquotLog <- log(dt$TAquot)
  
  #### Compute GGM overlap to adjacency ----
  
  # define significance threshold for GGM p-value cutoff
  alpha <- 0.01
  
  # initialize result list
  fis_p <- list()
  
  # loop over methods
  norm_fis <- sapply(names(dt) %>% {names(.)=.;.}, function(x){
    printf("  %s\n", x)
    # initialize vector
    fis_p[[x]] <- array(data = NA, dim = (nboot+1))
    
    # bootstrap data
    for (i in 0:nboot) {
      if (i > 0) {
        data_re <- dt[[x]][sample(nrow(dt[[x]]), nrow(dt[[x]]), replace = TRUE), ]
      } else{
        data_re <- dt[[x]]
      }
      
      # compute partial correlations
      corr <- ggm.estimate.pcor(as.matrix(data_re), method = "dynamic", verbose=F)
      # compute pvalues
      gn_pvalues <- network.test.edges.silent(corr, verbose=F, plot=F)
      # adjust pvalue
      gn_pvalues$p.adj <- p.adjust(gn_pvalues$pval, method = "fdr")
      # find correlation cutoff
      gn_pvalues$pval[gn_pvalues$p.adj > alpha] <- NaN
      cutoff <- abs(gn_pvalues$pcor[which.max(gn_pvalues$pval)])
      
      # create GGM adjacency
      adja_data <- corr
      adja_data[abs(adja_data) >= cutoff] <- 1
      adja_data[abs(adja_data) < cutoff] <- 0
      adja_data[lower.tri(adja_data, diag = TRUE)] <- NA
      
      # compare GGM with prior knowledge
      contin <- contab(adja, adja_data)
      
      # compute Fisher's test pvalue
      fis_p[[x]][i+1] <- fisher.test(contin)$p.value 
    }
    # compute fisher's test pvalues confidence intervals (log-scale)
    confin(fis_p[[x]][2:nboot+1] %>% log10)
  })
  
  #### Create barplot ----
  
  # arrange data for plotting
  xy <- norm_fis %>% as.data.frame %>% t
  colnames(xy) <- c("fpvalmin","fpval","fpvalmax")
  xy %<>% as.data.frame
  xy$name <- rownames(xy)
  
  # define grouping variable to make colors same as in paper
  xy$group <- 1
  xy$group[which(xy$name %in% c("quant","quantLog","quantSub","quantSubLog"))] <- 2
  xy$group[which(xy$name %in% c("quot","quotLog","quotSub","quotSubLog"))] <- 3
  xy$group[which(xy$name %in% c("rank","rankLog"))] <- 4
  xy$group[which(xy$name %in% c("med"))] <- 5
  xy$group[which(xy$name %in% c("raw","rawLog"))] <- 6
  xy$group[which(xy$name %in% c("TAquot","TAquotLog","TAquotSub","TAquotSubLog"))] <- 7
  
  # generate plot
  p <- ggplot(xy, aes(x=reorder(name, fpval),y=-fpval,fill=as.factor(group)))+
    geom_bar(stat="identity",colour="black")+
    coord_flip() +
    geom_errorbar(aes(ymin=-fpvalmax,ymax=-fpvalmin), width=.5)+
    theme(axis.text.x=element_text(angle=90,hjust=1,size=14)) +
    theme(axis.title.x=element_text(angle=0,hjust=1,size=14)) +
    theme(axis.text.y=element_text(angle=0,hjust=1,size=14)) +
    theme(legend.text=element_text(size=14)) + 
    theme(plot.title = element_text(size=14)) +
    xlab("") +
    ylab("-log10(Fisher's p-value)")+
    ggtitle("Normalizations")+
    scale_fill_manual(values = brewer.pal(n = 7, name = "Set3"),
                      name="Normalizations", 
                      labels = c("Total Area", 
                                 "Quantile",
                                 "Quotient",
                                 "Rank",
                                 "Median",
                                 "Raw",
                                 "Total Area Quotient"))
  
  # save plot to file  
  pdf(sprintf("Barplot_%s.pdf",cohort), width = 10.5, height= 10)
  print(p)
  dev.off()
  
  #### Check association with Age ----
  
  # loop over methods
  xx <- sapply(dt %>% names %>% {names(.)=.;.}, function(x){
    # loop over glycans
    sapply(1:dim(dt[[x]])[2], function(k) {
      fit <- lm(dt[[x]][,k]~age)
      # get p-value of the association with age
      coef(summary(fit))["age","Pr(>|t|)"]
    })
  }) %>% as.data.frame()
  rownames(xx) <- colnames(data)
  # adjust pvalues
  zz <- apply(xx, 2, function(x) {p.adjust(x, method="fdr")})
  
  # melt dataframe
  yy <- zz %>% log10 %>% as.matrix %>% 
    reshape2::melt()
  
  # define grouping variable to make colors same as in paper
  yy$group <- 1
  yy$group[which(yy$Var2 %in% c("quant","quantLog","quantSub","quantSubLog"))] <- 2
  yy$group[which(yy$Var2 %in% c("quot","quotLog","quotSub","quotSubLog"))] <- 3
  yy$group[which(yy$Var2 %in% c("rank","rankLog"))] <- 4
  yy$group[which(yy$Var2 %in% c("med"))] <- 5
  yy$group[which(yy$Var2 %in% c("raw","rawLog"))] <- 6
  yy$group[which(yy$Var2 %in% c("TAquot","TAquotLog","TAquotSub","TAquotSubLog"))] <- 7
  
  # define grouping variable for Log normalizations
  yy$shape <- 21
  yy$shape[which(yy$Var2 %in% c("TALog","quantLog","quotLog","rankLog","rawLog","TAquotLog","TASubLog","quantSubLog","quotSubLog","TAquotSubLog"))] <- 22
  
  nrow=5
  ncol=2
  pg <- ceiling(length(unique(yy$Var1))/(nrow*ncol))
  # for "LLS" cohort somehow facet_wrap_paginate crashes on the last page
  # therefore printing last page separate
  if(cohort=="LLS") pg <- pg-1

  # save plot to file
  pdf(sprintf("AgeAssociation_%s.pdf", cohort), width = 10.5, height= 18)
  for (i in 1:pg) {
    print(
      yy %>%
        ggplot(aes(x=Var2, y=-value)) +
        geom_point(aes(shape=as.factor(shape)), color="black", size=4) +
        geom_point(aes(color= as.factor(group), shape=as.factor(shape)), size=3) +
        geom_hline(yintercept = -log10(alpha), color="black", size=0.5, linetype="dashed") +
        theme_bw() +
        # theme(axis.text.x = element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_continuous(limits = c(0, NA)) +
        xlab("Normalization methods") +
        ylab("-log10(p.adj)") +
        scale_color_manual(values = brewer.pal(n = 7, name = "Set3"),
                           name="Normalizations",
                           labels = c("Total Area",
                                      "Quantile",
                                      "Quotient",
                                      "Rank",
                                      "Median",
                                      "Raw",
                                      "Total Area Quotient")) +
        scale_shape_discrete(name  ="Log",
                             labels=c("Non-Log", "Log")) +
        ggtitle(sprintf("%s: Glycan ~ Age", cohort)) +
        facet_wrap_paginate(~Var1,scales='free', ncol = ncol, nrow=nrow, page=i)
    )
  }
  if(cohort=="LLS") {
    p <- yy %>% 
      # selecting only glycans not included in the first pg-1 pages
      dplyr::filter(Var1 %in% levels(yy$Var1)[(pg*ncol*nrow+1):length(levels(yy$Var1))]) %>%
      # plot as before
      ggplot(aes(x=Var2, y=-value)) +
      geom_point(aes(shape=as.factor(shape)), color="black", size=4) +
      geom_point(aes(color= as.factor(group), shape=as.factor(shape)), size=3) +
      geom_hline(yintercept = -log10(alpha), color="black", size=0.5, linetype="dashed") +
      theme_bw() +
      # theme(axis.text.x = element_blank()) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(limits = c(0, NA)) +
      xlab("Normalization methods") +
      ylab("-log10(p.adj)") +
      scale_color_manual(values = brewer.pal(n = 7, name = "Set3"),
                         name="Normalizations",
                         labels = c("Total Area",
                                    "Quantile",
                                    "Quotient",
                                    "Rank",
                                    "Median",
                                    "Raw",
                                    "Total Area Quotient")) +
      scale_shape_discrete(name  ="Log",
                           labels=c("Non-Log", "Log")) +
      ggtitle(sprintf("%s: Glycan ~ Age", cohort))+
      facet_wrap(~Var1)
    print(
      ggarrange(p,
                ncol = ncol, nrow = nrow,
                common.legend = T)
    )
  }
  dev.off()
  
  list(fis=norm_fis, res=table(yy[yy$value<log10(alpha),"Var2"]), lm.pval=xx, lm.padj=zz)
})

val <- revert_list_str_4(results)

#### Write results to file ----

# get all result table
wb = createWorkbook()
lapply(cohort_names, function(x){
  sheet = addWorksheet(wb, sprintf("%s_pval",x))
  writeData(wb, sheet=sheet, val$lm.pval[[x]], rowNames=TRUE)
  sheet = addWorksheet(wb, sprintf("%s_padj",x))
  writeData(wb, sheet=sheet, val$lm.padj[[x]], rowNames=TRUE)
}) %>% invisible()
# remove file, if it exists
unlink("AgeAssociation_Results.xlsx")
# write out
saveWorkbook(wb, "AgeAssociation_Results.xlsx")

#### Plot summary table ----

res <- do.call(cbind, val$res[1:4]) %>% as.data.frame %>%
  dplyr::mutate(norm=rownames(.)) %>%
  left_join(val$res[["CRC"]] %>% as.data.frame %>% dplyr::rename(CRC=Freq), by=c("norm"="Var1"))
if("LLS" %in% cohort_names){
  res <- left_join(res, val$res[["LLS"]] %>% as.data.frame %>% dplyr::rename(LLS=Freq), by=c("norm"="Var1"))
}
rownames(res) <- res$norm
res$norm <- NULL

d <- c(50, 50, 50, 50,24,61)
f <- sapply(1:length(cohort_names),function(x){
  res[,x]/d[x]
}) %>% as.data.frame
colnames(f) <- cohort_names
rownames(f) <- rownames(res)
# computing average of LC-ESI-MS datasets
f$average_LCESIMS <- rowSums(f[1:4])/4
# compute weighted average across platforms
f$weighted_average <- rowSums(f[,c(5:(ncol(f)-1),ncol(f))], na.rm = T)/length(c(5:(ncol(f)-1), ncol(f)))

# remove normalizations per subclass (only valid for LC-ESI-MS)
f <-  f[-grep(pattern = "Sub", rownames(f)),] %>%
  dplyr::arrange(-weighted_average)
# reorder columns
f <- f[,c(2,1,3,4,(ncol(f)-1),5:(ncol(f)-2), ncol(f))]

# save table as image
ggtexttable(format(f, digits = 3 ,scientific = FALSE), 
            theme = ttheme(
              colnames.style = colnames_style(color = "white", fill = "lightsteelblue2"),
              rownames.style = colnames_style(color = "white", fill = "lightsteelblue2"),
              tbody.style = tbody_style(color = "black")
            )) +
  theme_void() +
  (if("LLS" %in% colnames(f)) {labs(title = "Model: Glycan ~ Age",subtitle = "Fraction of significant associations (FDR 0.01)")}else{labs(title = "Model: Glycan ~ Age (NOTE: data does not include LLS cohort)",subtitle = "Fraction of significant associations (FDR 0.01)")})
ggsave("Validation_Summary.pdf", width = 21, height = 13, units = "cm")

# save summary table to fille
wb = createWorkbook()
sheet = addWorksheet(wb, "Table")
writeData(wb, sheet=sheet, x= f, rowNames=TRUE)
if(!("LLS" %in% colnames(f))) {
  writeData(wb, sheet=sheet,x="(NOTE: data does not include LLS cohort)",startCol=ncol(f)+2, startRow=1)
}
# remove file, if it exists
unlink("Validation_Summary.xlsx")
# write out
saveWorkbook(wb, "Validation_Summary.xlsx")

toc()