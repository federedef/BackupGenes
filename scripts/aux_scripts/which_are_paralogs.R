#!/usr/bin/env Rscript
library(biomaRt)
library(optparse)

#################### FUCNTIONS ####################
###################################################

get_paralogs <- function(gene_hgnc,identity_filter=30){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  results <- getBM(attributes = c("hsapiens_paralog_ensembl_gene",
                                "hsapiens_paralog_perc_id",
                                "hsapiens_paralog_perc_id_r1"),
                 filters = "hgnc_id",
                 values = gene_hgnc,
                 mart = human)
  # Stay with the paralogous with at least identity_filter %
  #print(results)
  ensemble_paralogs_names <- results[which(results$hsapiens_paralog_perc_id >= identity_filter | results$hsapiens_paralog_perc_id_r1 >= identity_filter),1]
  #print(ensemble_paralogs_names)
  ensemble_paralogs_names <- ensemble_paralogs_names[!is.na(ensemble_paralogs_names)]
  ensemble_paralogs_names <- ensemble_paralogs_names[ensemble_paralogs_names != "" ]
  #print(ensemble_paralogs_names)
  # Obtain the HGNCs from those paralogs.
  if (length(ensemble_paralogs_names)>0){
  	results2 <- getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name",
                                  "hgnc_id"
                                  ),
                   filters = "ensembl_gene_id",
                   values = ensemble_paralogs_names,
                   mart = human)
    hgnc_paralogs_names <- results2$hgnc_id
  	} else {
        hgnc_paralogs_names <- c()
  	}
  #print(hgnc_paralogs_names)
  return(hgnc_paralogs_names)
}

################## OPTPARSE #####################
#################################################

# which_are_paralogues.R -i /backupgens/backup_gens -o / -O backup_feautes 

option_list <- list(
  make_option(c("-i", "--input_table"), type="character",
              help="insert name of the table of pairs genes"),
  make_option(c("-o", "--output_path"), default=".", type="character",
              help="Output path"),
  make_option(c("-O", "--output_name"), type="character",
              help="Output name ")
)

opt <- parse_args(OptionParser(option_list=option_list))

gene_pairs_data <- read.table(opt$input_table, sep="\t", header=F, check.names=F)
gene_pairs_data <- cbind(gene_pairs_data, rep("Not_Paralogs",nrow(gene_pairs_data))) 
query_genes <- unique(gene_pairs_data[,1])

for (query_gene in query_genes) {
	paralogs <- get_paralogs(gene_hgnc=query_gene,identity_filter=1)
	#print(gene_pairs_data[which(gene_pairs_data[,1] == query_gene & gene_pairs_data[,2] %in% paralogs),])
	gene_pairs_data[which(gene_pairs_data[,1] == query_gene & gene_pairs_data[,2] %in% paralogs),3] <- "Paralogs"
}

full_output <- file.path(opt$output_path, opt$output_name)
write.table(gene_pairs_data,full_output, sep="\t", col.names=F, row.names=F, quote=F)