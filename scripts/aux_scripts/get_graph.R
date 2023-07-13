#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)

option_list <- list(
        make_option(c("-d", "--data_file"), type="character",
                help="indicating the path to the input file"),
        make_option(c("-x", "--x_name"), type="character", default=NULL,
                help="Set the plot height"),
        make_option(c("-y", "--y_name"), type="character", default=NULL,
                help="Set the plot height"),
        make_option(c("-g", "--group_name"), type="character", default=NULL,
                help="Set the plot height"),
        make_option(c("-w", "--wrap_name"), type="character", default=NULL,
                help="Set the plot height"),
        make_option(c("-o", "--output"), type="character", default="output",
                help="Output figure path"),
        make_option(c("-O", "--output_file"), type="character", default="matrices_correlation",
                help="Output figure names (without '.png' extension"),
        make_option(c("-W", "--width"), type="integer", default=11,
                help="Set the plot width"),
        make_option(c("-H", "--height"), type="integer", default=11,
                help="Set the plot height"),
        make_option(c("-n", "--names"), type="character", default = NULL,
                help="comma-sepparated string, indicating the names to the input matrices at the same order as input matrices.")
        )

opt <- parse_args(OptionParser(option_list=option_list))

data_table <- read.table(opt$data_file, header = T)
ggp <- ggplot(data_table, aes_string(x=opt$x_name, y=opt$y_name, group= opt$group_name)) + 
  geom_line(aes_string(col=opt$group_name),size=0.5) + theme_minimal()  + theme(text = element_text(size = 6),legend.key.size = unit(0.3, 'cm')) + facet_wrap(as.formula(paste("~",opt$wrap_name,""))) + scale_color_brewer(palette = "Dark2")

output_path_name <- file.path(opt$output, paste0(opt$output_file, ".png"))
png(output_path_name, width = opt$width, height = opt$height, units = "cm", res = 300)
print(ggp)     
dev.off() 