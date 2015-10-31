library(ggplot2)
library(rjson)
library(reshape2)


# read data
fname <- "results/gene_codon_usages.json"
data <- fromJSON(file=fname, method='C')

# preprocess data
df <- data.frame(marker=numeric(0), usage=numeric(0))
for(cur in data) {
    mrkr <- cur$marker
    usg <- cur$usage

    df <- rbind(df, data.frame(marker=mrkr, usage=usg))
}

# plot data
ggplot(df, aes(x=usage, fill=..count..)) +
    geom_histogram(binwidth=0.005) +
    xlab("codon usage") +
    scale_y_log10() +
    facet_grid(marker ~ .) +
    ggsave(filename=paste0("codon_usage_hist.png"))
