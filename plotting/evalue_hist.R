library(ggplot2)


e.vals <- as.numeric(scan("e_values.dat", what="", sep="\n"))

pair_vec <- list()
pair_vec[["0"]] <- list(threshold=1, binwidth=0.02)
pair_vec[["1"]] <- list(threshold=0.25, binwidth=0.002)
pair_vec[["2"]] <- list(threshold=0.1, binwidth=0.001)

for(i in names(pair_vec)) {
    cur = pair_vec[[i]]

    tmp <- e.vals[e.vals<=cur$threshold]
    df <- data.frame(e=tmp)

    if(length(df$e) == 0) {
        next
    }

    ggplot(df, aes(x=e, fill=..count..)) +
        geom_histogram(binwidth=cur$binwidth) +
        xlab("e value") +
        ggsave(filename=paste0("evalue_hist", i, ".png"))
}
