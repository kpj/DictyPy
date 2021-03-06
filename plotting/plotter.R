library(ggplot2)
library(rjson)
library(reshape2)


fname <- "results/grouped_genes.json"
data <- fromJSON(file=fname, method='C')

df <- data.frame(group=character(0), AAA=numeric(0), AAG=numeric(0))
for(cur in data) {
  aaa <- cur$cumulative_codon_usage$AAA
  aag <- cur$cumulative_codon_usage$AAG
  sd <- round(sd(aaa), 3)

  df <- rbind(df, data.frame(group=paste0(cur$group, " ", "(sd: ", sd, ")"), AAA=mean(aaa), AAG=mean(aag)))
}
df <- df[order(df$AAA),]
df$group <- factor(df$group, levels=df$group, ordered=TRUE)

melted.df <- melt(df, id.var="group")
ggplot(melted.df, aes(x=group, y=value, fill=variable)) +
  geom_bar(stat="identity", width=.5) +
  geom_text(aes(label=ifelse(variable == "AAA", round(value, 2), ""), vjust=-0.25), size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4)) +
  ggsave(filename="codon_usage.png")
