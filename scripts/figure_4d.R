
 
library(ggplot2)
data <- read.table("figures/figure_4d.txt", header=T)

data$InformativeLoci <- factor(data$InformativeLoci, levels = data$InformativeLoci)

ggplot(data, aes(x=InformativeLoci, y=Count)) + geom_bar(stat="identity") + labs(x="Number of Informative MIRU loci") +
   theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
ggsave("figures/figure_4d.png")
