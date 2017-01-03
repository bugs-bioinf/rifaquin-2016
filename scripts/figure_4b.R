
 
library(ggplot2)
data <- read.table("figures/figure_4b.txt", header=T)

ggplot(data, aes(x=MIRUchanges, y=Count)) + geom_bar(stat="identity") + labs(x="MIRU loci differences") +
   theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
ggsave("figures/figure_4b.png")
