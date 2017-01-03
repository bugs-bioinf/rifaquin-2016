
library(ggplot2)
data <- read.table("figures/figure_4c.txt", header=T)

ggplot(data, aes(x=SNPs, y=MIRUdifferences)) + geom_point() + labs(x="SNP number", y="MIRU differences") +
   theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) + 
   geom_smooth(method=lm, se=FALSE) 
ggsave("figures/figure_4c.png")

