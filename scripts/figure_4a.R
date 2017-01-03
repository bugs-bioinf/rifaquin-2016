
 
library(ggplot2)
data <- read.table("figures/figure_4a.txt", colClasses=("character"))
colnames(data)<-c("name","count","snps")

# remove mixed infection
data <- data[-which(data$name =='035'),]

data2 <- data.frame(
	SNPs = c(0,1,2,3,4,5,6,7,8,9,10,"","737",">1200"), 
	Counts = c(
		length(which(data$snps==0)),
		length(which(data$snps==1)),
		length(which(data$snps==2)),
		length(which(data$snps==3)),
		length(which(data$snps==4)),
		length(which(data$snps==5)),
		length(which(data$snps==6)),
		length(which(data$snps==7)),
		length(which(data$snps==8)),
		length(which(data$snps==9)),
		length(which(data$snps==10)),
		0,
		length(which(data$snps==737)),
		length(which(as.numeric(data$snps)>1200)))
	)

data2$SNPs <- factor(data2$SNPs, levels = data2$SNPs)
ggplot(data2, aes(x=SNPs, y=Counts)) + geom_bar(stat="identity") + labs(x="SNP Differences Between Isolates", y="Isolate Pairs Count") +
   theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
ggsave("figures/figure_4a.png")
