rm(list=ls())

library("ggplot2")
library("RColorBrewer")

wtcc <- read.table("MS.pvals.out", header=T, as.is=T)

wtcc$Start <- wtcc$Start/1000
for(j in 2:22){
  wtcc[wtcc$Chr==j,"Start"] <- max(wtcc[wtcc$Chr==(j-1),"Start"])+wtcc[wtcc$Chr==j,"Start"]
  wtcc[wtcc$Chr==j, "Tick"] <- (min(wtcc[wtcc$Chr==j,"Start"]) + max(wtcc[wtcc$Chr==j,"Start"]))/2
}
wtcc[wtcc$Chr==1,"Tick"] <- (min(wtcc[wtcc$Chr==1,"Start"]) + max(wtcc[wtcc$Chr==1,"Start"]))/2
wtcc$Discovery_log <- -log10(wtcc$GenePvalue)
wtcc$Color_Dis <- wtcc$Chr %% 2
wtcc$Color_Dis <-ifelse(wtcc$GenePvalue< 0.05, (wtcc$Chr %% 2)+2, wtcc$Color_Dis)

colours <- c("#D3D3D3","#808080",brewer.pal(n = 3, name = "Set1"))

pdf("Manhattan_plot.pdf",width=12,height=6)
p <- ggplot(wtcc, aes(Start, Discovery_log))+geom_point(size=1.5,alpha=0.6,aes(colour=as.factor(Color_Dis)))+
  scale_colour_manual(values = colours) +
  #scale_color_brewer(palette="Set1")+
  geom_hline(yintercept=-log10(0.05),size=0.5, colour="gray")+
  ylab(expression(paste(-log[10]~'P value')))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(name = "Chromosome", breaks = unique(wtcc$Tick), labels = unique(wtcc$Chr),expand=c(0.01,0))
print(p)
dev.off()
