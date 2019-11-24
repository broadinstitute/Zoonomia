library(tidyverse)
library(ggplot2)

outfile <- "diversity_analysis.correlation.output.txt"
file.remove(outfile)

pantheria <- read.delim("diversity_analysis.pantheria.phenotypes.txt",header=T)
phenos <- read.delim("diversity_analysis.pantheria.phenotype_types.txt",header=T)
div <- read.delim("diversity_analysis.input.txt",header=T)
d <- merge(pantheria,div,by="Species",all=T)
d <- d[!is.na(d$Species),]
all <- d

### use Anova test for categorical pantheria phenotypes
cols <- phenos[phenos$type=="categorical",]
cols <- as.vector(cols$phenotype)
for (col in cols) {
  tmp <- d[c("het","soh",col)]
  colnames(tmp) <- c("het","soh","val")
  tmp <- tmp[!is.na(tmp$val),]
  
  write(paste("\n##### ANOVA: het vs",col,"#####",sep=" "),file=outfile,append=T)
  capture.output(length(na.omit(tmp$het)),file=outfile,append=T)
  fit <- aov(het~val, data=tmp)
  capture.output(summary(fit),file=outfile,append=T)
  write(paste("\n##### ANOVA: soh vs",col,"#####",sep=" "),file=outfile,append=T)
  capture.output(length(na.omit(tmp$soh)),file=outfile,append=T)
  fit <- aov(soh~val, data=tmp)
  capture.output(summary(fit),file=outfile,append=T)
}
d[cols] <- NULL
### use Linear regression for continuous pantheria phenotypes
cols <- phenos[phenos$type=="continuous",]
cols <- as.vector(cols$phenotype)
for (col in cols) {
  tmp <- d[c("het","soh",col)]
  colnames(tmp) <- c("het","soh","value")
  tmp <- tmp[!is.na(tmp$value),]
  write(paste("\n##### LM: het vs",col,"#####",sep=" "),file=outfile,append=T)
  capture.output(length(na.omit(tmp$het)),file=outfile,append=T)
  capture.output(summary(lm(het~value,data=tmp)),file=outfile,append=T)
  write(paste("\n##### LM: soh vs",col,"#####",sep=" "),file=outfile,append=T)
  capture.output(length(na.omit(tmp$soh)),file=outfile,append=T)
  capture.output(summary(lm(soh~value,data=tmp)),file=outfile,append=T)
}
d[cols] <- NULL
d <- d[d$OK==TRUE&d$IUCN=="Least Concern",]

count <- d %>% count(Family)
count <- data.frame("Family"=count$Family,"NFAM"=count$n)
d <- merge(d,count,by="Family",all=T)
count <- d %>% count(Order)
count <- data.frame("Order"=count$Order,"NORDER"=count$n)
d <- merge(d,count,by="Order",all=T)
d <- unique(d)

tmp <- d[d$NFAM>1,]
write(paste("\n##### ANOVA: het vs family #####",sep=" "),file=outfile,append=T)
capture.output(length(na.omit(tmp$het)),file=outfile,append=T)
fit <- aov(het~Family, data=tmp)
capture.output(summary(fit),file=outfile,append=T)

write(paste("\n##### ANOVA: soh vs family #####",sep=" "),file=outfile,append=T)
capture.output(length(na.omit(tmp$soh)),file=outfile,append=T)
fit <- aov(soh~Family, data=tmp)
capture.output(summary(fit),file=outfile,append=T)

tmp <- d[d$NORDER>3,]
write(paste("\n##### ANOVA: het vs order #####",sep=" "),file=outfile,append=T)
capture.output(length(na.omit(tmp$het)),file=outfile,append=T)
fit <- aov(het~Order, data=tmp)
capture.output(summary(fit),file=outfile,append=T)

write(paste("\n##### ANOVA: soh vs order #####",sep=" "),file=outfile,append=T)
capture.output(length(na.omit(tmp$soh)),file=outfile,append=T)
fit <- aov(soh~Order, data=tmp)
capture.output(summary(fit),file=outfile,append=T)

d <- all[all$OK==TRUE,]
tmp <- d[d$POP=="wild"|d$POP=="captive",]

write(paste("\n##### T-TEST: het vs pop #####",sep=" "),file=outfile,append=T)
tmp <- tmp[!is.na(tmp$het),]
cnt <- tmp %>% count(POP)
capture.output(print(cnt),file=outfile,append=T)
capture.output(t.test(tmp$het~tmp$POP),file=outfile,append=T)

write(paste("\n##### T-TEST: SoH vs pop #####",sep=" "),file=outfile,append=T)
tmp <- tmp[!is.na(tmp$soh),]
cnt <- tmp %>% count(POP)
capture.output(print(cnt),file=outfile,append=T)
capture.output(t.test(tmp$soh~tmp$POP),file=outfile,append=T)

tmp <- tmp[tmp$IUCN=="Least Concern",]

write(paste("\n##### T-TEST: LC het vs pop #####",sep=" "),file=outfile,append=T)
tmp <- tmp[!is.na(tmp$het),]
cnt <- tmp %>% count(POP)
capture.output(print(cnt),file=outfile,append=T)
capture.output(t.test(tmp$het~tmp$POP),file=outfile,append=T)


write(paste("\n##### T-TEST: LC SoH vs pop #####",sep=" "),file=outfile,append=T)
tmp <- tmp[!is.na(tmp$soh),]
cnt <- tmp %>% count(POP)
capture.output(print(cnt),file=outfile,append=T)
capture.output(t.test(tmp$soh~tmp$POP),file=outfile,append=T)

colors <- c("#000000","#b2182b","#525252","#1f78b4")
conv <- data.frame("IUCN"=c("Critically Endangered","Endangered","Near Threatened","Vulnerable","Least Concern","Data deficient"),"nIUCN"=c(5,4,2,3,1,0),"shIUCN"=c("EN","EN","VU","VU","LC","DD"))
all <- merge(all,conv,by="IUCN",all=T)
all$Rhet <- (round(all$het*5000))/5000
all$Rsoh <- (round(all$soh*40))/40

d <- all[all$OK==TRUE,]
d <- d[!is.na(d$IUCN),]
#d <- d[d$IUCN!="Data deficient"]
d$IUCN <- factor(d$IUCN)

write(paste("\n##### het COUNTS: IUCN #####",sep=" "),file=outfile,append=T)
tmp <- d
tmp <- tmp[!is.na(tmp$het),]
cnt <- tmp %>% count(IUCN)
capture.output(print(cnt),file=outfile,append=T)
tmp <- tmp[order(tmp$nIUCN,tmp$Rhet,tmp$shIUCN,tmp$het),]
tmp$N <- c(1:length(tmp$Species))
p <- ggplot(tmp,aes(x=N,y=Rhet)) + geom_point(aes(colour=shIUCN,shape=POP),alpha=0.6) +  facet_grid(.~nIUCN, scales = "free_x", space = "free_x")  + scale_shape_manual(values=c(1,2,3,0)) + scale_colour_manual(values=colors) + scale_y_continuous(breaks=c(0,0.002,0.004,0.006,0.008),limits=c(0,0.01)) + theme_bw()
ggsave(plot=p,filename="diversity_analysis.figure3A.het.pdf",width=10.5,height=6)


write(paste("\n##### soh COUNTS: IUCN #####",sep=" "),file=outfile,append=T)
tmp <- d
tmp <- tmp[!is.na(tmp$soh),]
cnt <- tmp %>% count(IUCN)
capture.output(print(cnt),file=outfile,append=T)
tmp <- tmp[order(tmp$nIUCN,tmp$Rsoh,tmp$shIUCN,tmp$soh),]
tmp$N <- c(1:length(tmp$Species))
p <- ggplot(tmp,aes(x=N,y=Rsoh)) + geom_point(aes(colour=shIUCN,shape=POP),alpha=0.6) +  facet_grid(.~nIUCN, scales = "free_x", space = "free_x")  + scale_shape_manual(values=c(1,2,3,0)) + scale_colour_manual(values=colors) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1)) + theme_bw()
ggsave(plot=p,filename="diversity_analysis.figure3B.soh.pdf",width=10.5,height=6)

d <- all[all$OK==TRUE,]
d <- d[!is.na(d$POP),]
d <- d[d$POP=="captive"|d$POP=="wild",]
d$POP <- factor(d$POP)
d <- d[order(d$POP,d$Rhet,d$shIUCN,d$het),]
d$N <- c(1:length(d$Species))

## make figure 3C
cap <- d[d$POP=="captive",]
wild <- d[d$POP=="wild",]
mean <- data.frame("POP"=c("wild","captive"),"het"=c(mean(wild$het),mean(cap$het)),"soh"=c(mean(na.omit(wild$soh)),mean(na.omit(cap$soh))))
p <- ggplot(d,aes(x=N,y=Rhet)) + geom_point(aes(colour=shIUCN,shape=POP),alpha=0.6) + geom_hline(aes(yintercept=het),data=mean) + facet_wrap(~POP,nrow=1,scales="free_x") + scale_shape_manual(values=c(1,0)) + scale_colour_manual(values=colors) + scale_y_continuous(breaks=c(0,0.002,0.004,0.006,0.008),limits=c(0,0.01)) + theme_bw()
ggsave(plot=p,filename="diversity_analysis.figure3C.het.pdf",width=10.5,height=6)

## make figure 3D
d <- d[order(d$POP,d$Rsoh,d$shIUCN,d$soh),]
d$N <- c(1:length(d$Species))
p <- ggplot(d,aes(x=N,y=Rsoh)) + geom_point(aes(colour=shIUCN,shape=POP),alpha=0.6) + geom_hline(aes(yintercept=soh),data=mean) + facet_wrap(~POP,nrow=1,scales="free_x") + scale_shape_manual(values=c(1,0)) + scale_colour_manual(values=colors) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1)) +  theme_bw()
ggsave(plot=p,filename="diversity_analysis.figure3D.soh.pdf",width=10.5,height=6)

## make figure 3E
d <- all
lc <- d[d$OK==TRUE&d$shIUCN=="LC",]
cr <- d[d$IUCN=="Critically Endangered",]
p <- ggplot(d,aes(x=het,y=soh)) + geom_point(aes(colour=shIUCN,shape=POP),alpha=0.6)
p <- p + geom_text(aes(label=Common.Name),color="#b2182b",size=2,hjust="left",nudge_x=0.00012,data=cr)
p <- p + geom_hline(yintercept=median(na.omit(lc$soh)),linetype=2) + geom_vline(xintercept=median(na.omit(lc$het)),linetype=2)
p <- p + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1)) + scale_x_continuous(breaks=c(0,0.002,0.004,0.006,0.008),limits=c(0,0.01))
p <- p + scale_shape_manual(values=c(1,2,3,0)) + scale_colour_manual(values=colors) + theme_bw() 
ggsave(plot=p,filename="diversity_analysis.figure3E.het_v_soh.pdf",width=9.5,height=8)

