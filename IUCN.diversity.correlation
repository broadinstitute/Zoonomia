z=read.table("diversity_analysis.input.txt", header=T)

library(ordPens)

###MATRIX PREP
for.het<-data.frame(as.numeric(z[,9]), as.factor(as.integer(z[,2])),as.factor(as.integer(z[,3])))
for.soh<-data.frame(as.numeric(z[,10]), as.factor(as.integer(z[,2])),as.factor(as.integer(z[,3])))


#het and soh can't be NA, so remove any species missing value for the relevant analysis
for.het<-for.het[!is.na(for.het[, 9]), ]
for.soh<-for.soh[!is.na(for.soh[, 10]), ]

#exclude species for which IUCN is missing
for.het.iucn<-for.het[!is.na(for.het[, 9]), ]
for.soh.iucn<-for.soh[!is.na(for.soh[, 10]), ]


#infer relationship between IUCN status and het, soh | these data
het.iucn.coeff<- ordSmooth(as.numeric(for.het.iucn[,2]), for.het.iucn[,1], offset = rep(0,length(for.het.iucn[,2])), lambda=1000, model = c("linear"))$coef
soh.iucn.coeff<-ordSmooth(as.numeric(for.soh.iucn[,2]), for.soh.iucn[,1], offset = rep(0,length(for.soh.iucn[,2])), lambda=1000, model = c("linear"))$coef

####create null to determine relationship if conservation categories instead assigned at random
real.iucn<-for.het.iucn[,2]
it<-10000
reorder.iucn<-matrix(0, length(real.iucn),it)
randomize.iucn<- for (i in 1:it)

{
	#print(i)
	random.iucn<-sample(real.iucn)
	#print(random.iucn)
#reorder.iucn[,i]<-random.iucn
test.coeff<-ordSmooth(as.numeric(random.iucn), for.het.iucn[,1], offset = rep(0,length(for.het.iucn[,2])), lambda=1000, model = c("linear"))$coef
reorder.iucn[,i]<-mean(test.coeff[3:6])
}
