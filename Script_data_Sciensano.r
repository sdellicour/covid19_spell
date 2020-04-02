library(lubridate)
library(raster)
library(RColorBrewer)

periods = c("15-20/03/2020","22-27/03/2020","25-30/03/2020")
selectedDays1 = ymd(c("2020-03-20","2020-03-27","2020-03-30"))
firstDay = ymd("2020-01-30"); D = 5 # time interval
selectedDays2 = as.numeric(selectedDays1-firstDay)

provinces = getData("GADM", country="BEL", level=2)
provinces@data$NAME_3 = c("Brussels","Antwerpen","Limburg","OostVlaanderen","VlaamsBrabant",
				 "WestVlaanderen","BrabantWallon","Hainaut","Li\xe8ge","Luxembourg","Namur")
data = read.csv("Data_Sciensano_3103/COVID19BE_HOSP.csv")
data = data[!is.na(data[,"DATE"]),]
data$DAYS = as.numeric(ymd(data[,"DATE"])-firstDay)
for (i in 1:dim(provinces@data)[1])
	{
		provincesID = provinces@data[i,"NAME_3"]
		lines = which(data[,"PROVINCE"] == provincesID)
		temp = data[lines,c("DAYS","TOTAL_IN","TOTAL_IN_ICU")]
		for (j in 1:length(selectedDays2))
			{
				if ((i == 1)&(j == 1)) provinces@data$DT1 = rep(NA,dim(provinces@data)[1])
				if ((i == 1)&(j == 2)) provinces@data$DT2 = rep(NA,dim(provinces@data)[1])
				if ((i == 1)&(j == 3)) provinces@data$DT3 = rep(NA,dim(provinces@data)[1])
				index1 = which(temp[,"DAYS"]==(selectedDays2[j]-D))
				index2 = which(temp[,"DAYS"]==selectedDays2[j])
				QTD1 = temp[index1,"TOTAL_IN"]
				QTD2 = temp[index2,"TOTAL_IN"]
				DT = (D*log(2))/(log(QTD2/QTD1))
				provinces@data[i,paste0("DT",j)] = DT
			}
	}
DTmax = ceiling(max(provinces@data[,c("DT1","DT2","DT3")]))
colourScale = colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101]; cols = list()
dev.new(width=3.2,height=7); legendRast = raster(as.matrix(seq(0,DTmax,1)))
par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:length(selectedDays2))
	{
		cols[[i]] = colourScale[1+((provinces@data[,paste0("DT",i)]/DTmax)*100)]
		plot(provinces, border="gray30", col=cols[[i]])
		mtext(paste0("Doubling time - ",periods[i]), cex=0.54, col="gray30", at=3.55, line=-14.2)
		plot(legendRast, legend.only=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
	 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
	 		 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.2,0), at=seq(0,DTmax,1)))
	}
tab = provinces@data[,c("DT1","DT2","DT3")]
row.names(tab) = provinces@data$NAME_2; colnames(tab) = periods
write.csv(round(tab,2), "Doubling_times_provinces.csv", quote=F)

