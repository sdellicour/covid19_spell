library(ade4)
library(ape)
library(dismo)
library(exactextractr)
library(fields)
library(gbm)
library(gplots)
library(lubridate)
library(mapi)
library(mgcv)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(spdep)

# 1. Analyses of Sciensano data (generating tables)
# 2. Analyses at the province levels (hospitalisations)
# 3. Analyses of hospital catchment areas

writingFiles = FALSE
showingPlots = FALSE

zTransformation = function(x)
	{ 
		x = (x-mean(x))/sqrt(var(x))
	}

# 1. Analyses of Sciensano data (generating tables)

data_cases = read.csv("https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv")
data_hosps = read.csv("https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv")
data_death = read.csv("https://epistat.sciensano.be/Data/COVID19BE_MORT.csv")
data_tests = read.csv("https://epistat.sciensano.be/Data/COVID19BE_tests.csv")
dates1 = unique(data_hosps$DATE)[!is.na(unique(data_hosps$DATE))]
dates2 = unique(data_cases$DATE)[!is.na(unique(data_cases$DATE))]
dates = unique(c(as.character(dates1),as.character(dates2)))
dates = as.character(dates[order(dates)])
dates = dates[which(dates=="2020-03-13"):length(dates)]
columns = c("Date","Day","Tests","CasesN","CasesT","C5DDT","CasesTL","CasesRt",
			"HospNet","HospT","H5DT","HospTL","HospRt","HospN1","HospTN1",
			"ICUN","ICUT","I5DT","ICUTL","ICURt",
			"DeathN1","DeathT1","D5DT1","DeathTL1","DeathRt")
tab = matrix(nrow=length(dates), ncol=length(columns))
colnames(tab) = columns; tab[,"Date"] = dates; tab[,"Day"] = c(1:dim(tab)[1])
for (i in 1:dim(tab)[1])
	{
		index = which(as.character(data_tests[,"DATE"])==tab[i,"Date"])
		tab[i,"Tests"] = sum(data_tests[index,"TESTS_ALL"])
	}
dates = ymd(as.character(data_cases[,"DATE"]))
for (i in 1:dim(tab)[1])
	{
		indices1 = which(dates==ymd(tab[i,"Date"]))
		indices2 = which(dates<=ymd(tab[i,"Date"]))
		indices3 = which(dates<=(ymd(tab[i,"Date"])-5))
		tab[i,"CasesN"] = sum(as.numeric(data_cases[indices1,"CASES"]))
		tab[i,"CasesT"] = sum(as.numeric(data_cases[indices2,"CASES"]))
		QTD1 = sum(as.numeric(data_cases[indices3,"CASES"]))
		D = 5; QTD2 = as.numeric(tab[i,"CasesT"]); DT = (D*log(2))/(log(QTD2/QTD1))
		tab[i,"C5DDT"] = DT; tab[i,"CasesTL"] = log10(as.numeric(tab[i,"CasesT"]))
		if (i != 1)
			{
				tab[i,"CasesRt"] = (as.numeric(tab[i-1,"CasesT"])+as.numeric(tab[i,"CasesN"]))/as.numeric(tab[i-1,"CasesT"])
			}
	}
dates = ymd(as.character(data_hosps[,"DATE"]))
for (i in 1:dim(tab)[1])
	{
		indices1 = which(dates==ymd(tab[i,"Date"]))
		indices2 = which(dates<=(ymd(tab[i,"Date"])-5))
		if (length(indices1) > 0)
			{
				tab[i,"HospT"] = sum(as.numeric(data_hosps[indices1,"TOTAL_IN"]))
				tab[i,"HospTL"] = log10(as.numeric(tab[i,"HospT"]))
				tab[i,"HospN1"] = sum(as.numeric(data_hosps[indices1,"NEW_IN"]))
				tab[i,"ICUT"] = sum(as.numeric(data_hosps[indices1,"TOTAL_IN_ICU"]))
				tab[i,"ICUTL"] = log10(as.numeric(tab[i,"ICUT"]))
				if ((i != 1)&(!is.na(tab[i-1,"HospT"])))
					{
						tab[i,"HospNet"] = as.numeric(tab[i,"HospT"])-as.numeric(tab[i-1,"HospT"])
						tab[i,"HospRt"] = (as.numeric(tab[i-1,"HospT"])+as.numeric(tab[i,"HospNet"]))/as.numeric(tab[i-1,"HospT"])
						tab[i,"HospTN1"] = as.numeric(tab[i-1,"HospTN1"])+as.numeric(tab[i,"HospN1"])
						tab[i,"ICUN"] = as.numeric(tab[i,"ICUT"])-as.numeric(tab[i-1,"ICUT"])
						tab[i,"ICURt"] = (as.numeric(tab[i-1,"ICUT"])+as.numeric(tab[i,"ICUN"]))/as.numeric(tab[i-1,"ICUT"])
					}	else	{
						tab[i,"HospTN1"] = as.numeric(tab[i,"HospN1"])
					}
				if (length(indices2) > 0)
					{
						QTD1 = sum(as.numeric(data_hosps[indices2,"TOTAL_IN"]))
						D = 5; QTD2 = as.numeric(tab[i,"HospT"]); DT1 = (D*log(2))/(log(QTD2/QTD1))
						QTD1 = sum(as.numeric(data_hosps[indices2,"TOTAL_IN_ICU"]))
						D = 5; QTD2 = as.numeric(tab[i,"ICUT"]); DT2 = (D*log(2))/(log(QTD2/QTD1))
						tab[i,"H5DT"] = DT1; tab[i,"I5DT"] = DT2
					}
			}
	}
dates = ymd(as.character(data_death[,"DATE"]))
for (i in 1:dim(tab)[1])
	{
		indices1 = which(dates==ymd(tab[i,"Date"]))
		indices2 = which(dates<=ymd(tab[i,"Date"]))
		indices3 = which(dates<=(ymd(tab[i,"Date"])-5))
		tab[i,"DeathN1"] = sum(as.numeric(data_death[indices1,"DEATHS"]))
		tab[i,"DeathT1"] = sum(as.numeric(data_death[indices2,"DEATHS"]))
		QTD1 = sum(as.numeric(data_death[indices3,"DEATHS"]))
		D = 5; QTD2 = as.numeric(tab[i,"DeathT1"]); DT = (D*log(2))/(log(QTD2/QTD1))
		tab[i,"D5DT1"] = DT; tab[i,"DeathTL1"] = log10(as.numeric(tab[i,"DeathT1"]))
		if (i != 1)
			{
				tab[i,"DeathRt"] = (as.numeric(tab[i-1,"DeathT1"])+as.numeric(tab[i,"DeathN1"]))/as.numeric(tab[i-1,"DeathT1"])
			}
	}
write.csv(tab, "Last_data_Sciensano.csv", row.names=F, quote=F)
provinceNames = unique(data_hosps[,"PROVINCE"])
for (h in 1:length(provinceNames))
	{
		dates1 = unique(data_hosps$DATE)[!is.na(unique(data_hosps$DATE))]
		dates2 = unique(data_cases$DATE)[!is.na(unique(data_cases$DATE))]
		dates = unique(c(as.character(dates1),as.character(dates2)))
		dates = as.character(dates[order(dates)])
		dates = dates[which(dates=="2020-03-13"):length(dates)]
		columns = c("Date","Day","CasesN","CasesT","C5DDT","CasesTL","CasesRt",
					"HospNet","HospT","H5DT","HospTL","HospRt","HospN1","HospTN1",
					"ICUN","ICUT","I5DT","ICUTL","ICURt")
		tab = matrix(nrow=length(dates), ncol=length(columns))
		colnames(tab) = columns; tab[,"Date"] = dates; tab[,"Day"] = c(1:dim(tab)[1])
		dates = ymd(as.character(data_cases[,"DATE"]))
		for (i in 1:dim(tab)[1])
			{
				indices1 = which((dates==ymd(tab[i,"Date"]))&(data_cases[,"PROVINCE"]==provinceNames[h]))
				indices2 = which((dates<=ymd(tab[i,"Date"]))&(data_cases[,"PROVINCE"]==provinceNames[h]))
				indices3 = which((dates<=(ymd(tab[i,"Date"])-5))&(data_cases[,"PROVINCE"]==provinceNames[h]))
				tab[i,"CasesN"] = sum(as.numeric(data_cases[indices1,"CASES"]))
				tab[i,"CasesT"] = sum(as.numeric(data_cases[indices2,"CASES"]))
				QTD1 = sum(as.numeric(data_cases[indices3,"CASES"]))
				D = 5; QTD2 = as.numeric(tab[i,"CasesT"]); DT = (D*log(2))/(log(QTD2/QTD1))
				tab[i,"C5DDT"] = DT; tab[i,"CasesTL"] = log10(as.numeric(tab[i,"CasesT"]))
				if (i != 1)
					{
						tab[i,"CasesRt"] = (as.numeric(tab[i-1,"CasesT"])+as.numeric(tab[i,"CasesN"]))/as.numeric(tab[i-1,"CasesT"])
					}
			}
		dates = ymd(as.character(data_hosps[,"DATE"]))
		for (i in 1:dim(tab)[1])
			{
				indices1 = which((dates==ymd(tab[i,"Date"]))&(data_hosps[,"PROVINCE"]==provinceNames[h]))
				indices2 = which((dates<=(ymd(tab[i,"Date"])-5))&(data_hosps[,"PROVINCE"]==provinceNames[h]))
				if (length(indices1) > 0)
					{
						tab[i,"HospT"] = sum(as.numeric(data_hosps[indices1,"TOTAL_IN"]))
						tab[i,"HospTL"] = log10(as.numeric(tab[i,"HospT"]))
						tab[i,"HospN1"] = sum(as.numeric(data_hosps[indices1,"NEW_IN"]))
						tab[i,"ICUT"] = sum(as.numeric(data_hosps[indices1,"TOTAL_IN_ICU"]))
						tab[i,"ICUTL"] = log10(as.numeric(tab[i,"ICUT"]))
						if ((i != 1)&(!is.na(tab[i-1,"HospT"])))
							{
								tab[i,"HospNet"] = as.numeric(tab[i,"HospT"])-as.numeric(tab[i-1,"HospT"])
								tab[i,"HospRt"] = (as.numeric(tab[i-1,"HospT"])+as.numeric(tab[i,"HospNet"]))/as.numeric(tab[i-1,"HospT"])
								tab[i,"HospTN1"] = as.numeric(tab[i-1,"HospTN1"])+as.numeric(tab[i,"HospN1"])
								tab[i,"ICUN"] = as.numeric(tab[i,"ICUT"])-as.numeric(tab[i-1,"ICUT"])
								tab[i,"ICURt"] = (as.numeric(tab[i-1,"ICUT"])+as.numeric(tab[i,"ICUN"]))/as.numeric(tab[i-1,"ICUT"])
							}	else	{
								tab[i,"HospTN1"] = as.numeric(tab[i,"HospN1"])
							}
						if (length(indices2) > 0)
							{
								QTD1 = sum(as.numeric(data_hosps[indices2,"TOTAL_IN"]))
								D = 5; QTD2 = as.numeric(tab[i,"HospT"]); DT1 = (D*log(2))/(log(QTD2/QTD1))
								QTD1 = sum(as.numeric(data_hosps[indices2,"TOTAL_IN_ICU"]))
								D = 5; QTD2 = as.numeric(tab[i,"ICUT"]); DT2 = (D*log(2))/(log(QTD2/QTD1))
								tab[i,"H5DT"] = DT1; tab[i,"I5DT"] = DT2
							}
					}
			}
		dates = ymd(as.character(data_death[,"DATE"]))
		write.csv(tab, paste0("Last_data_",gsub("Li\xe8ge","Liege",provinceNames[h]),".csv"), row.names=F, quote=F)
	}

# 2. Analyses at the province levels (hospitalisations)

selectedDays1 = ymd(c("2020-07-10","2020-07-17","2020-07-24","2020-07-31",
					  "2020-08-07","2020-08-14","2020-08-21","2020-08-28",
					  "2020-09-04","2020-09-11","2020-09-18","2020-09-25",
					  "2020-10-07","2020-10-14","2020-10-21","2020-10-28",
					  "2020-11-03"))
firstDay = ymd("2020-01-30"); D = 7 # time interval
selectedDays2 = as.numeric(selectedDays1-firstDay)
provinces = raster::getData("GADM", country="BEL", level=2)
provinces@data$NAME_3 = c("Brussels","Antwerpen","Limburg","OostVlaanderen","VlaamsBrabant",
				 "WestVlaanderen","BrabantWallon","Hainaut","Liège","Luxembourg","Namur")
data = read.csv("https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv")
data = data[!is.na(data[,"DATE"]),]; firstDay = ymd("2020-01-30")
data$DAYS = as.numeric(ymd(data[,"DATE"])-firstDay)
cumulatedHosCases = matrix(nrow=dim(provinces@data)[1]+1, ncol=selectedDays2[length(selectedDays2)])
cumulatedICUCases = matrix(nrow=dim(provinces@data)[1]+1, ncol=selectedDays2[length(selectedDays2)])
doublingTHosCases = matrix(nrow=dim(provinces@data)[1]+1, ncol=selectedDays2[length(selectedDays2)])
doublingTICUCases = matrix(nrow=dim(provinces@data)[1]+1, ncol=selectedDays2[length(selectedDays2)])
for (i in 1:dim(provinces@data)[1])
	{
		provincesID = provinces@data[i,"NAME_3"]
		lines = which(data[,"PROVINCE"] == provincesID)
		temp = data[lines,c("DAYS","TOTAL_IN","TOTAL_IN_ICU")]
		for (j in 1:length(selectedDays2))
			{
				index1 = which(temp[,"DAYS"]==(selectedDays2[j]-D))
				index2 = which(temp[,"DAYS"]==selectedDays2[j])
				if (length(index1) == 0)	
					{
						QTD1 = 0
					}	else	{
						QTD1 = temp[index1,"TOTAL_IN"]
					}
				QTD2 = temp[index2,"TOTAL_IN"]
				if ((QTD1 >= QTD2)|(QTD1 == 0)|(QTD2 < 10))
					{
						DT = NA
					}	else		{
						DT = (D*log(2))/(log(QTD2/QTD1))
					}
				provinces@data[i,paste0("DT",j)] = DT
			}
		for (j in 1:selectedDays2[length(selectedDays2)])
			{
				index = which(temp[,"DAYS"]==j)
				if (length(index) > 0)
					{
						cumulatedHosCases[i,j] = temp[index,"TOTAL_IN"]
						cumulatedICUCases[i,j] = temp[index,"TOTAL_IN_ICU"]
					}	else	{
						cumulatedHosCases[i,j] = 0
						cumulatedICUCases[i,j] = 0
					}
			}
	}
for (i in 1:dim(cumulatedHosCases)[2])
	{
		cumulatedHosCases[dim(cumulatedHosCases)[1],i] = sum(cumulatedHosCases[1:(dim(cumulatedHosCases)[1]-1),i])
		cumulatedICUCases[dim(cumulatedICUCases)[1],i] = sum(cumulatedICUCases[1:(dim(cumulatedICUCases)[1]-1),i])
	}

if (showingPlots)
	{
		periods = c("","","","","","","","","","","","","","","15/10-21/10/2020","22/10-28/10/2020","29/10-03/11/2020")
		DTmax = ceiling(max(provinces@data[,c("DT10","DT11","DT12")])); DTmax = 50
		colourScale = colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101]; cols = list()
		dev.new(width=3.2,height=7); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 15:length(selectedDays2))
			{
				values = provinces@data[,paste0("DT",i)]
				values[which(values>DTmax)] = DTmax
				cols[[i]] = colourScale[1+((values/DTmax)*100)]
				plot(provinces, border="gray30", col=cols[[i]])
				mtext(paste0("Doubling time - ",periods[i]), cex=0.50, col="gray30", at=3.55, line=-14.2)
				plot(legendRast, legend.only=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
	 				 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
	 				 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.2,0), at=seq(0,DTmax,10), labels=c("0","10","20","30","40","50")))
			}
		DT_values = c(0,1,2,4,8,16,32); DTmax = max(DT_values); cols = list()
		colourScale1 = colorRampPalette(brewer.pal(9,"YlGn"))(12)[1:length(DT_values)]
		colourScale1 = colorRampPalette(brewer.pal(11,"RdYlGn"))(12)[3:(length(DT_values)+2)]
		colourScale2 = c("gray90",colourScale1)
		dev.new(width=9,height=2.7); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 15:length(selectedDays2))
			{
				values1 = provinces@data[,paste0("DT",i)]; values2 = provinces@data[,paste0("DT",i)]
				values2[values1[]>32] = 7; values2[values1[]<=32] = 6; values2[values1[]<=16] = 5
				values2[values1[]<=8] = 4; values2[values1[]<=4] = 3; values2[values1[]<=2] = 2
				values2[values1[]<=1] = 1; values2[is.na(values1)] = 0; cols[[i]] = colourScale2[values2+1]
				plot(provinces, border="gray30", col=cols[[i]])
				mtext(paste0("Doubling time hospitalisations"), cex=0.5, col="gray30", at=3.55, line=-14.2)
				mtext(paste0(periods[i]), cex=0.6, col="gray30", at=3.55, line=-15.2)
			}
		legend(6.3, 51.5, c("0 - 1 day","1 - 2 days","2 - 4 days","4 - 8 days","8 - 16 days","16 - 32 days",">32 days"),
			   col=colourScale1, text.col="gray30", pch=16, pt.cex=1.5, box.lty=0, cex=0.9, y.intersp=1.2)
	}

for (i in 1:dim(cumulatedHosCases)[1])
	{
		for (j in (D+1):dim(cumulatedHosCases)[2])
			{
				QTD1 = cumulatedHosCases[i,j-D]
				QTD2 = cumulatedHosCases[i,j]
				DT = (D*log(2))/(log(QTD2/QTD1))
				doublingTHosCases[i,j] = DT
				QTD1 = cumulatedICUCases[i,j-D]
				QTD2 = cumulatedICUCases[i,j]
				DT = (D*log(2))/(log(QTD2/QTD1))
				doublingTICUCases[i,j] = DT
			}
	}
row.names(cumulatedHosCases) = c(provinces@data$NAME_2,"Belgium")
row.names(cumulatedICUCases) = c(provinces@data$NAME_2,"Belgium")
colnames(cumulatedHosCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
colnames(cumulatedICUCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
cumulatedHosCases = cumulatedHosCases[,45:selectedDays2[length(selectedDays2)]]
cumulatedICUCases = cumulatedICUCases[,45:selectedDays2[length(selectedDays2)]]
row.names(doublingTHosCases) = c(provinces@data$NAME_2,"Belgium")
row.names(doublingTICUCases) = c(provinces@data$NAME_2,"Belgium")
colnames(doublingTHosCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
colnames(doublingTICUCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
doublingTHosCases = doublingTHosCases[,52:selectedDays2[length(selectedDays2)]]
doublingTICUCases = doublingTICUCases[,52:selectedDays2[length(selectedDays2)]]
if (writingFiles) write.csv(cumulatedHosCases, "Hospitalisations_provinces.csv")
if (writingFiles) write.csv(cumulatedICUCases, "SoinsIntensifs_provinces.csv")
if (writingFiles) write.csv(doublingTHosCases, "Hospitalisations_dT_provinces.csv")
if (writingFiles) write.csv(doublingTICUCases, "SoinsIntensifs_dT_provinces.csv")

if (showingPlots)
	{
		selectedDays = as.numeric(dmy(c(paste0(c(21:31),"-10-2020"),paste0(c(01:03),"-11-2020")))-firstDay)
		doublingTHosCases_selected = doublingTHosCases[,paste0("day_",selectedDays)]
		xLabels = c(paste0(c(21:31),"-10"),paste0(c(01:03),"-11")); dates = c(1:length(xLabels))
		DTmax = 50; cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#d3d3d3","#4d4d4d")
		dev.new(width=7,height=5); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(1,1), mar=c(2.9,3.1,1,1), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:dim(doublingTHosCases_selected)[1])
			{
				if (i == 1)
					{
						plot(doublingTHosCases_selected[i,], col=cols[i], lwd=1, ylim=c(1.8,33), axes=F, ann=F, type="l")
					}	else	{
						lines(doublingTHosCases_selected[i,], col=cols[i], lwd=1)
					}
				points(dates,doublingTHosCases_selected[i,], col=cols[i], cex=0.7, pch=16)
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.013, col.axis="gray30", mgp=c(0,0.05,0), at=c(-1,dates), labels=c(-1,xLabels))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-2,32,2))
		title(ylab="doubling time hospitalisation (time window = 7 days)", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		title(xlab="day", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
		legend(1, 33, provinces@data$NAME_2, col=cols, text.col="gray30", pch=16, pt.cex=1.2, box.lty=0, cex=0.7, y.intersp=1.3)
	}

# 3. Analyses of hospital catchment areas

communes1 = shapefile("Shapefile_communes/Shapefile_NIS5_codes.shp")
communes2 = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
communes2 = spTransform(communes2, proj4string(communes1))
equivalence = read.csv("Shapefile_communes/Postal_codes_vs_NIS.csv", header=T)
catchmentAreas = shapefile("Hosp_catchmentArea/Hospital_catchment_areas_080420.shp")
catchmentAreas@data[,"X_ID"] = gsub(" ","",catchmentAreas@data[,"X_ID"])
catchmentAreas@data$area = as.numeric(catchmentAreas@data$area)

	# 3.1. Establishing the link between catchment areas and communes

computeSharedAreas = FALSE
if (computeSharedAreas)
	{
		sharedAreas = matrix(nrow=dim(communes1@data)[1], ncol=length(catchmentAreas@polygons))
		row.names(sharedAreas) = communes1@data[,"NIS5"]; row.names(populations) = communes1@data[,"NIS5"]
		totalArea1 = 0; totalArea2 = 0 # N.B.: sharedAreas not used anymore
		for (i in 1:length(catchmentAreas@polygons))
			{
				for (j in 1:length(catchmentAreas@polygons[[i]]@Polygons))
					{
						pol = catchmentAreas@polygons[[i]]@Polygons[[j]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol1 = sps; proj4string(pol1) = crs(catchmentAreas)
						totalArea1 = totalArea1 + pol1@polygons@area
					}
				for (j in 1:dim(communes1@data)[1])
					{
						area = 0
						for (k in 1:length(communes1@polygons[[j]]@Polygons))
							{
								pol = communes1@polygons[[j]]@Polygons[[k]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol2 = sps; proj4string(pol2) = crs(communes1)
								totalArea2 = totalArea2 + pol@area
								if (!is.null(raster::intersect(pol1, pol2)))
									{
										pol3 = intersect(pol1, pol2)
										for (k in 1:length(pol3@polygons))
											{
												area = area + pol3@polygons[[k]]@area
											}
									}				
							}
						sharedAreas[j,i] = area
					}
			}
	}
population_WP = raster("WorldPop_pop_raster.tif")
communes1_pop = spTransform(communes1, population_WP@crs)
communes2_pop = spTransform(communes2, population_WP@crs)
catchmentAreas_pop = spTransform(catchmentAreas, population_WP@crs)
catchmentAreas@data$population = rep(0,dim(catchmentAreas@data)[1])
for (i in 1:length(catchmentAreas_pop@polygons))
	{
		population = 0
		for (j in 1:length(catchmentAreas_pop@polygons[[i]]@Polygons))
			{
				pol = catchmentAreas_pop@polygons[[i]]@Polygons[[j]]
				p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol1 = sf::st_as_sfc(sps); st_crs(pol1) = 4326; crs(population_WP)=crs(pol1)
				population = population + exact_extract(population_WP, pol1, fun='sum')
			}
		catchmentAreas@data[i,"population"] = population
	}
if ((!file.exists("Catchment_pop_1.csv"))|(!file.exists("Catchment_pop_2.csv")))
	{
		populations1 = matrix(nrow=dim(communes1@data)[1], ncol=length(catchmentAreas@polygons))
		populations2 = matrix(nrow=dim(communes2@data)[1], ncol=length(catchmentAreas@polygons))
		row.names(populations1) = communes1@data[,"NIS5"]
		for (i in 1:length(catchmentAreas_pop@polygons))
			{
				population = 0
				for (j in 1:length(catchmentAreas_pop@polygons[[i]]@Polygons))
					{
						pol = catchmentAreas_pop@polygons[[i]]@Polygons[[j]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol1 = sps; proj4string(pol1) = crs(catchmentAreas_pop)
						population = population + exact_extract(population_WP, sf::st_as_sfc(pol1), fun='sum')
					}
				catchmentAreas@data[i,"population"] = population
				for (j in 1:dim(communes1@data)[1])
					{
						population = 0
						for (k in 1:length(catchmentAreas_pop@polygons[[i]]@Polygons))
							{
								pol = catchmentAreas_pop@polygons[[i]]@Polygons[[k]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol1 = sps; proj4string(pol1) = crs(catchmentAreas_pop)
								for (l in 1:length(communes1_pop@polygons[[j]]@Polygons))
									{
										pol = communes1_pop@polygons[[j]]@Polygons[[l]]
										p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
										pol2 = sps; proj4string(pol2) = crs(communes1_pop)
										if (!is.null(raster::intersect(pol1,pol2)))
											{
												pol3 = intersect(pol1,pol2)
												for (k in 1:length(pol3@polygons))
													{
														population = population + exact_extract(population_WP, sf::st_as_sfc(pol3), fun='sum')
													}
											}				
									}
							}
						populations1[j,i] = population
					}
				for (j in 1:dim(communes2@data)[1])
					{
						population = 0
						for (k in 1:length(catchmentAreas_pop@polygons[[i]]@Polygons))
							{
								pol = catchmentAreas_pop@polygons[[i]]@Polygons[[k]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol1 = sps; proj4string(pol1) = crs(catchmentAreas_pop)
								for (l in 1:length(communes2_pop@polygons[[j]]@Polygons))
									{
										pol = communes2_pop@polygons[[j]]@Polygons[[l]]
										p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
										pol2 = sps; proj4string(pol2) = crs(communes2_pop)
										if (!is.null(raster::intersect(pol1,pol2)))
											{
												pol3 = intersect(pol1,pol2)
												for (k in 1:length(pol3@polygons))
													{
														population = population + exact_extract(population_WP, sf::st_as_sfc(pol3), fun='sum')
													}
											}				
									}
							}
						populations2[j,i] = population
					}
			}
		write.csv(populations1, "Catchment_pop_1.csv", quote=F)
		write.csv(populations2, "Catchment_pop_2.csv", quote=F)
	}
populations1 = read.csv("Catchment_pop_1.csv", head=T)
populations2 = read.csv("Catchment_pop_2.csv", head=T)
proportions1 = matrix(nrow=dim(communes1@data)[1], ncol=length(catchmentAreas@polygons))
proportions2 = matrix(nrow=dim(communes2@data)[1], ncol=length(catchmentAreas@polygons))
for (i in 1:dim(proportions1)[1])
	{
		for (j in 1:dim(proportions1)[2]) proportions1[i,j] = populations1[i,j]/sum(populations1[i,])
	}
for (i in 1:dim(proportions2)[1])
	{
		for (j in 1:dim(proportions2)[2]) proportions2[i,j] = populations2[i,j]/sum(populations2[i,])
	}
row.names(proportions1) = communes1@data[,"NIS5"]

	# 3.2. Extracting and assigning covariate values to each commune

data = read.csv("https://epistat.sciensano.be/Data/COVID19BE_CASES_MUNI_CUM.csv")
communes1@data$cases = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which(data[,"NIS5"]==communes1@data[i,"NIS5"])
		if (length(index) != 1)
			{
				# cat(i,"\n")
			}	else		{
				if (as.character(data[index,"CASES"]) != "<5")
					{
						communes1@data[i,"cases"] = as.numeric(as.character(data[index,"CASES"]))
					}
			}
	}
data = read.csv("Data_SPF_Economie/SPF_total_population.csv")
communes1@data$population = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else	{
				communes1@data[i,"population"] = data[index,"TOTAL"]
			}
	}
communes1@data$popDensity = communes1@data$population/(communes1@data$Shape_Area/(10^6))
communes1@data$populationLog = log(communes1@data$population)
communes1@data$popDensityLog = log(communes1@data$popDensity)
data = read.csv("Data_SPF_Economie/SPF_pop_median_age.csv")
communes1@data$medianAge = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes1@data[i,"medianAge"] = data[index,"AGE_MEDIAN"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_more_than_65yrs.csv")
communes1@data$moreThan65 = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which((data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])&(data[,"MS_SEX"]=="TOTAL"))
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes1@data[i,"moreThan65"] = data[index,"X..65year"]/data[index,"TOTAL"]
			}
	}
data = read.csv("MRs_BE_17-04-20.csv", head=T)
communes1@data$maisonsDeRepos = rep(0,dim(communes1@data)[1])
communes1@data$bedsInMRs = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which(data[,"NIS5"]==communes1@data[i,"NISCode"])
		if (length(index) == 1)
			{
				communes1@data[i,"maisonsDeRepos"] = data[index,"n_mr"]
				communes1@data[i,"bedsInMRs"] = data[index,"n_tot_beds"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_median_incomes.csv")
communes1@data$medianIncome = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes1@data[i,"medianIncome"] = data[index,"MEDIAN_DECL"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_working_sectors.csv")
communes1@data$sectorP = rep(0,dim(communes1@data)[1])
communes1@data$sectorS = rep(0,dim(communes1@data)[1])
communes1@data$sectorT = rep(0,dim(communes1@data)[1])
for (i in 1:dim(communes1@data)[1])
	{
		index = which((data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])&(data[,"CD_SECT"]=="P"))
		if (length(index) == 1) communes1@data[i,"sectorP"] = data[index,"MS_PROP"]*100
		index = which((data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])&(data[,"CD_SECT"]=="S"))
		if (length(index) == 1) communes1@data[i,"sectorS"] = data[index,"MS_PROP"]*100
		index = which((data[,"CD_REFNIS"]==communes1@data[i,"NISCode"])&(data[,"CD_SECT"]=="T"))
		if (length(index) == 1) communes1@data[i,"sectorT"] = data[index,"MS_PROP"]*100
	}

	# 3.3. 	Extracting and assigning covariate values to each catchment area

catchmentAreas@data$xCentroid = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$yCentroid = rep(0,dim(catchmentAreas@data)[1])
for (i in 1:dim(catchmentAreas@data)[1])
	{
		maxArea = 0; polIndex = 0
		for (j in 1:length(catchmentAreas@polygons[[i]]@Polygons))
			{
				if (maxArea < catchmentAreas@polygons[[i]]@Polygons[[j]]@area)
					{
						maxArea = catchmentAreas@polygons[[i]]@Polygons[[j]]@area; polIndex = j
					}
			}
		pol = catchmentAreas@polygons[[i]]@Polygons[[polIndex]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes1@proj4string
		centroidCoordinates = coordinates(pol)
		catchmentAreas@data[i,"xCentroid"] = centroidCoordinates[1,1]
		catchmentAreas@data[i,"yCentroid"] = centroidCoordinates[1,2]
	}
population_WP = raster("WorldPop_pop_raster.tif")
catchmentAreas_pop = spTransform(catchmentAreas, population_WP@crs)
catchmentAreas@data$population = rep(0,dim(catchmentAreas@data)[1])
for (i in 1:length(catchmentAreas_pop@polygons))
	{
		population = 0
		for (j in 1:length(catchmentAreas_pop@polygons[[i]]@Polygons))
			{
				pol = catchmentAreas_pop@polygons[[i]]@Polygons[[j]]
				p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol1 = sf::st_as_sfc(sps); st_crs(pol1) = 4326; crs(population_WP) = crs(pol1)
				population = population + exact_extract(population_WP, pol1, fun='sum')
			}
		catchmentAreas@data[i,"population"] = population
	}
catchmentAreas@data$popDensity = catchmentAreas@data$population/(catchmentAreas@data$area/(10^6))
catchmentAreas@data$popDensityLog = log(catchmentAreas@data$popDensity)
catchmentAreas@data$medianAge = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$moreThan65 = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$maisonsDeRepos = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$bedsInMRs = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$medianIncome = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$sectorP = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$sectorS = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$sectorT = rep(0,dim(catchmentAreas@data)[1])
for (i in 1:dim(catchmentAreas@data)[1])
	{
		catchmentAreas@data[i,"medianAge"] = sum(communes1@data[,"medianAge"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"moreThan65"] = sum(communes1@data[,"moreThan65"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"medianIncome"] = sum(communes1@data[,"medianIncome"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"maisonsDeRepos"] = sum(communes1@data[,"maisonsDeRepos"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"bedsInMRs"] = sum(communes1@data[,"bedsInMRs"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"sectorP"] = sum(communes1@data[,"sectorP"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"sectorS"] = sum(communes1@data[,"sectorS"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"sectorT"] = sum(communes1@data[,"sectorT"]*populations1[,i])/catchmentAreas@data[i,"population"]
		catchmentAreas@data[i,"cases"] = sum(communes1@data[,"cases"]*populations1[,i])/catchmentAreas@data[i,"population"]
		if (catchmentAreas@data[i,"medianAge"] == 0) print(i)
	}
pm10 = raster("PM10_mean_2017.asc")
pm25 = raster("PM25_mean_2017.asc")
parks = raster("Urban_parcs_UNamur.tif")
catchmentAreas@data$pm10 = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$pm25 = rep(0,dim(catchmentAreas@data)[1])
catchmentAreas@data$ratioParkPop = rep(0,dim(catchmentAreas@data)[1])
for (i in 1:dim(catchmentAreas@data)[1])
	{
		maxArea = 0; polIndex = 0
		for (j in 1:length(catchmentAreas@polygons[[i]]@Polygons))
			{
				if (maxArea < catchmentAreas@polygons[[i]]@Polygons[[j]]@area)
					{
						maxArea = catchmentAreas@polygons[[i]]@Polygons[[j]]@area; polIndex = j
					}
			}
		pol = catchmentAreas@polygons[[i]]@Polygons[[polIndex]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sf::st_as_sfc(sps); st_crs(pol) = communes1@proj4string
		crs(pm10) = crs(pol); crs(pm25) = crs(pol); crs(parks) = crs(pol)
		catchmentAreas@data[i,"pm10"] = exact_extract(pm10, pol, fun='mean')
		catchmentAreas@data[i,"pm25"] = exact_extract(pm25, pol, fun='mean')
		catchmentAreas@data[i,"ratioParkPop"] = exact_extract(parks, pol, fun='sum')/catchmentAreas@data[i,"population"]
	}
if (!file.exists("Hosp_catchmentArea/Proportion_of_urban_areas_080420.csv"))
	{
		clc = raster("CorineLandCover.asc")
		catchmentAreas_clc = spTransform(catchmentAreas, raster("CorineLandCover18.tif")@crs)
		propUrbanArea = matrix(nrow=dim(catchmentAreas@data)[1], ncol=3)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				maxArea = 0; polIndex = 0
				for (j in 1:length(catchmentAreas_clc@polygons[[i]]@Polygons))
					{
						if (maxArea < catchmentAreas_clc@polygons[[i]]@Polygons[[j]]@area)
							{
								maxArea = catchmentAreas_clc@polygons[[i]]@Polygons[[j]]@area; polIndex = j
							}
					}
				pol = catchmentAreas_clc@polygons[[i]]@Polygons[[polIndex]]
				p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol = sps; proj4string(pol) = catchmentAreas_clc@proj4string
				pol_light = gSimplify(pol, 100)
				rast = mask(crop(clc,pol_light),pol_light)
				greenAreas = sum(rast[]==141, na.rm=T)
				urbanAreas = sum(rast[]==111, na.rm=T)+sum(rast[]==112, na.rm=T)
				propUrbanArea[i,1] = i
				if (greenAreas == 0)
					{
						propUrbanArea[i,2] = 0
					}	else	{
						propUrbanArea[i,2] = catchmentAreas@data[i,"population"]/greenAreas
					}
				propUrbanArea[i,3] = urbanAreas/length(rast[!is.na(rast[])])
			}
		colnames(propUrbanArea) = c("NIS","popGreenArea","propUrbanArea")
		write.csv(propUrbanArea, "Hosp_catchmentArea/Proportion_of_urban_areas_080420.csv", row.names=F, quote=F)
	}
catchmentAreas@data$propUrbanArea = read.csv("Hosp_catchmentArea/Proportion_of_urban_areas_080420.csv")[,3]

	# 3.4. Cumputing doubling times for hospitalisations and ICU

data = read.csv("Raw_data_Sciensano/Hosp_surge_overview_19092020_1500.csv", head=T, sep=";")
days1 = c(paste0("2020-07-0",c(1:9)),paste0("2020-07-",c(10:31)),paste0("2020-08-0",c(1:9)),paste0("2020-08-",c(10:31)),
		  paste0("2020-09-0",c(1:9)),paste0("2020-09-",c(10:18))); D = 7 # time interval
firstDay = ymd("2020-01-30"); days2 = ymd(days1); days3 = as.numeric(days2-firstDay)
data$date = rep(NA, dim(data)[1])
for (i in 1:dim(data)[1])
	{
		data[i,"date"] = as.character(dmy(unlist(strsplit(as.character(data[i,"Date"])," "))[1]))
	}
data$days = as.numeric(ymd(data[,"date"])-firstDay)
newHospitaliCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
cumulatedHosCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
cumulatedICUCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
incidenceHosCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
incidenceICUCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
doublingTHosCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
doublingTICUCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
growthRHosCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
growthRICUCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
for (i in 1:dim(catchmentAreas@data)[1])
	{
		hospitalIDs = unlist(strsplit(catchmentAreas@data[i,"X_ID"],"-"))
		hospitalIDs = as.numeric(hospitalIDs); lines = c()
		for (j in 1:length(hospitalIDs))
			{
				lines = c(lines,which(as.numeric(data[,"ERK.AGR"])==hospitalIDs[j]))
			}
		temp1 = data[lines,c("days","NewPatientsNotReferredHospital","NewPatientsNotReferredHospitalNursingHome",
							 "Confirmed.patients.in.hospital","Confirmed.patients.in.ICU")]
		temp1 = temp1[order(temp1[,"days"]),]; temp2 = temp1[1,]
		for (j in 2:dim(temp1)[1])
			{
				if (temp1[j,"days"] == temp2[dim(temp2)[1],"days"])
					{
						for (k in 2:dim(temp2)[2])
							{
								temp2[dim(temp2)[1],k] = sum(c(temp2[dim(temp2)[1],k],temp1[j,k]), na.rm=T)
							}
					}	else	{
						temp2 = rbind(temp2, temp1[j,])
					}
			}
		temp3 = temp2[1,]; d1 = temp2[1,"days"]+1; d2 = temp2[dim(temp2)[1],"days"]
		for (d in d1:d2)
			{
				index = which(temp2[,"days"]==d)
				if (length(index) == 1)
					{
						temp3 = rbind(temp3, temp2[index,])
					}	else	{
						temp3 = rbind(temp3, temp3[dim(temp3)[1],])
						temp3[dim(temp3)[1],"days"] = d 
					}
			}
		temp = temp3
		for (j in 1:length(days3))
			{
				index1 = which(temp[,"days"]==(days3[j]-D))				
				index2 = which(temp[,"days"]==days3[j])
				if (length(index2) > 0)
					{
						newHospitaliCases[i,j] = temp[index2,"NewPatientsNotReferredHospital"]+temp[index2,"NewPatientsNotReferredHospitalNursingHome"]
						if (length(index1) == 0)	
							{
								QTD1 = 0
							}	else	{
								QTD1 = temp[index1,"Confirmed.patients.in.hospital"]
							}
						QTD2 = temp[index2,"Confirmed.patients.in.hospital"]
						DT = (D*log(2))/(log(QTD2/QTD1))
						cumulatedHosCases[i,j] = QTD2
						incidenceHosCases[i,j] = (QTD2/(catchmentAreas@data[i,"population"]))*1000
						doublingTHosCases[i,j] = DT
						growthRHosCases[i,j] = QTD2/QTD1
						if (length(index1) == 0)	
							{
								QTD1 = 0
							}	else	{
								QTD1 = temp[index1,"Confirmed.patients.in.ICU"]
							}
						QTD2 = temp[index2,"Confirmed.patients.in.ICU"]
						DT = (D*log(2))/(log(QTD2/QTD1))
						cumulatedICUCases[i,j] = QTD2
						incidenceICUCases[i,j] = (QTD2/(catchmentAreas@data[i,"population"]))*1000
						doublingTICUCases[i,j] = DT
						growthRICUCases[i,j] = QTD2/QTD1
					}
			}
	}
colnames(newHospitaliCases) = paste0("newHospitaliCases_",days1)
colnames(cumulatedHosCases) = paste0("cumulatedHosCases_",days1)
colnames(incidenceHosCases) = paste0("incidenceHosCases_",days1)
colnames(doublingTHosCases) = paste0("doublingTHosCases_",days1)
colnames(growthRHosCases) = paste0("growthRHosCases_",days1)
colnames(cumulatedICUCases) = paste0("cumulatedICUCases_",days1)
colnames(incidenceICUCases) = paste0("incidenceICUCases_",days1)
colnames(doublingTICUCases) = paste0("doublingTICUCases_",days1)
colnames(growthRICUCases) = paste0("growthRICUCases_",days1)
catchmentAreas@data = cbind(catchmentAreas@data, newHospitaliCases)
catchmentAreas@data = cbind(catchmentAreas@data, cumulatedHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, incidenceHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, doublingTHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, growthRHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, cumulatedICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, incidenceICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, doublingTICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, growthRICUCases)

	# 3.5. Saving and plotting the variables assigned to each area

if (writingFiles)
	{
		df = catchmentAreas@data
		indices = which(grepl("cumulated",colnames(df)))
		for (i in 1:length(indices)) df[which(is.na(df[,indices[i]])),indices[i]] = 0
		indices = which(grepl("incidence",colnames(df)))
		for (i in 1:length(indices)) df[which(is.na(df[,indices[i]])),indices[i]] = 0
		colnames(df) = gsub("cumulatedH","h",colnames(df)); colnames(df) = gsub("cumulatedI","I",colnames(df))
		write.csv(df, "All_figures_&_outputs/Covariables_catchment_areas_210920.csv", row.names=F, quote=F)
	}
if (showingPlots)
	{
		selectedDays1 = ymd(c("2020-09-04","2020-09-11","2020-09-18"))
		periods = c("29/08-04/09/2020","05/09-11/09/2020","12/09-18/09/2020")
		variables = c("incidenceHosCases_2020-09-04","incidenceHosCases_2020-09-11","incidenceHosCases_2020-09-18",
					  "doublingTHosCases_2020-09-04","doublingTHosCases_2020-09-11","doublingTHosCases_2020-09-18",
					  "incidenceICUCases_2020-09-04","incidenceICUCases_2020-09-11","incidenceICUCases_2020-09-18",
					  "doublingTICUCases_2020-09-04","doublingTICUCases_2020-09-11","doublingTICUCases_2020-09-18")
		variableNames1 = c("Hosp. incidence 04/09/20","Hosp. incidence 11/09/20","Hosp. incidence 18/09/20",
						  "Hosp. doubling time","Hosp. doubling time","Hosp. doubling time",
						  "ICU incidence 04/09/20","ICU incidence 11/09/20","ICU incidence 18/09/20",
						  "ICU doubling time","ICU doubling time","ICU doubling time")
		variableNames2 = c("(# cases/1000 persons)","(# cases/1000 persons)","(# cases/1000 persons)",periods[1],periods[2],periods[3],
						  "(# cases/1000 persons)","(# cases/1000 persons)","(# cases/1000 persons)",periods[1],periods[2],periods[3])
		colourScales = list()
		colourScales[[1]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[2]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[3]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[4]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[5]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[6]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[7]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(121)[1:101])
		colourScales[[8]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(121)[1:101])
		colourScales[[9]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(121)[1:101])
		colourScales[[10]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		colourScales[[11]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		colourScales[[12]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				values = catchmentAreas@data[,variables[i]]
				minV = 0; maxV = max(values)
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]
				if (i%in%c(4,5,6,10,11,12))
					{
						values[is.infinite(values)] = 0; values[is.na(values)] = 0
						minV = 0; values[values[]<minV] = minV
						maxV = 30; values[values[]>maxV] = maxV
					}
				if (i%in%c(1,2,3))
					{
						maxV = 2; values[values[]>maxV] = maxV
					}
				if (i%in%c(7,8,9))
					{
						maxV = 0.3; values[values[]>maxV] = maxV
					}
				legendCols = legendCols[2:length(legendCols)]
				legendRast = raster(as.matrix(c(minV,maxV)))
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(catchmentAreas, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames1[i], cex=0.54, col="gray30", at=92000, line=-11.9)
				mtext(variableNames2[i], cex=0.54, col="gray30", at=92000, line=-12.6)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}

		variables = c("growthRHosCases_2020-09-04","growthRHosCases_2020-09-11","growthRHosCases_2020-09-18",
					  "doublingTHosCases_2020-09-04","doublingTHosCases_2020-09-11","doublingTHosCases_2020-09-18",
					  "growthRICUCases_2020-09-04","growthRICUCases_2020-09-11","growthRICUCases_2020-09-18",
					  "doublingTICUCases_2020-09-04","doublingTICUCases_2020-09-11","doublingTICUCases_2020-09-18")
		variableNames1 = c("Hosp. growth rate 04/09/20","Hosp. growth rate 11/09/20","Hosp. growth rate 18/09/20",
						  "Hosp. doubling time","Hosp. doubling time","Hosp. doubling time",
						  "ICU growthR 04/09/20","ICU growthR 11/09/20","ICU growthR 18/09/20",
						  "ICU doubling time","ICU doubling time","ICU doubling time")
		variableNames2 = c("(# cases/1000 persons)","(# cases/1000 persons)","(# cases/1000 persons)",periods[1],periods[2],periods[3],
						  "(# cases/1000 persons)","(# cases/1000 persons)","(# cases/1000 persons)",periods[1],periods[2],periods[3])
		colourScales = list()
		colourScales[[1]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdBu"))(161)[31:131])
		colourScales[[2]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdBu"))(161)[31:131])
		colourScales[[3]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdBu"))(161)[31:131])
		colourScales[[4]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[5]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[6]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"RdYlGn"))(161)[31:131])
		colourScales[[7]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"PRGn"))(161)[31:131])
		colourScales[[8]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"PRGn"))(161)[31:131])
		colourScales[[9]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"PRGn"))(161)[31:131])
		colourScales[[10]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		colourScales[[11]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		colourScales[[12]] = c("#E5E5E5",colorRampPalette(brewer.pal(11,"BrBG"))(161)[31:131])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				values = catchmentAreas@data[,variables[i]]; minV = 0
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]
				if (i%in%c(4,5,6,10,11,12))
					{
						values[is.infinite(values)] = 0; values[is.na(values)] = 0
						minV = 0; values[values[]<minV] = minV
						maxV = 30; values[values[]>maxV] = maxV
					}
				if (i%in%c(1,2,3))
					{
						values[is.infinite(values)] = 0; values[is.na(values)] = 0
						maxV = 10; values[values[]>maxV] = maxV
					}
				if (i%in%c(7,8,9))
					{
						values[is.infinite(values)] = 0; values[is.na(values)] = 0
						maxV = 10; values[values[]>maxV] = maxV
					}
				legendCols = legendCols[2:length(legendCols)]
				legendRast = raster(as.matrix(c(minV,maxV)))
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(catchmentAreas, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames1[i], cex=0.54, col="gray30", at=92000, line=-11.9)
				mtext(variableNames2[i], cex=0.54, col="gray30", at=92000, line=-12.6)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}
	}

	# 3.6. Performing and plotting the first axes of an exploratory PCA

if (showingPlots)
	{
		variableNames = c("population","popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
						  "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea")
		df = catchmentAreas@data[,variableNames]
		pca = dudi.pca(df, scannf=F, nf=length(variableNames)); lis = pca$li[,1:2]; cos = pca$co
		dev.new(width=6, height=6); par(mar=c(3,3,1.5,1.5), lwd=0.2, col="gray30")
		plot(lis, col="gray50", cex=0.3, pch=16, ann=F, axes=F, xlim=c(-10.0,6.5), ylim=c(-4.0,5.0))
		points(lis, col="gray30", cex=0.75, pch=1, lwd=0.3); points(lis, col="gray70", cex=0.70, pch=16); 
		s.corcircle(2*cos, xax=1, yax=2, box=F, sub="", csub=0.7, clabel=0.7, possub="topleft", grid=F, cgrid=1, full=F, add.plot=T)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-12,10,2))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-7,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
	}

	# 3.7. Assessing spatial autocorrelation with the Moran's I test

variableNames = c("xCentroid","yCentroid","popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
				  "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea",
				  "incidenceHosCases_2020-09-04","incidenceHosCases_2020-09-11","incidenceHosCases_2020-09-18",
				  "growthRHosCases_2020-09-04","growthRHosCases_2020-09-11","growthRHosCases_2020-09-18")
df = catchmentAreas@data[,variableNames]; colnames(df) = gsub("-","",colnames(df))
responseVariables = c("incidenceHosCases_20200904","incidenceHosCases_20200911","incidenceHosCases_20200918")
# responseVariables = c("growthRHosCases_20200904","growthRHosCases_20200911","growthRHosCases_20200918")
for (i in 1:length(responseVariables))
	{
		values = df[,responseVariables[i]]
		indices = which((!is.na(values))&(!is.infinite(values)))
		geoDists = as.matrix(dist(df[indices,c("xCentroid","yCentroid")]))
		weights = 1/geoDists; diag(weights) = 0
		print(Moran.I(values[indices],weights)$p.value)
	}

	# 3.8. Univariate (LR) followed by multivariate regression (GLM) analyses

selectedVariables = list()
for (i in 1:length(responseVariables))
	{
		predictors = c("popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
			   		   "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea")
		buffer = c()
		for (j in 1:length(predictors))
			{
				values = df[,responseVariables[i]]; tmp = df
				tmp = tmp[which((!is.na(values))&(!is.infinite(values))),]
				formula = paste0(responseVariables[i]," ~ ",predictors[j])
				lr = glm(formula, data=tmp)
				pValue = summary(lr)$coefficients[2,4]
				if (pValue < 0.05)
					{
						buffer = c(buffer, predictors[j])
					}
			}
		selectedVariables[[i]] = buffer
		if (is.null(buffer)) selectedVariables[[i]] = NA
	}
for (i in 1:length(responseVariables))
	{
		if (!is.na(selectedVariables[[i]][1]))
			{
				values = df[,responseVariables[i]]; tmp = df
				tmp = tmp[which((!is.na(values))&(!is.infinite(values))),]
				formula = paste0(responseVariables[i]," ~ ",selectedVariables[[i]][1])
				if (length(selectedVariables[[i]]) > 1)
					{
						for (j in 2:length(selectedVariables[[i]]))
							{
								formula = paste0(formula," + ",selectedVariables[[i]][j])
							}
					}
				glm = glm(formula, data=tmp); print(summary(glm)); cat("\n\n")
			}
	}

	# 3.9. GAM (generalised additive model) analyses (NOT ACTUALISED YET)

gams = list(); zTransformations = FALSE
for (i in 1:length(responseVariables))
	{
		values = df[,responseVariables[i]]; tmp = df
		tmp = tmp[which((!is.na(values))&(!is.infinite(values))),]
		colnames(tmp) = gsub(responseVariables[i],"responseVariable",colnames(tmp))
		if (zTransformations == TRUE)
			{
				for (j in 1:dim(dfs[[i]])[2]) df[,j] = zTransformation(df[,j])
			}
		gam = gam(responseVariable ~ s(maisonsDeRepos) + s(bedsInMRs) + s(sectorT) + s(cases) + s(pm10) + s(pm25), data=tmp, method="REML")
		print(summary(gam)); cat("\n\n"); gams[[i]] = gam
		if (showingPlots)
			{
				dev.new(); plot(gam, pages=1)
			}
	}
if (showingPlots)
	{
		gam = gams[[3]]; responseCurves = list()
		curves = plot(gam, pages=1); dev.off()
		selectedVariables = c("maisonsDeRepos")
		variableNames = c("nurshing homes")
		dev.new(width=3.5,height=3)
		par(mfrow=c(1,1), mar=c(3,3,1,2), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, col="gray30", bty="o")
		for (i in 1:length(selectedVariables))
			{
				index = NA
				for (j in 1:length(curves))
					{
						if (curves[[j]]$xlab==selectedVariables[i]) index = j
					}
				lower_l = curves[[index]]$fit-curves[[index]]$se
				upper_l = curves[[index]]$fit+curves[[index]]$se
				yLim = c(min(c(lower_l,upper_l)),max(c(lower_l,upper_l)))
				xx_l = c(curves[[index]]$x,rev(curves[[index]]$x)); yy_l = c(lower_l,rev(upper_l))
				plot(curves[[index]]$x, curves[[index]]$fit, ylim=yLim, ann=F, axes=F, type="l", col="gray30", lwd=1.0)
				polygon(xx_l, yy_l, col=rgb(100,100,100,100,maxColorValue=255), border=0)
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				title(xlab=variableNames[i], cex.lab=0.7, mgp=c(0.9,0,0), col.lab="gray30")
				title(ylab=paste0("s(",variableNames[i],")"), cex.lab=0.7, mgp=c(1.1,0,0), col.lab="gray30")
			}
		dev.new(width=3.5,height=3)
		par(mfrow=c(1,1), mar=c(3,3,1,2), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, col="gray30", bty="o")
		for (i in 1:length(selectedVariables))
			{
				tmp = df[,c("sectorT","pm10")]
				# tmp = cbind(df[5:dim(df)[2]],df[,c("xCentroid","yCentroid")])
				for (j in 1:dim(tmp)[2])
					{
						if (colnames(tmp)[j] == selectedVariables[i])
							{
								tmp[,j] = seq(min(tmp[,j]),max(tmp[,j]),(max(tmp[,j])-min(tmp[,j]))/(dim(tmp)[1]-1))
							}	else	{
								tmp[,j] = median(tmp[,j])
							}
					}
				plot(tmp[,selectedVariables[i]], predict(gam, tmp), ann=F, axes=F, type="l", col="gray30", lwd=1.0)
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				title(xlab=variableNames[i], cex.lab=0.7, mgp=c(0.9,0,0), col.lab="gray30")
				title(ylab=paste0("response"), cex.lab=0.7, mgp=c(1.1,0,0), col.lab="gray30")
			}
	}

	# 3.10. Classic correlation analyses

variableNames = c("incidenceHosCases_2020-09-04","incidenceHosCases_2020-09-11","incidenceHosCases_2020-09-18",
				  "xCentroid","yCentroid","popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
				  "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea")
df = catchmentAreas@data[,variableNames]; colnames(df) = gsub("-","",colnames(df))
correlations = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
pValues = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
for (i in 2:dim(correlations)[2])
	{
		for (j in 1:(i-1))
			{
				vS = df[,c(i,j)]
				vS[is.infinite(vS[,1]),1] = NA
				vS[is.infinite(vS[,2]),1] = NA
				vS = vS[which((!is.na(vS[,1]))&(!is.na(vS[,2]))),]	
				test = cor.test(vS[,1],vS[,2])
				correlations[j,i] = test$estimate
				pValues[j,i] = test$p.value
			}
	}
texts = correlations
for (i in 2:dim(texts)[1])
	{
		for (j in 1:(i-1))
			{
				if (pValues[j,i] < 0.05)
					{
						texts[j,i] = as.character(round(correlations[j,i],1))
					}	else	{
						texts[j,i] = NA
					}
			}
	}
if (showingPlots)
	{
		cols = colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:111]
		minV = min(correlations, na.rm=T); maxV = max(correlations, na.rm=T)
		index1 = (50+(minV*50))+1; index2 = (50+(maxV*50))+1; cols = cols[index1:index2]
		heatmap.2(correlations, cellnote=texts, notecex=0.6, notecol="gray30", main=NULL, density.info="none",
				  trace="none", margins=c(12,9), col=cols, Rowv=NULL, Colv="NA", dendrogram="none", key=F,
				  labRow=variableNames, labCol=variableNames, cexRow=0.8, cexCol=0.8, offsetRow=0.0,
				  colRow=rep("gray30",length(variableNames)), colCol=rep("gray30",length(variableNames)))
	}

