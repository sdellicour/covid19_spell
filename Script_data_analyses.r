library(lubridate)
library(raster)
library(RColorBrewer)
library(rgeos)

writingFiles = FALSE

# 1. Doubling time analyses and maps

periods = c("15-20/03/2020","22-27/03/2020","26-31/03/2020")
selectedDays1 = ymd(c("2020-03-20","2020-03-27","2020-03-31"))
firstDay = ymd("2020-01-30"); D = 5 # time interval
selectedDays2 = as.numeric(selectedDays1-firstDay)

provinces = getData("GADM", country="BEL", level=2)
provinces@data$NAME_3 = c("Brussels","Antwerpen","Limburg","OostVlaanderen","VlaamsBrabant",
				 "WestVlaanderen","BrabantWallon","Hainaut","Li\xe8ge","Luxembourg","Namur")
data = read.csv("Data_Sciensano_0104/COVID19BE_HOSP_20200401.csv")
data = data[!is.na(data[,"DATE"]),]; firstDay = ymd("2020-01-30")
data$DAYS = as.numeric(ymd(data[,"DATE"])-firstDay)
cumulatedHosCases = matrix(nrow=dim(provinces@data)[1], ncol=selectedDays2[length(selectedDays2)])
cumulatedICUCases = matrix(nrow=dim(provinces@data)[1], ncol=selectedDays2[length(selectedDays2)])
doublingTHosCases = matrix(nrow=dim(provinces@data)[1], ncol=selectedDays2[length(selectedDays2)])
doublingTICUCases = matrix(nrow=dim(provinces@data)[1], ncol=selectedDays2[length(selectedDays2)])
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
if (writingFiles) write.csv(round(tab,2), "Doubling_times_provinces.csv", quote=F)
for (i in 1:dim(cumulatedHosCases)[1])
	{
		for (j in (7+1):dim(cumulatedHosCases)[2])
			{
				QTD1 = cumulatedHosCases[i,j-7]
				QTD2 = cumulatedHosCases[i,j]
				DT = (7*log(2))/(log(QTD2/QTD1))
				doublingTHosCases[i,j] = DT
				QTD1 = cumulatedICUCases[i,j-7]
				QTD2 = cumulatedICUCases[i,j]
				DT = (7*log(2))/(log(QTD2/QTD1))
				doublingTICUCases[i,j] = DT
			}
	}
row.names(cumulatedHosCases) = provinces@data$NAME_2
row.names(cumulatedICUCases) = provinces@data$NAME_2
colnames(cumulatedHosCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
colnames(cumulatedICUCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
cumulatedHosCases = cumulatedHosCases[,45:selectedDays2[length(selectedDays2)]]
cumulatedICUCases = cumulatedICUCases[,45:selectedDays2[length(selectedDays2)]]
row.names(doublingTHosCases) = provinces@data$NAME_2
row.names(doublingTICUCases) = provinces@data$NAME_2
colnames(doublingTHosCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
colnames(doublingTICUCases) = paste0("day_",seq(1,selectedDays2[length(selectedDays2)]))
doublingTHosCases = doublingTHosCases[,52:selectedDays2[length(selectedDays2)]]
doublingTICUCases = doublingTICUCases[,52:selectedDays2[length(selectedDays2)]]
if (writingFiles) write.csv(cumulatedHosCases, "Hospitalisations_provinces.csv")
if (writingFiles) write.csv(cumulatedICUCases, "SoinsIntensifs_provinces.csv")
if (writingFiles) write.csv(doublingTHosCases, "Hospitalisations_dT_provinces.csv")
if (writingFiles) write.csv(doublingTICUCases, "SoinsIntensifs_dT_provinces.csv")

cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#d3d3d3")
DTmax = ceiling(max(c(max(doublingTHosCases))))
dev.new(width=5,height=3.9); legendRast = raster(as.matrix(seq(0,DTmax,1)))
par(mfrow=c(1,1), mar=c(2.9,3.1,1,1), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
dates = c(22:31); xLabels = paste0(c(21:31),"-03")
for (i in 1:dim(doublingTHosCases)[1])
	{
		if (i == 1)
			{
				plot(dates,doublingTHosCases[i,], col=cols[i], lwd=1, ylim=c(1.8,10), axes=F, ann=F, type="l")
			}	else	{
				lines(dates,doublingTHosCases[i,], col=cols[i], lwd=1)
			}
		points(dates,doublingTHosCases[i,], col=cols[i], cex=0.7, pch=16)
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.013, col.axis="gray30", mgp=c(0,0.05,0), at=c(21:31), labels=xLabels)
axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.30,0), at=c(-1:10))
title(ylab="doubling time", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
title(xlab="day", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
legend(22.2, 10.28, provinces@data$NAME_2, col=cols, text.col="gray30", pch=16, pt.cex=1.0, box.lty=0, cex=0.65, y.intersp=1.1)

# 2. Preparation of different covariates to test

	# 2.1. Preparation at the commune level

communes = shapefile("Shapefile_communes/Shapefile_communes.shp")
data = read.csv("Data_Sciensano_0404/COVID19BE_CASES_MUNI_CUM.csv")
communes@data$cases = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which(data[,"NIS5"]==communes@data[i,"NIS5"])
		if (length(index) != 1)
			{
				# cat(i,"\n")
			}	else		{
				if (as.character(data[index,"CASES"]) != "<5")
					{
						communes@data[i,"cases"] = as.numeric(as.character(data[index,"CASES"]))
					}
			}
	}
data = read.csv("Data_SPF_Economie/SPF_total_population.csv")
communes@data$population = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else	{
				communes@data[i,"population"] = data[index,"TOTAL"]
			}
	}
communes@data$incidences = communes@data$cases/(communes@data$population/1000)
communes@data$popDensity = communes@data$population/(communes@data$Shape_Area/(10^6))
communes@data$populationLog = log(communes@data$population)
communes@data$popDensityLog = log(communes@data$popDensity)
data = read.csv("Data_SPF_Economie/SPF_pop_median_age.csv")
communes@data$medianAge = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes@data[i,"medianAge"] = data[index,"AGE_MEDIAN"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_more_than_65yrs.csv")
communes@data$moreThan65 = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which((data[,"CD_REFNIS"]==communes@data[i,"NISCode"])&(data[,"MS_SEX"]=="TOTAL"))
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes@data[i,"moreThan65"] = data[index,"X..65year"]/data[index,"TOTAL"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_median_incomes.csv")
communes@data$medianIncome = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which(data[,"CD_REFNIS"]==communes@data[i,"NISCode"])
		if (length(index) != 1)
			{
				cat(i,"\n")
			}	else		{
				communes@data[i,"medianIncome"] = data[index,"MEDIAN_DECL"]
			}
	}
data = read.csv("Data_SPF_Economie/SPF_working_sectors.csv")
communes_light = gSimplify(communes, 100)
communes@data$sectorP = rep(0,dim(communes@data)[1])
communes@data$sectorS = rep(0,dim(communes@data)[1])
communes@data$sectorT = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		index = which((data[,"CD_REFNIS"]==communes@data[i,"NISCode"])&(data[,"CD_SECT"]=="P"))
		if (length(index) == 1) communes@data[i,"sectorP"] = data[index,"MS_PROP"]*100
		index = which((data[,"CD_REFNIS"]==communes@data[i,"NISCode"])&(data[,"CD_SECT"]=="S"))
		if (length(index) == 1) communes@data[i,"sectorS"] = data[index,"MS_PROP"]*100
		index = which((data[,"CD_REFNIS"]==communes@data[i,"NISCode"])&(data[,"CD_SECT"]=="T"))
		if (length(index) == 1) communes@data[i,"sectorT"] = data[index,"MS_PROP"]*100
	}
if (!file.exists("PM10_anmean_2017.asc"))
	{
		pm10 = aggregate(raster("Rasters_de_irceline_be/PM10_anmean_2017_v564.tif"),100)
		pm25 = aggregate(raster("Rasters_de_irceline_be/PM25_anmean_2017_v564.tif"),100)
		writeRaster(pm10,"PM10_anmean_2017.asc"); writeRaster(pm25,"PM25_anmean_2017.asc")		
	}
pm10 = raster("PM10_anmean_2017.asc")
pm25 = raster("PM10_anmean_2017.asc")
communes@data$pm10 = rep(0,dim(communes@data)[1])
communes@data$pm25 = rep(0,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		maxArea = 0; polIndex = 0
		for (j in 1:length(communes@polygons[[i]]@Polygons))
			{
				if (maxArea < communes@polygons[[i]]@Polygons[[j]]@area)
					{
						maxArea = communes@polygons[[i]]@Polygons[[j]]@area; polIndex = j
					}
			}
		pol = communes@polygons[[i]]@Polygons[[polIndex]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pol_light = gSimplify(pol, 100)
		communes@data[i,"pm10"] = mean(raster::extract(pm10,pol_light)[[1]], na.rm=T)
		communes@data[i,"pm25"] = mean(raster::extract(pm25,pol_light)[[1]], na.rm=T)
	}
if (!file.exists("CorineLandCover.asc"))
	{
		clc = crop(raster("CorineLandCover18.tif"), extent(3500000,4500000,2500000,4000000))
		clc_crs = clc@crs; communes_clc = spTransform(communes, clc_crs)
		clc = mask(crop(clc, communes_clc), communes_clc)
		clc@crs = clc_crs; writeRaster(clc, "CorineLandCover.asc")
	}
if (!file.exists("CLC_propUrbanArea.csv"))
	{
		clc = raster("CorineLandCover.asc")
		communes_clc = spTransform(communes, raster("CorineLandCover18.tif")@crs)
		propUrbanArea = matrix(nrow=dim(communes@data)[1], ncol=3)
		for (i in 1:dim(communes@data)[1])
			{
				maxArea = 0; polIndex = 0
				for (j in 1:length(communes_clc@polygons[[i]]@Polygons))
					{
						if (maxArea < communes_clc@polygons[[i]]@Polygons[[j]]@area)
							{
								maxArea = communes_clc@polygons[[i]]@Polygons[[j]]@area; polIndex = j
							}
					}
				pol = communes_clc@polygons[[i]]@Polygons[[polIndex]]
				p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol = sps; proj4string(pol) = communes_clc@proj4string
				pol_light = gSimplify(pol, 100)
				rast = mask(crop(clc,pol_light),pol_light)
				greenAreas = sum(rast[]==141, na.rm=T)
				urbanAreas = sum(rast[]==111, na.rm=T)+sum(rast[]==112, na.rm=T)
				propUrbanArea[i,1] = communes@data[i,"NIS5"]
				if (greenAreas == 0)
					{
						propUrbanArea[i,2] = 0
					}	else	{
						propUrbanArea[i,2] = communes@data[i,"population"]/greenAreas
					}
				propUrbanArea[i,3] = urbanAreas/length(rast[!is.na(rast[])])
			}
		colnames(propUrbanArea) = c("NIS","popGreenArea","propUrbanArea")
		write.csv(propUrbanArea, "CLC_propUrbanArea.csv", row.names=F, quote=F)
	}
communes@data$propUrbanArea = read.csv("CLC_propUrbanArea.csv")[,3]

variables = c("incidences","popDensityLog","medianIncome","sectorP","sectorS","sectorT",
			  "medianAge","moreThan65","pm25")
variableNames = c("# cases per 1000 persons","Population density (log)","Median declared income (€)",
				  "% in primary sector","% in secundary sector","% in tertiary sector",
				  "Median age (years)",">= 65 years (proportion)","PM 2.5 emission")
communes_light = gSimplify(communes, 100); colourScales = list()
colourScales[[1]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101])
colourScales[[2]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(151)[1:101])
colourScales[[3]] = c(colorRampPalette(brewer.pal(9,"RdPu"))(151)[1:101])
colourScales[[4]] = c(colorRampPalette(brewer.pal(9,"Greens"))(151)[1:101])
colourScales[[5]] = c(colorRampPalette(brewer.pal(9,"Oranges"))(151)[1:101])
colourScales[[6]] = c(colorRampPalette(brewer.pal(9,"Blues"))(151)[1:101])
colourScales[[7]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
colourScales[[8]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
colourScales[[9]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(151)[1:101])
dev.new(width=7,height=6); par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:length(variables))
	{
		minV = min(communes@data[,variables[i]]); maxV = max(communes@data[,variables[i]])
		legendCols = colourScales[[i]][1:length(colourScales[[i]])]; legendRast = raster(as.matrix(c(minV,maxV)))		
		cols = colourScales[[i]][(((communes@data[,variables[i]]-minV)/(maxV-minV))*100)+1]
		plot(communes_light, border="gray30", col=cols, lwd=0.1)
		mtext(variableNames[i], cex=0.54, col="gray30", at=92000, line=-12.2)
		plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
	 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
			 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
	}

