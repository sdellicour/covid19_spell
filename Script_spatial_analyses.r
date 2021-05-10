library(ade4)
library(ape)
library(caTools)
library(dismo)
library(exactextractr)
library(fields)
library(gbm)
library(gdistance)
library(gplots)
library(lubridate)
library(malariaAtlas)
library(mapi)
library(mgcv)
library(randomForest)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(sf)
library(spdep)

# 2. Analyses of Sciensano data (generating tables)
# 3. Analyses at the province levels (hospitalisations)
# 4. Analyses of hospital catchment areas (spatial)
# 5. Analyses of hospital catchment areas (temporal)

writingFiles = FALSE
showingPlots = FALSE

zTransformation = function(x)
	{ 
		x = (x-mean(x))/sqrt(var(x))
	}

# 1. Delimitating the hospital catchment areas (HCAs)

frictionR = readGDAL("Friction_raster_Belgium.tif")
hospitals = read_sf("All_hospitals_Belgium/All_hospitals_Belgium.shp")

	# 1.1. Computing the travel time to the closest hospital for the whole country

tr = transition(raster(frictionR), function(x) 1/mean(x), 8)
trG = geoCorrection(tr)
hospitals_xy = st_coordinates(hospitals)
hospitals_access = accCost(trG, hospitals_xy)
writeRaster(hospitals_access, "Access_to_hospitals.tif")

	# 1.2. For each hospital, estimating the travel time from each pixel of the map

s = stack()
for (i in 1:nrow(hospitals))
	{
		coord = st_coordinates(hospitals[i,])
		rast = accCost(trG, coord); s = stack(s, rast)
		print(paste0("Computation for hospital ",i))
	}
vect = getValues(hospitals_access)[!is.infinite(getValues(hospitals_access))]
mat = matrix(nrow=nrow(hospitals), ncol=length(vect))
for (i in 1:nrow(hospitals))
	{
		mat[i,] = getValues(s[[i]])[!is.infinite(getValues(s[[i]]))]
	}

	# 1.3. For each pixel, identifying which hospital is associated with the lowest travel time

vectsel = c()
for (i in 1:ncol(mat))
	{
		vectsel[i] = which(mat[,i]==min(mat[,i]))
	}
vectInf = values(hospitals_access); j = 1
for (i in 1:(ncol(hospitals_access)*nrow(hospitals_access)))
	{
		if (is.infinite(vectInf[i]) == FALSE)
			{
		 		vectInf[i] = vectsel[j]
				j = j+1
			}
	}
catchmentAreas = raster(vals=(vectInf), ext=extent(hospitals_access), crs=crs(hospitals_access), nrows=dim(hospitals_access)[1], ncols=dim(hospitals_access)[2])

	# 1.4. Converting and saving the resulting raster into a shapefile gathering HCAs

catchmentAreas = rasterToPolygons(catchmentAreas, dissolve=T, na.rm=T)
catchmentAreas = subset(catchmentAreas, is.finite(HCAs@data[,"layer"]))
metadata = matrix(nrow=dim(catchmentAreas@data)[1], ncol=2); colnames(metadata) = c("area","X_ID")
xS = as(hospitals,"Spatial")@coords[,1]; yS = as(hospitals,"Spatial")@coords[,2]
for (i in 1:dim(metadata)[1])
	{
		area = 0
		for (j in 1:length(catchmentAreas@polygons[[i]]@Polygons))
			{
				pol = catchmentAreas@polygons[[i]]@Polygons[[j]]
				indices = which(point.in.polygon(xS, yS, pol@coords[,1], pol@coords[,2]) == 1)
				if (length(indices) > 0)
					{
						metadata[i,"X_ID"] = hospitals$ID[indices[1]]
						if (length(indices) > 1)
							{
								for (k in 2:length(indices))
									{
										metadata[i,"X_ID"] = paste(metadata[i,"X_ID"],hospitals$ID[indices[k]],sep="-")
									}
							}
					}
				p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol = sps; proj4string(pol) = crs(catchmentAreas)
				area = area + raster::area(pol)
			}
		metadata[i,"area"] = round(area)
	}
catchmentAreas@data = as.data.frame(metadata)
writeOGR(catchmentAreas, dsn="Hosp_catchmentArea", "Hospital_catchment_areas_NEW", driver="ESRI Shapefile")

# 2. Analyses of Sciensano data (generating tables)

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

# 3. Analyses at the province levels (hospitalisations)

selectedDays1 = ymd(c("2020-07-10","2020-07-17","2020-07-24","2020-07-31",
					  "2020-08-07","2020-08-14","2020-08-21","2020-08-28",
					  "2020-09-04","2020-09-11","2020-09-18","2020-09-25",
					  "2020-10-07","2020-10-14","2020-10-21","2020-10-28",
					  "2020-11-04"))
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
		periods = c("","","","","","","","","","","","","","","15/10-21/10/2020","22/10-28/10/2020","29/10-04/11/2020")
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
		selectedDays = as.numeric(dmy(c(paste0(c(21:31),"-10-2020"),paste0(c(01:04),"-11-2020")))-firstDay)
		doublingTHosCases_selected = doublingTHosCases[,paste0("day_",selectedDays)]
		xLabels = c(paste0(c(21:31),"-10"),paste0(c(01:04),"-11")); dates = c(1:length(xLabels))
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

# 4. Analyses of hospital catchment areas (spatial)

communes1 = shapefile("Shapefile_communes/Shapefile_NIS5_codes.shp")
communes2 = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
communes2 = spTransform(communes2, proj4string(communes1))
equivalence = read.csv("Shapefile_communes/Postal_codes_vs_NIS.csv", header=T)
catchmentAreas = shapefile("Hosp_catchmentArea/Hospital_catchment_areas_080420.shp")
catchmentAreas@data[,"X_ID"] = gsub(" ","",catchmentAreas@data[,"X_ID"])
catchmentAreas@data$area = as.numeric(catchmentAreas@data$area)

	# 4.1. Establishing the link between catchment areas and communes

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
												pol3 = raster::intersect(pol1, pol2)
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
												pol3 = raster::intersect(pol1,pol2)
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
		write.csv(populations1, "Catchment_pop_1.csv", quote=F, row.names=F)
		write.csv(populations2, "Catchment_pop_2.csv", quote=F, row.names=F)
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

	# 4.2. Extracting and assigning covariate values to each commune

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

	# 4.3. 	Extracting and assigning covariate values to each catchment area

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
catchmentAreas@data[,"ratioBedsInMRsPopulation"] = catchmentAreas@data[,"bedsInMRs"]/catchmentAreas@data[,"population"]
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
		catchmentAreas@data[i,"pm10"] = exact_extract(pm10, pol, fun="mean")
		catchmentAreas@data[i,"pm25"] = exact_extract(pm25, pol, fun="mean")
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

	# 4.4. Assigning total numbers of new entries to each catchment area

selected_dates = dmy(c("31-05-2020","31-08-2020","30-11-2020"))
data = read.csv("Raw_data_Sciensano/Hosp_surge_overview_03122020_1500.csv", head=T, sep=";")
hospitalisations = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(selected_dates))
colnames(hospitalisations) = paste0("hospitalisation_",as.character(selected_dates))
hospitalisationsMR = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(selected_dates))
colnames(hospitalisationsMR) = paste0("hospitalisationMR_",as.character(selected_dates))
dates = dmy(gsub("\\/","-",data[,"Date"]))
for (i in 1:dim(hospitalisations)[1])
	{
		hospitalIDs = unlist(strsplit(catchmentAreas@data[i,"X_ID"],"-"))
		hospitalIDs = as.numeric(hospitalIDs); lines = c()
		for (j in 1:length(hospitalIDs))
			{
				lines = c(lines,which(as.numeric(data[,"ERK.AGR"])==hospitalIDs[j]))
			}
		sub = data[lines,]; sub_dates = dates[lines]
		for (j in 1:dim(hospitalisations)[2])
			{
				indices = which(sub_dates<=selected_dates[j])
				hospitalisations[i,j] = sum(sub[indices,"NewPatientsNotReferredHospital"], na.rm=T)
										# (ERROR:) + sum(sub[indices,"NewPatientsNotReferredHospitalNursingHome"], na.rm=T)
				hospitalisationsMR[i,j] = sum(sub[indices,"NewPatientsNotReferredHospitalNursingHome"], na.rm=T)
			}
	}
catchmentAreas@data = cbind(catchmentAreas@data, hospitalisations, hospitalisationsMR)
catchmentAreas@data[,paste0("hospitalisation_",as.character(selected_dates))] = hospitalisations
catchmentAreas@data[,paste0("hospitalisationMR_",as.character(selected_dates))] = hospitalisationsMR
catchmentAreas@data[,"hospitalisation_2020-06-01_2020-08-31"] = catchmentAreas@data[,"hospitalisation_2020-08-31"]-catchmentAreas@data[,"hospitalisation_2020-05-31"]
catchmentAreas@data[,"hospitalisation_2020-09-01_2020-11-30"] = catchmentAreas@data[,"hospitalisation_2020-11-30"]-catchmentAreas@data[,"hospitalisation_2020-08-31"]
catchmentAreas@data[,"hosp_10^5_habitants_2020-05-31"] = (catchmentAreas@data[,"hospitalisation_2020-05-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitants_2020-08-31"] = (catchmentAreas@data[,"hospitalisation_2020-08-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitants_2020-11-30"] = (catchmentAreas@data[,"hospitalisation_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitants_2020-03-01_2020-05-31"] = (catchmentAreas@data[,"hospitalisation_2020-05-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitants_2020-06-01_2020-08-31"] = (catchmentAreas@data[,"hospitalisation_2020-06-01_2020-08-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitants_2020-09-01_2020-11-30"] = (catchmentAreas@data[,"hospitalisation_2020-09-01_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hospitalisationMR_2020-06-01_2020-08-31"] = catchmentAreas@data[,"hospitalisationMR_2020-08-31"]-catchmentAreas@data[,"hospitalisationMR_2020-05-31"]
catchmentAreas@data[,"hospitalisationMR_2020-09-01_2020-11-30"] = catchmentAreas@data[,"hospitalisationMR_2020-11-30"]-catchmentAreas@data[,"hospitalisationMR_2020-08-31"]
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-05-31"] = (catchmentAreas@data[,"hospitalisationMR_2020-05-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-08-31"] = (catchmentAreas@data[,"hospitalisationMR_2020-08-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-11-30"] = (catchmentAreas@data[,"hospitalisationMR_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-03-01_2020-05-31"] = (catchmentAreas@data[,"hospitalisationMR_2020-05-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-06-01_2020-08-31"] = (catchmentAreas@data[,"hospitalisationMR_2020-06-01_2020-08-31"]/catchmentAreas@data[,"population"])*(10^5)
catchmentAreas@data[,"hosp_10^5_habitantsMR_2020-09-01_2020-11-30"] = (catchmentAreas@data[,"hospitalisationMR_2020-09-01_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)

if (writingFiles) write.csv(catchmentAreas@data, "Catchment_areas_1.csv", quote=F, row.names=F)

if (showingPlots)
	{
		colourScales = list()
		colourScales[[1]] = c(colorRampPalette(brewer.pal(9,"Reds"))(121)[1:101])
		colourScales[[2]] = c(colorRampPalette(brewer.pal(9,"Reds"))(121)[1:101])
		colourScales[[3]] = c(colorRampPalette(brewer.pal(9,"Reds"))(121)[1:101])
		colourScales[[4]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(151)[1:101])
		colourScales[[5]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(151)[1:101])
		colourScales[[6]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
		colourScales[[7]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
		colourScales[[8]] = c(colorRampPalette(brewer.pal(9,"Greys"))(151)[21:121])
		colourScales[[9]] = c(colorRampPalette(brewer.pal(9,"Greens"))(151)[1:101])
		colourScales[[10]] = c(colorRampPalette(brewer.pal(9,"Oranges"))(151)[1:101])
		colourScales[[11]] = c(colorRampPalette(brewer.pal(9,"Blues"))(151)[1:101])
		colourScales[[12]] = c(colorRampPalette(brewer.pal(9,"Purples"))(151)[1:101])
		variables = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30","pm25",
					  "popDensityLog","ratioBedsInMRsPopulation","medianAge","propUrbanArea","sectorP","sectorS","sectorT","medianIncome")
		variableNames = c("CNHP (01/03-31/05/2020)","CNHP (01/09-30/11/2020)","CNHP (01/03-30/11/2020)","Ratio of urban areas","PM 2.5 emission",
						  "Population density (log)","ratio NH beds/population","Median age",">= 65 years (proportion)",
						  "% in primary sector","% in secundary sector","% in tertiary sector","Median declared income")
		dev.new(width=8, height=4.85); par(mfrow=c(3,4), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				values = catchmentAreas@data[,variables[i]]
				minV = min(values); maxV = max(values)
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]
				if (i <= 3)
					{
						maxV = max(catchmentAreas@data[,variables[1:3]])
						maxV = 2500; values[values[]>maxV] = maxV
					}
				legendRast = raster(as.matrix(c(minV,maxV)))
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(catchmentAreas, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames[i], cex=0.54, col="gray30", at=92000, line=-10.5)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.125),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.65, lwd=0,
					 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.07,0)))
			}
		variables = c("hospitalisation_2020-05-31","hospitalisation_2020-06-01_2020-08-31","hospitalisation_2020-09-01_2020-11-30","population",
					  "sectorP","sectorS","sectorT","medianIncome","moreThan65","ratioBedsInMRsPopulation","propUrbanArea","pm25")
		variableNames = c("Cumul. new hosp. 01/03-31/05","Cumul. new hosp. 01/06-31/08","Cumul. new hosp. 01/09-30/11","Population",
						  "% in primary sector","% in secundary sector","% in tertiary sector","Median declared income",
						  ">= 65 years (proportion)","ratio NH beds/population","Ratio of urban areas","PM 2.5 emission")
		dev.new(width=8, height=4.85); par(mfrow=c(3,4), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				values = catchmentAreas@data[,variables[i]]
				minV = min(values); maxV = max(values)
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]
				if (i <= 3)
					{
						maxV = max(catchmentAreas@data[,variables[1:3]]); values[values[]>maxV] = maxV
					}
				legendRast = raster(as.matrix(c(minV,maxV)))
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				if (i <= 4) cols = "gray90"
				plot(catchmentAreas, border="gray30", col=cols, lwd=0.1)
				if (i <= 4)
					{
						if (i <= 3)
							{
								mtext(variableNames[i], cex=0.54, col="gray30", at=110000, line=-12)
								points(catchmentAreas@data[,c("xCentroid","yCentroid")], pch=16, cex=(values/maxV)*5.7, col=rgb(255,0,0,0.3*255,maxColorValue=255))
								points(catchmentAreas@data[,c("xCentroid","yCentroid")], pch=1, cex=(values/maxV)*5.7, col=rgb(255,0,0,255,maxColorValue=255), lwd=0.5)
								legendPoints = rbind(cbind(48000,65000),cbind(87000,56000),cbind(107000,47000))
								points(legendPoints, pch=16, cex=(c(700,400,100)/maxV)*5.7, col=rgb(1,0,0,0.3,1))
								points(legendPoints, pch=1, cex=(c(700,400,100)/maxV)*5.7, col=rgb(1,0,0,1,1), lwd=0.5)
								mtext("700", cex=0.54, col="gray30", at=48000, line=-8.45)
								mtext("400", cex=0.54, col="gray30", at=87000, line=-9.4)
								mtext("100", cex=0.54, col="gray30", at=115000, line=-10.35)
							}
						if (i == 4)
							{
								mtext(variableNames[i], cex=0.54, col="gray30", at=110000, line=-12)
								points(catchmentAreas@data[,c("xCentroid","yCentroid")], pch=16, cex=(values/max(values))*5.7, col=rgb(54,160,84,0.3*255,maxColorValue=255))
								points(catchmentAreas@data[,c("xCentroid","yCentroid")], pch=1, cex=(values/max(values))*5.7, col=rgb(54,160,84,255,maxColorValue=255), lwd=0.5)
								legendPoints = rbind(cbind(45000,57000),cbind(82000,52000),cbind(107000,47000))
								points(legendPoints, pch=16, cex=(c(300000,200000,100000)/max(values))*5.7, col=rgb(54,160,84,0.3*255,maxColorValue=255))
								points(legendPoints, pch=1, cex=(c(300000,200000,100000)/max(values))*5.7, col=rgb(54,160,84,255,maxColorValue=255), lwd=0.5)
								mtext("3*10^5", cex=0.54, col="gray30", at=48000, line=-9.2)
								mtext("2*10^5", cex=0.54, col="gray30", at=82000, line=-9.7)
								mtext("1*10^5", cex=0.54, col="gray30", at=115000, line=-10.35)
							}
					}	else	{
						mtext(variableNames[i], cex=0.54, col="gray30", at=92000, line=-10.5)
						plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.125),
			 				 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.65, lwd=0,
							 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.07,0)))
					}
			}
	}

	# 4.5. Performing and plotting the first axes of an exploratory PCA

if (showingPlots)
	{
		values = catchmentAreas@data[,"hosp_10^5_habitants_2020-11-30"]
		dev.new(width=9, height=2.5); par(mfrow=c(1,3), mar=c(2,3,1,1), lwd=0.2, col="gray30")
		x = log(catchmentAreas@data[,"hosp_10^5_habitants_2020-03-01_2020-05-31"])
		y = log(catchmentAreas@data[,"hosp_10^5_habitants_2020-09-01_2020-11-30"])
		lr = lm(as.formula(paste0("y ~ x"))); R2 = summary(lr)$r.squared
		plot(x, y, col=NA, cex=0.2, pch=16, ann=F, axes=F, xlim=c(1.0,8.0), ylim=c(3.0,8.0))
		abline(lr, lwd=1.5, col="gray80", lty=2); print(c(round(R2,2),cor(x,y,method="spearman")))
		points(x, y, cex=(values/max(values))*8, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
		points(x, y, cex=(values/max(values))*8, pch=1, lwd=0.5, col=rgb(1,0,0,1,1))
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.05,0), at=seq(0,9,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.25,0), at=seq(2,9,1))
		title(xlab="CNHP (01/03-31/05/20, log)", cex.lab=0.8, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="CNHP (01/09-30/11/20, log)", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")		
		legendPoints = rbind(cbind(2,7.73),cbind(2.57,7.65),cbind(3.5,7.5))
		points(legendPoints, pch=16, cex=(c(500,1000,2000)/max(values))*8, col=rgb(0,0,1,0.3,1))
		points(legendPoints, pch=1, cex=(c(500,1000,2000)/max(values))*8, col=rgb(0,0,1,1,1), lwd=0.5)
		x = log(catchmentAreas@data[,"hosp_10^5_habitants_2020-09-01_2020-11-30"])-log(catchmentAreas@data[,"hosp_10^5_habitants_2020-03-01_2020-05-31"])
		y = log(catchmentAreas@data[,"hosp_10^5_habitants_2020-09-01_2020-11-30"])
		lr = lm(as.formula(paste0("y ~ x"))); R2 = summary(lr)$r.squared
		plot(x, y, col=NA, cex=0.2, pch=16, ann=F, axes=F, xlim=c(-1.8,4.2), ylim=c(3.0,7.8))
		abline(lr, lwd=1.5, col="gray80", lty=2); print(c(round(R2,2),cor(x,y,method="spearman")))
		points(x, y, cex=(values/max(values))*8, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
		points(x, y, cex=(values/max(values))*8, pch=1, lwd=0.5, col=rgb(1,0,0,1,1))
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-2,5,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.25,0), at=seq(1,8,1))
		title(xlab="CNHP (01/09-30/11/20, log) - CNHP (01/03-31/05/20, log)", cex.lab=0.8, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="CNHP (01/09-30/11/20, log)", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")		
		x = catchmentAreas@data[,"ratioBedsInMRsPopulation"]
		y = log(catchmentAreas@data[,"hosp_10^5_habitants_2020-11-30"])
		lr = lm(as.formula(paste0("y ~ x"))); R2 = summary(lr)$r.squared
		plot(x, y, col=NA, cex=0.2, pch=16, ann=F, axes=F, xlim=c(0,0.03), ylim=c(4.0,8.5))
		abline(lr, lwd=1.5, col="gray80", lty=2); print(c(round(R2,2),cor(x,y,method="spearman")))
		points(x, y, cex=(values/max(values))*8, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
		points(x, y, cex=(values/max(values))*8, pch=1, lwd=0.5, col=rgb(1,0,0,1,1))
		axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-0.005,0.035,0.005))
		axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.25,0), at=seq(3,9,1))
		title(xlab="ratio NH beds/population", cex.lab=0.8, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="CNHP (30/11/2020, log)", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")		
	}
if (showingPlots)
	{
		values = (catchmentAreas@data[,"hospitalisation_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)
		variableNames = c("popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation",
						  "medianIncome","sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")
		df = catchmentAreas@data[,variableNames]
		colnames(df) = c("Pop. density","median age","> 65 yrs","MR beds/population","median income",
						 "sector P","sector S","sector T","PM 10","PM 2.5","prop. urban area")
		pca = dudi.pca(df, scannf=F, nf=length(variableNames)); lis = pca$li[,1:2]; cos = pca$co
		dev.new(width=6, height=6); par(mar=c(3,3,1.5,1.5), lwd=0.2, col="gray30")
		plot(lis, col=NA, cex=0.2, pch=16, ann=F, axes=F, xlim=c(-4.3,8), ylim=c(-5,3.5))
		points(lis, cex=(values/max(values))*8, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
		points(lis, cex=(values/max(values))*8, pch=1, lwd=0.5, col=rgb(1,0,0,1,1))
		s.corcircle(2*cos, xax=1, yax=2, box=F, sub="", csub=0.7, clabel=0.7, possub="topleft", grid=F, cgrid=1, full=F, add.plot=T)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-12,10,2))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-7,9,1))
		title(xlab=paste0("PCA axis 1 (",round((pca$eig[1]/(sum(pca$eig))*100),1),"%)"), cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab=paste0("PCA axis 2 (",round((pca$eig[2]/(sum(pca$eig))*100),1),"%)"), cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		legendPoints = rbind(cbind(6,-4.73),cbind(6.57,-4.65),cbind(7.5,-4.5))
		points(legendPoints, pch=16, cex=(c(500,1000,2000)/max(values))*8, col=rgb(1,0,0,0.3,1))
		points(legendPoints, pch=1, cex=(c(500,1000,2000)/max(values))*8, col=rgb(1,0,0,1,1), lwd=0.5)
		mtext("500", cex=0.7, col="gray30", at=5.8, line=-23.58)
		mtext("1000", cex=0.7, col="gray30", at=6.57, line=-23.10)
		mtext("2000", cex=0.7, col="gray30", at=7.5, line=-22.25)
		values = (catchmentAreas@data[,"hospitalisation_2020-11-30"]/catchmentAreas@data[,"population"])*(10^5)
		variableNames = c("hosp_10^5_habitants_2020-03-01_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30",
						  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation",
						  "medianIncome","sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")
		df = catchmentAreas@data[,variableNames]
		colnames(df) = c("CNHP 01/03-31/05","CNHP 01/09-30/11",
						 "Pop. density","median age","> 65 yrs","MR beds/population","median income",
						 "sector P","sector S","sector T","PM 10","PM 2.5","prop. urban area")
		pca = dudi.pca(df, scannf=F, nf=length(variableNames)); lis = pca$li[,2:1]; cos = pca$co
		dev.new(width=9, height=4.5); par(mar=c(3,3,1.5,1.5), lwd=0.2, col="gray30")
		plot(lis, col=NA, cex=0.2, pch=16, ann=F, axes=F, xlim=c(-4,5), ylim=c(-4.5,4))
		points(lis, cex=(values/max(values))*8, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
		points(lis, cex=(values/max(values))*8, pch=1, lwd=0.5, col=rgb(1,0,0,1,1))
		s.corcircle(2*cos, xax=2, yax=1, box=F, sub="", csub=0.7, clabel=0.7, possub="topleft", grid=F, cgrid=1, full=F, add.plot=T)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,8,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-8,8,2))
		title(xlab=paste0("PCA axis 2 (",round((pca$eig[2]/(sum(pca$eig))*100),1),"%)"), cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab=paste0("PCA axis 1 (",round((pca$eig[1]/(sum(pca$eig))*100),1),"%)"), cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
	}

	# 4.6. Classic correlation analyses (correlogram)

variableNames = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30",
				  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation",
				  "medianIncome","sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")
df = catchmentAreas@data[,variableNames]; colnames(df) = gsub("-","",colnames(df))
colnames(df) = c("CNHP - 01/03-31/05/20","CNHP - 01/09-30/11/20","CNHP - 01/03-30/11/20",
				 "population density","median age","prop. >65 years old","ratio NH beds/population",
				 "median income","% in primary sector","% in secundary sector","% in tertiary sector",
				 "PM 10 emission","PM 2.5 emission","prop. urban areas")
correlations_spearman = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
correlations_ranking = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
pValues_spearman = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
pValues_ranking = matrix(nrow=dim(df)[2], ncol=dim(df)[2])
for (i in 2:dim(correlations_spearman)[2])
	{
		for (j in 1:(i-1))
			{
				vS = df[,c(i,j)]
				vS[is.infinite(vS[,1]),1] = NA
				vS[is.infinite(vS[,2]),1] = NA
				vS = vS[which((!is.na(vS[,1]))&(!is.na(vS[,2]))),]	
				test = cor.test(vS[,1], vS[,2], method="spearman")
				correlations_spearman[j,i] = test$estimate
				pValues_spearman[j,i] = test$p.value
				ranked_1 = vS[order(vS[,1]),1]
				ranked_2 = vS[order(vS[,2]),2]
				rS = vS; rS[] = NA
				for (k in 1:dim(rS)[1])
					{
						rS[k,1] = which(ranked_1==vS[k,1])
						rS[k,2] = which(ranked_2==vS[k,2])
					}
				test = cor.test(rS[,1], rS[,2], method="pearson")
				correlations_ranking[j,i] = test$estimate
				pValues_ranking[j,i] = test$p.value
			}
	}
if (showingPlots)
	{
		texts = correlations_ranking
		for (i in 2:dim(texts)[1])
			{
				for (j in 1:(i-1))
					{
						if (pValues[j,i] < 0.05)
							{
								texts[j,i] = as.character(round(correlations_spearman[j,i],1))
							}	else	{
								texts[j,i] = NA
							}
					}
			}
		cols = colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:111]
		minV = min(correlations_ranking, na.rm=T); maxV = max(correlations_ranking, na.rm=T)
		index1 = (50+(minV*50))+1; index2 = (50+(maxV*50))+1; cols = cols[index1:index2]
		heatmap.2(correlations_ranking, cellnote=texts, notecex=0.6, notecol="gray30", main=NULL, density.info="none",
				  trace="none", margins=c(12,9), col=cols, Rowv=NULL, Colv="NA", dendrogram="none", key=F,
				  labRow=colnames(df), labCol=colnames(df), cexRow=0.8, cexCol=0.8, offsetRow=0.0,
				  colRow=rep("gray30",length(variableNames)), colCol=rep("gray30",length(variableNames)))
		cols = colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:111]; par(lwd=0.2)
		plot(raster(as.matrix(c(-1,1))), legend.only=T, col=cols, legend.lwd=0.2, legend.width=0.2, legend.shrink=0.3, smallplot=c(0.32,0.78,0.80,0.81),
			 	  alpha=1, horizontal=T, lwd=0.2, legend.args=list(text="", cex=0.7, line=0.2, col="gray30"), axis.args=list(cex.axis=0.65, lwd=0,
				  lwd.tick=0.2, tck=-0.6, col.axis="gray30", col.tick="gray30", line=0, mgp=c(0,0.07,0)))
		texts = correlations_spearman
		for (i in 2:dim(texts)[1])
			{
				for (j in 1:(i-1))
					{
						if (pValues[j,i] < 0.05)
							{
								texts[j,i] = as.character(round(correlations_spearman[j,i],1))
							}	else	{
								texts[j,i] = NA
							}
					}
			}
		minV = min(correlations_spearman, na.rm=T); maxV = max(correlations_spearman, na.rm=T)
		index1 = (50+(minV*50))+1; index2 = (50+(maxV*50))+1; cols = cols[index1:index2]
		heatmap.2(correlations_spearman, cellnote=texts, notecex=0.6, notecol="gray30", main=NULL, density.info="none",
				  trace="none", margins=c(12,9), col=cols, Rowv=NULL, Colv="NA", dendrogram="none", key=F,
				  labRow=colnames(df), labCol=colnames(df), cexRow=0.8, cexCol=0.8, offsetRow=0.0,
				  colRow=rep("gray30",length(variableNames)), colCol=rep("gray30",length(variableNames)))
		cols = colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:111]; par(lwd=0.2)
		plot(raster(as.matrix(c(-1,1))), legend.only=T, col=cols, legend.lwd=0.2, legend.width=0.2, legend.shrink=0.3, smallplot=c(0.32,0.78,0.80,0.81),
			 	  alpha=1, horizontal=T, lwd=0.2, legend.args=list(text="", cex=0.7, line=0.2, col="gray30"), axis.args=list(cex.axis=0.65, lwd=0,
				  lwd.tick=0.2, tck=-0.6, col.axis="gray30", col.tick="gray30", line=0, mgp=c(0,0.07,0)))
	}

	# 4.7. Assessing spatial autocorrelation with the Moran's I test

variableNames = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30",
				  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation","medianIncome",
				  "sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea","xCentroid","yCentroid")
df = catchmentAreas@data[,variableNames]; colnames(df) = gsub("-","",colnames(df))
responseVariables = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30")
responseVariables = gsub("-","",responseVariables)
for (i in 1:length(responseVariables))
	{
		values = df[,responseVariables[i]]
		indices = which((!is.na(values))&(!is.infinite(values)))
		geoDists = as.matrix(dist(df[indices,c("xCentroid","yCentroid")]))
		weights = 1/geoDists; diag(weights) = 0
		print(Moran.I(values[indices],weights))
			# H0: no spatial autocorrelation
	}

	# 4.8. Univariate (LR) followed by multivariate regression (GLM) analyses

normal_analyses = FALSE; analyses_on_MR_cases = FALSE; analyses_without_ouliers_BXL = TRUE
normal_analyses = FALSE; analyses_on_MR_cases = TRUE; analyses_without_ouliers_BXL = FALSE
normal_analyses = TRUE; analyses_on_MR_cases = FALSE; analyses_without_ouliers_BXL = FALSE
if (analyses_without_ouliers_BXL == TRUE) outliersToExclude = c("689","077-150","110","723","076-079","087")

predictors = c("popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation","medianIncome",
			   "sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")	
table1 = matrix(nrow=11, ncol=9); row.names(table1) = predictors # all hospitalisation cases
tableS1 = matrix(nrow=11, ncol=9); row.names(tableS1) = predictors # only MR hospitalisations
tableS2 = matrix(nrow=11, ncol=9); row.names(tableS2) = predictors # without the outlier and Brussels

if (analyses_on_MR_cases != TRUE)
	{
		variableNames = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30",
						  "hosp_10^5_habitants_2020-11-30",predictors,"xCentroid","yCentroid")
		responseVariables = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30")
	}	else		{
		variableNames = c("hosp_10^5_habitantsMR_2020-05-31","hosp_10^5_habitantsMR_2020-09-01_2020-11-30",
						  "hosp_10^5_habitantsMR_2020-11-30",predictors,"xCentroid","yCentroid")
		responseVariables = c("hosp_10^5_habitantsMR_2020-05-31","hosp_10^5_habitantsMR_2020-09-01_2020-11-30","hosp_10^5_habitantsMR_2020-11-30")
	}
responseVariables = gsub("\\^","e",gsub("-","",responseVariables)); selectedVariables = list()
for (i in 1:length(responseVariables))
	{
		df = catchmentAreas@data[,variableNames]; buffer = c()
		colnames(df) = gsub("\\^","e",gsub("-","",colnames(df)))
		if (analyses_without_ouliers_BXL == TRUE)
			{
				df = df[which(!catchmentAreas@data[,"X_ID"]%in%outliersToExclude),]
			}
		for (j in 1:length(predictors))
			{
				values = df[,responseVariables[i]]; tmp = df
				tmp = tmp[which((!is.na(values))&(!is.infinite(values))),]
				formula = paste0(responseVariables[i]," ~ ",predictors[j])
				lr = lm(formula, data=tmp)
				R2 = summary(lr)$r.squared
				f = summary(lr)$fstatistic
				pValue = pf(f[1],f[2],f[3],lower.tail=F)
				if (pValue < 0.05)
					{
						buffer = c(buffer, predictors[j])
						if (normal_analyses == TRUE) table1[j,((i-1)*3)+1] = paste0(round(R2,2),"*")
						if (analyses_on_MR_cases == TRUE) tableS1[j,((i-1)*3)+1] = paste0(round(R2,2),"*") 
						if (analyses_without_ouliers_BXL == TRUE) tableS2[j,((i-1)*3)+1] = paste0(round(R2,2),"*")
					}	else	{
						if (normal_analyses == TRUE) table1[j,((i-1)*3)+1] = round(R2,2)
						if (analyses_on_MR_cases == TRUE) tableS1[j,((i-1)*3)+1] = round(R2,2) 
						if (analyses_without_ouliers_BXL == TRUE) tableS2[j,((i-1)*3)+1] = round(R2,2)
					}
			}
		selectedVariables[[i]] = buffer
		if (is.null(buffer)) selectedVariables[[i]] = NA
	}
for (i in 1:length(responseVariables))
	{
		df = catchmentAreas@data[,variableNames]
		colnames(df) = gsub("\\^","e",gsub("-","",colnames(df)))
		for (j in 1:dim(df)[2]) df[,j] = (df[,j]-min(df[,j]))/(max(df[,j])-min(df[,j]))
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
				lrm = lm(formula, data=tmp); print(round(summary(lrm)$r.squared,2)); cat("\n")
				for (j in 2:length(row.names(summary(lrm)$coefficients)))
					{
						pValue = summary(lrm)$coefficients[j,4]
						beta = summary(lrm)$coefficients[j,1]
						if (pValue < 0.05)
							{
								if (normal_analyses == TRUE) table1[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = paste0(round(beta,2),"*")
								if (analyses_on_MR_cases == TRUE) tableS1[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = paste0(round(beta,2),"*")
								if (analyses_without_ouliers_BXL == TRUE) tableS2[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = paste0(round(beta,2),"*")
							}	else		{
								if (normal_analyses == TRUE) table1[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = paste0(round(beta,2),"*")
								if (analyses_on_MR_cases == TRUE) tableS1[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = paste0(round(beta,2),"*")
								if (analyses_without_ouliers_BXL == TRUE) tableS2[row.names(summary(lrm)$coefficients)[j],((i-1)*3)+2] = round(beta,2)
							}
					}
			}
	}

	# 4.9. Multivariate GAM (generalised additive model) analyses

variableNames = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30",
				  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation","medianIncome",
				  "sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea","xCentroid","yCentroid")
responseVariables = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30")
responseVariables = gsub("\\^","e",gsub("-","",responseVariables)); selectedVariables = list(); gams = list(); zTransformations = FALSE
for (i in 1:length(responseVariables))
	{
		values = df[,responseVariables[i]]; tmp = df
		tmp = tmp[which((!is.na(values))&(!is.infinite(values))),]
		colnames(tmp) = gsub(responseVariables[i],"responseVariable",colnames(tmp))
		if (zTransformations == TRUE)
			{
				for (j in 1:dim(dfs[[i]])[2]) df[,j] = zTransformation(df[,j])
			}
		gam = gam(responseVariable ~ s(ratioBedsInMRsPopulation) + s(sectorP) + s(propUrbanArea), data=tmp, method="REML")
		print(summary(gam)); cat("\n\n"); gams[[i]] = gam
		if (showingPlots) dev.new(); plot(gam, pages=1)
	}
if (showingPlots)
	{
		gam = gams[[3]]; responseCurves = list()
		curves = plot(gam, pages=1); dev.off()
		selectedVariables = c("ratioBedsInMRsPopulation")
		variableNames = c("ratio NH beds/population")
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
				tmp = df[,c("ratioBedsInMRsPopulation","sectorP","propUrbanArea")]
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

	# 4.10. Multivariate analyses with the boosted regression trees approach

normal_analyses = FALSE; analyses_on_MR_cases = FALSE; analyses_without_ouliers_BXL = TRUE
normal_analyses = FALSE; analyses_on_MR_cases = TRUE; analyses_without_ouliers_BXL = FALSE
normal_analyses = TRUE; analyses_on_MR_cases = FALSE; analyses_without_ouliers_BXL = FALSE
if (analyses_on_MR_cases == FALSE)
	{
		variableNames = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30",
						  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation","medianIncome",
						  "sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")
		responseVariables = c("hosp_10^5_habitants_2020-05-31","hosp_10^5_habitants_2020-09-01_2020-11-30","hosp_10^5_habitants_2020-11-30")
	}	else		{
		variableNames = c("hosp_10^5_habitantsMR_2020-05-31","hosp_10^5_habitantsMR_2020-09-01_2020-11-30","hosp_10^5_habitantsMR_2020-11-30",
						  "popDensity","medianAge","moreThan65","ratioBedsInMRsPopulation","medianIncome",
						  "sectorP","sectorS","sectorT","pm10","pm25","propUrbanArea")
		responseVariables = c("hosp_10^5_habitantsMR_2020-05-31","hosp_10^5_habitantsMR_2020-09-01_2020-11-30","hosp_10^5_habitantsMR_2020-11-30")
	}
responseVariables = gsub("\\^","e",gsub("-","",responseVariables)); nberOfReplicates = 10
if (!file.exists("All_the_BRT_models.rds"))
	{
		brts_list = list()
		for (i in 1:length(responseVariables))
			{	
				df = catchmentAreas@data[,variableNames]; brts = list()
				if (analyses_without_ouliers_BXL == TRUE)
					{
						df = df[which(!catchmentAreas@data[,"X_ID"]%in%outliersToExclude),]
					}
				colnames(df) = gsub("\\^","e",gsub("-","",colnames(df)))
				sub = df[,c(responseVariables[i],colnames(df)[!colnames(df)%in%responseVariables])]
				sub[,responseVariables[i]] = round(sub[,responseVariables[i]]) # for family="possoin"
				gbm.x = colnames(df)[!colnames(df)%in%responseVariables]; gbm.y = responseVariables[i]
				offset = NULL; fold.vector = NULL; tree.complexity = 5; learning.rate = 0.005; bag.fraction = 0.75
				site.weights = rep(1,dim(sub)[1]); var.monotone = rep(0,length(gbm.x)); n.folds = 5; prev.stratify = TRUE
				family = "poisson"; step.size = 5; max.trees = 1000; tolerance.method = "auto"
				tolerance = 0.001; plot.main = TRUE; plot.folds = FALSE; verbose = TRUE; silent = FALSE
				keep.fold.models = FALSE; keep.fold.vector = FALSE; keep.fold.fit = FALSE; showingFoldsPlot = FALSE
				for (j in 1:nberOfReplicates)
					{
						n.trees = 10
						brts[[j]] = gbm.step(sub, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									 		 var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance,
											 plot.main, plot.folds, verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
					}
				brts_list[[i]] = brts
			}
		df = catchmentAreas@data[,variableNames]; brts = list()
		colnames(df) = gsub("\\^","e",gsub("-","",colnames(df)))
		sub = df[,c(responseVariables[2],responseVariables[1],colnames(df)[!colnames(df)%in%responseVariables])]
		sub[,responseVariables[2]] = round(sub[,responseVariables[2]]) # for family="possoin"
		sub[,responseVariables[1]] = round(sub[,responseVariables[1]]) # for family="possoin"
		gbm.x = c(responseVariables[1],colnames(df)[!colnames(df)%in%responseVariables]); gbm.y = responseVariables[2]
		offset = NULL; fold.vector = NULL; tree.complexity = 5; learning.rate = 0.005; bag.fraction = 0.75
		site.weights = rep(1,dim(sub)[1]); var.monotone = rep(0,length(gbm.x)); n.folds = 5; prev.stratify = TRUE
		family = "poisson"; step.size = 5; max.trees = 1000; tolerance.method = "auto"
		tolerance = 0.001; plot.main = TRUE; plot.folds = FALSE; verbose = TRUE; silent = FALSE
		keep.fold.models = FALSE; keep.fold.vector = FALSE; keep.fold.fit = FALSE; showingFoldsPlot = FALSE
		for (j in 1:nberOfReplicates)
			{
				n.trees = 10
				brts[[j]] = gbm.step(sub, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							 		 var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance,
									 plot.main, plot.folds, verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
			}
		brts_list[[4]] = brts
		if (normal_analyses == TRUE) saveRDS(brts_list, "All_the_BRT_models.rds")
	}	else		{
		if (normal_analyses == TRUE) brts_list = readRDS("All_the_BRT_models.rds")
	}
if (analyses_on_MR_cases != TRUE) responseVariables_mod = c(responseVariables, "hosp_10e5_habitants_20200901_20201130")
if (analyses_on_MR_cases == TRUE) responseVariables_mod = c(responseVariables, "hosp_10e5_habitantsMR_20200901_20201130")
meanCorrelations = matrix(nrow=length(responseVariables_mod), ncol=length(responseVariables))
for (i in 1:length(brts_list))
	{
		for (j in 1:(length(responseVariables_mod)-1))
			{
				df1 = catchmentAreas@data[,variableNames]; correlations = c()
				if (analyses_without_ouliers_BXL == TRUE)
					{
						df1 = df1[which(!catchmentAreas@data[,"X_ID"]%in%outliersToExclude),]
					}
				colnames(df1) = gsub("\\^","e",gsub("-","",colnames(df1)))
				if (i <= 3)
					{
						sub = df1[,c(responseVariables_mod[i],colnames(df1)[!colnames(df1)%in%responseVariables_mod])]
						sub[,responseVariables_mod[i]] = round(sub[,responseVariables_mod[i]]); df2 = sub
					}	else		{
						sub = df1[,c(responseVariables_mod[i],responseVariables_mod[1],colnames(df1)[!colnames(df1)%in%responseVariables_mod])]
						sub[,responseVariables_mod[i]] = round(sub[,responseVariables_mod[i]])
						sub[,responseVariables_mod[1]] = round(sub[,responseVariables_mod[1]]); df2 = sub
					}
				for (k in 1:length(brts_list[[i]]))
					{
						n.trees = brts_list[[i]][[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brts_list[[i]][[k]], newdata=df2, n.trees, type, single.tree)
						correlations = c(correlations, cor(df1[,responseVariables_mod[j]],prediction,method="spearman"))
					}
				meanCorrelations[i,j] = paste0(round(mean(correlations),2)," [",round(min(correlations),2),"-",round(max(correlations),2),"]")
			}
	}
row.names(meanCorrelations) = responseVariables_mod; colnames(meanCorrelations) = responseVariables
if (writingFiles) write.table(meanCorrelations, "BRT_correlations.csv", quote=F, sep=",")
for (i in 1:length(responseVariables))
	{
		variables = summary(brts_list[[i]][[1]])$var
		tab1 = matrix(nrow=length(variables), ncol=length(brts_list[[i]]))
		for (j in 1:length(brts_list[[i]]))
			{
				for (k in 1:length(variables)) tab1[k,j] = summary(brts_list[[i]][[j]])[variables[k],"rel.inf"]
			}
		tab2 = matrix(nrow=length(variables), ncol=1); row.names(tab2) = variables
		tab3 = matrix(nrow=length(variables), ncol=1); row.names(tab3) = variables
		for (j in 1:length(variables)) tab2[j,1] = paste0("[",round(min(tab1[j,]),1),"-",round(max(tab1[j,]),1),"]")
		for (j in 1:length(variables)) tab3[j,1] = round(mean(tab1[j,]),1)
		if (normal_analyses == TRUE) table1[variables,((i-1)*3)+3] = tab3
		if (analyses_on_MR_cases == TRUE) tableS1[variables,((i-1)*3)+3] = tab3
		if (analyses_without_ouliers_BXL == TRUE) tableS2[variables,((i-1)*3)+3] = tab3
	}
if (writingFiles)
	{
		if (normal_analyses == TRUE) write.table(table1, "Table_1_TEMP.txt", sep="\t", quote=F, row.names=F, col.names=F)
		if (analyses_on_MR_cases == TRUE) write.table(tableS1, "Table_S1_TEMP.txt", sep="\t", quote=F, row.names=F, col.names=F)
		if (analyses_without_ouliers_BXL == TRUE) write.table(tableS2, "Table_S2_TEMP.txt", sep="\t", quote=F, row.names=F, col.names=F)
	}
for (i in 1:length(responseVariables))
	{
		df1 = catchmentAreas@data[,variableNames]
		if (analyses_without_ouliers_BXL == TRUE)
			{
				df1 = df1[which(!catchmentAreas@data[,"X_ID"]%in%outliersToExclude),]
			}
		colnames(df1) = gsub("\\^","e",gsub("-","",colnames(df1)))
		sub = df1[,c(responseVariables[i],colnames(df1)[!colnames(df1)%in%responseVariables])]
		sub[,responseVariables[i]] = round(sub[,responseVariables[i]]); df2s = list()
		variableNames_to_print = c("","population density","median age","proportion of >65 years old","ratio NH beds/population","median income",
				 				   "% of workers in primary sector","% of workers in secundary sector","% of workers in tertiary sector",
				 				   "PM 10 emission","PM 2.5 emission","proportion of urban areas")
		pdf(paste0("BRT_RC_",i,".pdf"), width=7.5, height=4)
		par(mfrow=c(3,4), oma=c(1,1,1,1), mar=c(2,1.3,0.5,0.5), lwd=0.2, col="gray30")
		for (j in 2:dim(sub)[2])
			{
				valuesInterval = (max(sub[,j])-min(sub[,j]))/100
				df2 = data.frame(matrix(nrow=length(seq(min(sub[,j]),max(sub[,j]),valuesInterval)),ncol=dim(sub)[2]-1))
				colnames(df2) = colnames(sub)[2:dim(sub)[2]]
				for (k in 2:dim(sub)[2])
					{
						valuesInterval = (max(sub[,k])-min(sub[,k]))/100
						if (j == k) df2[,colnames(sub)[k]] = seq(min(sub[,k]),max(sub[,k]),valuesInterval)
						if (j != k) df2[,colnames(sub)[k]] = rep(median(sub[,k]),dim(df2)[1])
					}
				predictions = list()
				for (k in 1:length(brts_list[[i]]))
					{
						n.trees = brts_list[[i]][[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brts_list[[i]][[k]], newdata=df2, n.trees, type, single.tree)
						if (k == 1)
							{
								minX = min(df2[,colnames(sub)[j]]); maxX = max(df2[,colnames(sub)[j]])
								minY = min(prediction); maxY = max(prediction)
							}	else	{
								if (minX > min(df2[,colnames(sub)[j]])) minX = min(df2[,colnames(sub)[j]])
								if (maxX < max(df2[,colnames(sub)[j]])) maxX = max(df2[,colnames(sub)[j]])
								if (minY > min(prediction)) minY = min(prediction)
								if (maxY < max(prediction)) maxY = max(prediction)
							}
						predictions[[k]] = prediction
					}
				for (k in 1:length(brts_list[[i]]))
					{
						if (k == 1)
							{
								plot(df2[,colnames(sub)[j]],predictions[[k]],col="gray30",ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
							}	else	{
								lines(df2[,colnames(sub)[j]],predictions[[k]],col="gray30",lwd=0.2)
							}
					}
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", col.tick="gray30", mgp=c(0,0.07,0))
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", col.tick="gray30", mgp=c(0,0.30,0))
				title(xlab=variableNames_to_print[j], cex.lab=0.8, mgp=c(0.9,0,0), col.lab="gray30")
				box(lwd=0.2, col="gray30")
			}
		dev.off()
	}

# 5. Analyses of hospital catchment areas (temporal)

days1 = c(paste0("2020-03-0",c(1:9)),paste0("2020-03-",c(10:31)),paste0("2020-04-0",c(1:9)),paste0("2020-04-",c(10:30)),
		  paste0("2020-05-0",c(1:9)),paste0("2020-05-",c(10:31)),paste0("2020-06-0",c(1:9)),paste0("2020-06-",c(10:30)),
		  paste0("2020-07-0",c(1:9)),paste0("2020-07-",c(10:31)),paste0("2020-08-0",c(1:9)),paste0("2020-08-",c(10:31)),
		  paste0("2020-09-0",c(1:9)),paste0("2020-09-",c(10:30)),paste0("2020-10-0",c(1:9)),paste0("2020-10-",c(10:31)),
		  paste0("2020-11-0",c(1:9)),paste0("2020-11-",c(10:30)))

	# 5.1. Cumputing doubling times for hospitalisations and ICU

data = read.csv("Raw_data_Sciensano/Hosp_surge_overview_03122020_1500.csv", head=T, sep=";")
firstDay = ymd("2020-01-30"); days2 = ymd(days1); days3 = as.numeric(days2-firstDay)
data$date = rep(NA, dim(data)[1]); D = 7 # time interval
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
dailyRatioHosCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
dailyRatioICUCases = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
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
						newHospitaliCases[i,j] = temp[index2,"NewPatientsNotReferredHospital"]
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
						dailyRatioHosCases[i,j] = QTD2/QTD1
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
						dailyRatioICUCases[i,j] = QTD2/QTD1
					}
			}
	}
colnames(newHospitaliCases) = paste0("newHospitaliCases_",days1)
colnames(cumulatedHosCases) = paste0("cumulatedHosCases_",days1)
colnames(incidenceHosCases) = paste0("incidenceHosCases_",days1)
colnames(doublingTHosCases) = paste0("doublingTHosCases_",days1)
colnames(dailyRatioHosCases) = paste0("dailyRatioHosCases_",days1)
colnames(cumulatedICUCases) = paste0("cumulatedICUCases_",days1)
colnames(incidenceICUCases) = paste0("incidenceICUCases_",days1)
colnames(doublingTICUCases) = paste0("doublingTICUCases_",days1)
colnames(dailyRatioICUCases) = paste0("dailyRatioICUCases_",days1)
catchmentAreas@data = cbind(catchmentAreas@data, newHospitaliCases)
catchmentAreas@data = cbind(catchmentAreas@data, cumulatedHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, incidenceHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, doublingTHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, dailyRatioHosCases)
catchmentAreas@data = cbind(catchmentAreas@data, cumulatedICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, incidenceICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, doublingTICUCases)
catchmentAreas@data = cbind(catchmentAreas@data, dailyRatioICUCases)

	# 5.2. Extracting mobility data from mobile phone data

		# - out_per_capita = # trips out / # subscriptions
		# - in_per_capita = # trips in / # subscriptions
		# - in_out_per_capita = (# trips in + trips out) / # subscriptions
		# - trips out = journey of >15min outside the zip code
		# - trips in = journey of >15min inside the zip code (for someone that is not living in that zip code)
		# - subscriptions = number of SIM cards (Proximus or Telenet) "living" in the zip code, even for in_per_capita (!)
		# - n.b.: the negative in_per_capita result from not enough inhabitants, because this is not allowed to report counts < 30 
		# 		  (so they put these values on "-1"). These negative values should thus be treated as "NA"

if (!file.exists("Raw_data_Sciensano/proximus_catchment_areas_all.csv"))
	{
		days2 = gsub("-","",days1); days3 = ymd(days1)
		data = read.csv("Raw_data_Sciensano/proximus_in_out_zip.csv")
		proximusIndexCommunesIn = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		proximusIndexCommunesOut = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		proximusIndexCommunesAll = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		row.names(proximusIndexCommunesIn) = communes2@data[,"nouveau_PO"]; colnames(proximusIndexCommunesIn) = days1
		row.names(proximusIndexCommunesOut) = communes2@data[,"nouveau_PO"]; colnames(proximusIndexCommunesOut) = days1
		row.names(proximusIndexCommunesAll) = communes2@data[,"nouveau_PO"]; colnames(proximusIndexCommunesAll) = days1
		for (i in 1:dim(communes2@data)[1])
			{
				for (j in 1:length(days2))
					{
						indices = which((data[,"postalcode"]==communes2@data[i,"nouveau_PO"])&(data[,"date"]==days2[j]))
						if (length(indices) > 0)	
							{
								proximusIndexCommunesIn[i,j] = mean(data[indices,"in_per_capita"])
								proximusIndexCommunesOut[i,j] = mean(data[indices,"out_per_capita"])
								proximusIndexCommunesAll[i,j] = mean(data[indices,"in_out_per_capita"])
							}
					}
			}
		proximusIndexCatchmentsIn = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		proximusIndexCatchmentsOut = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		proximusIndexCatchmentsAll = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(proximusIndexCatchmentsIn) = paste0("proximusIndexIn_",days1)
		colnames(proximusIndexCatchmentsOut) = paste0("proximusIndexOut_",days1)
		colnames(proximusIndexCatchmentsAll) = paste0("proximusIndexAll_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						proximusIndexCatchmentsIn[i,j] = sum(proximusIndexCommunesIn[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
						proximusIndexCatchmentsOut[i,j] = sum(proximusIndexCommunesOut[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
						proximusIndexCatchmentsAll[i,j] = sum(proximusIndexCommunesAll[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(proximusIndexCatchmentsIn, "Raw_data_Sciensano/proximus_catchment_areas_in.csv", quote=F, row.names=F, sep=",")
		write.table(proximusIndexCatchmentsOut, "Raw_data_Sciensano/proximus_catchment_areas_out.csv", quote=F, row.names=F, sep=",")
		write.table(proximusIndexCatchmentsAll, "Raw_data_Sciensano/proximus_catchment_areas_all.csv", quote=F, row.names=F, sep=",")
		data = read.csv("Raw_data_Sciensano/telenet_zip_per_capita.csv")
		telenetIndexCommunesIn = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		telenetIndexCommunesOut = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		telenetIndexCommunesAll = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		row.names(telenetIndexCommunesIn) = communes2@data[,"nouveau_PO"]; colnames(telenetIndexCommunesIn) = days1
		row.names(telenetIndexCommunesOut) = communes2@data[,"nouveau_PO"]; colnames(telenetIndexCommunesOut) = days1
		row.names(telenetIndexCommunesAll) = communes2@data[,"nouveau_PO"]; colnames(telenetIndexCommunesAll) = days1
		for (i in 1:dim(communes2@data)[1])
			{
				for (j in 1:length(days2))
					{
						indices = which((data[,"postalcode"]==communes2@data[i,"nouveau_PO"])&(data[,"date"]==days2[j]))
						if (length(indices) > 0)	
							{
								telenetIndexCommunesIn[i,j] = mean(data[indices,"in_per_capita"])
								telenetIndexCommunesOut[i,j] = mean(data[indices,"out_per_capita"])
								telenetIndexCommunesAll[i,j] = mean(data[indices,"in_out_per_capita"])
							}
					}
			}
		telenetIndexCatchmentsIn = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		telenetIndexCatchmentsOut = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		telenetIndexCatchmentsAll = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(telenetIndexCatchmentsIn) = paste0("telenetIndexIn_",days1)
		colnames(telenetIndexCatchmentsOut) = paste0("telenetIndexOut_",days1)
		colnames(telenetIndexCatchmentsAll) = paste0("telenetIndexAll_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						telenetIndexCatchmentsIn[i,j] = sum(telenetIndexCommunesIn[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
						telenetIndexCatchmentsOut[i,j] = sum(telenetIndexCommunesOut[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
						telenetIndexCatchmentsAll[i,j] = sum(telenetIndexCommunesAll[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(telenetIndexCatchmentsIn, "Raw_data_Sciensano/telenet_catchment_areas_in.csv", quote=F, row.names=F, sep=",")
		write.table(telenetIndexCatchmentsOut, "Raw_data_Sciensano/telenet_catchment_areas_out.csv", quote=F, row.names=F, sep=",")
		write.table(telenetIndexCatchmentsAll, "Raw_data_Sciensano/telenet_catchment_areas_all.csv", quote=F, row.names=F, sep=",")
	}
proximusIndexCatchmentsAll = read.csv("Raw_data_Sciensano/proximus_catchment_areas_all.csv")
colnames(proximusIndexCatchmentsAll) = gsub("\\.","-",colnames(proximusIndexCatchmentsAll))
catchmentAreas@data = cbind(catchmentAreas@data, proximusIndexCatchmentsAll)

	# 5.3. Extracting mobility data from B-post postal data

if (!file.exists("B-post_avisages_data/B-post_data_catchment_areas.csv"))
	{
		days2 = gsub("-","",days1[1:73]); days3 = ymd(days1) # TO ACTUALISE !!
		data = read.csv("B-post_avisages_data/B-post_data2_Daniele_081220.csv")
		colnames(data) = gsub("\\.","",gsub("X","",colnames(data)))
		bpostIndexCommunes = matrix(nrow=dim(communes2@data)[1], ncol=length(days1[1:73]))
		row.names(bpostIndexCommunes) = communes2@data[,"nouveau_PO"]; colnames(bpostIndexCommunes) = days1[1:73]
		for (i in 1:dim(communes2@data)[1])
			{
				for (j in 1:length(days2))
					{
						if ((length(which(data[,"PO"]==communes2@data[i,"nouveau_PO"])) > 0)&(length(which((colnames(data)==days2[j]))) > 0))
							{
								bpostIndexCommunes[i,j] = data[which(data[,"PO"]==communes2@data[i,"nouveau_PO"]),which((colnames(data)==days2[j]))]
							}
					}
			}
		bpostIndexCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1[1:73]))
		colnames(bpostIndexCatchments) = paste0("bpostIndexAll_",days1[1:73])
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						bpostIndexCatchments[i,j] = sum(bpostIndexCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(bpostIndexCatchments, "B-post_avisages_data/B-post_data_catchment_areas.csv", quote=F, row.names=F, sep=",")
	}

	# 5.4. Extracting temperature and precipitation data from CDS Copernicus

if (!file.exists("CDS_Copernicus_data/Catchment_areas_temp.csv"))
	{
		days2 = gsub("-","",days1); days3 = ymd(days1); communes3 = spTransform(communes2, CRS("+init=epsg:4326"))
		temperatureCommunes = matrix(nrow=dim(communes3@data)[1], ncol=length(days1))
		row.names(temperatureCommunes) = communes2@data[,"nouveau_PO"]; colnames(temperatureCommunes) = days1
		for (i in 1:length(days2))
			{
				month = substr(days2[i], 5, 6)
				rasts = brick(paste0("CDS_Copernicus_data/CDS_temperature_",month,"20.nc"))
				rast = rasts[[which(gsub("\\.","",gsub("X","",gsub(".12.00.00","",names(rasts))))==days2[i])]]
				for (j in 1:dim(communes3@data)[1])
					{
						maxArea = 0; polIndex = 0
						for (k in 1:length(communes3@polygons[[j]]@Polygons))
							{
								if (maxArea < communes3@polygons[[j]]@Polygons[[k]]@area)
									{
										maxArea = communes3@polygons[[j]]@Polygons[[k]]@area; polIndex = k
									}
							}
						pol = communes3@polygons[[j]]@Polygons[[polIndex]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sf::st_as_sfc(sps); st_crs(pol) = communes3@proj4string
						temperatureCommunes[j,i] = exact_extract(rast, pol, fun="mean")
					}
			}
		temperatureCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(temperatureCatchments) = paste0("temperature_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						temperatureCatchments[i,j] = sum(temperatureCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(temperatureCatchments, "CDS_Copernicus_data/Catchment_areas_temp.csv", quote=F, row.names=F, sep=",")
	}
if (!file.exists("CDS_Copernicus_data/Catchment_areas_rhum.csv"))
	{
		days2 = gsub("-","",days1); days3 = ymd(days1); communes3 = spTransform(communes2, CRS("+init=epsg:4326"))
		relativeHumidityCommunes = matrix(nrow=dim(communes3@data)[1], ncol=length(days1))
		row.names(relativeHumidityCommunes) = communes2@data[,"nouveau_PO"]; colnames(relativeHumidityCommunes) = days1
		for (i in 1:length(days2))
			{
				month = substr(days2[i], 5, 6)
				rasts = brick(paste0("CDS_Copernicus_data/CDS_relativeHumidity_",month,"20.nc"))
				rast = rasts[[which(gsub("\\.","",gsub("X","",gsub(".12.00.00","",names(rasts))))==days2[i])]]
				for (j in 1:dim(communes3@data)[1])
					{
						maxArea = 0; polIndex = 0
						for (k in 1:length(communes3@polygons[[j]]@Polygons))
							{
								if (maxArea < communes3@polygons[[j]]@Polygons[[k]]@area)
									{
										maxArea = communes3@polygons[[j]]@Polygons[[k]]@area; polIndex = k
									}
							}
						pol = communes3@polygons[[j]]@Polygons[[polIndex]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sf::st_as_sfc(sps); st_crs(pol) = communes3@proj4string
						relativeHumidityCommunes[j,i] = exact_extract(rast, pol, fun="mean")
					}
			}
		relativeHumidityCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(relativeHumidityCatchments) = paste0("relativeHumidity_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						relativeHumidityCatchments[i,j] = sum(relativeHumidityCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(relativeHumidityCatchments, "CDS_Copernicus_data/Catchment_areas_rhum.csv", quote=F, row.names=F, sep=",")
	}
if (!file.exists("CDS_Copernicus_data/Catchment_areas_prec.csv"))
	{
		days2 = gsub("-","",days1); days3 = ymd(days1); communes3 = spTransform(communes2, CRS("+init=epsg:4326"))
		precipitationCommunes = matrix(nrow=dim(communes3@data)[1], ncol=length(days1))
		row.names(precipitationCommunes) = communes2@data[,"nouveau_PO"]; colnames(precipitationCommunes) = days1
		for (i in 1:length(days2))
			{
				month = substr(days2[i], 5, 6)
				rasts = brick(paste0("CDS_Copernicus_data/CDS_precipitation_",month,"20.nc"))
				rast = rasts[[which(gsub("\\.","",gsub("X","",gsub(".12.00.00","",names(rasts))))==days2[i])]]
				for (j in 1:dim(communes3@data)[1])
					{
						maxArea = 0; polIndex = 0
						for (k in 1:length(communes3@polygons[[j]]@Polygons))
							{
								if (maxArea < communes3@polygons[[j]]@Polygons[[k]]@area)
									{
										maxArea = communes3@polygons[[j]]@Polygons[[k]]@area; polIndex = k
									}
							}
						pol = communes3@polygons[[j]]@Polygons[[polIndex]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sf::st_as_sfc(sps); st_crs(pol) = communes3@proj4string
						precipitationCommunes[j,i] = exact_extract(rast, pol, fun="mean")
					}
			}
		precipitationCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(precipitationCatchments) = paste0("precipitation_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						precipitationCatchments[i,j] = sum(precipitationCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(precipitationCatchments, "CDS_Copernicus_data/Catchment_areas_prec.csv", quote=F, row.names=F, sep=",")
	}
if (!file.exists("CDS_Copernicus_data/Catchment_areas_radia.csv"))
	{
		days2 = gsub("-","",days1); days3 = ymd(days1); communes3 = spTransform(communes2, CRS("+init=epsg:4326"))
		solarRadiationCommunes = matrix(nrow=dim(communes3@data)[1], ncol=length(days1))
		row.names(solarRadiationCommunes) = communes2@data[,"nouveau_PO"]; colnames(solarRadiationCommunes) = days1
		for (i in 1:length(days2))
			{
				month = substr(days2[i], 5, 6)
				rasts = brick(paste0("CDS_Copernicus_data/CDS_solarRadiation_",month,"20.nc"))
				rast = rasts[[which(gsub("\\.","",gsub("X","",gsub(".12.00.00","",names(rasts))))==days2[i])]]
				for (j in 1:dim(communes3@data)[1])
					{
						maxArea = 0; polIndex = 0
						for (k in 1:length(communes3@polygons[[j]]@Polygons))
							{
								if (maxArea < communes3@polygons[[j]]@Polygons[[k]]@area)
									{
										maxArea = communes3@polygons[[j]]@Polygons[[k]]@area; polIndex = k
									}
							}
						pol = communes3@polygons[[j]]@Polygons[[polIndex]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sf::st_as_sfc(sps); st_crs(pol) = communes3@proj4string
						solarRadiationCommunes[j,i] = exact_extract(rast, pol, fun="mean")
					}
			}
		solarRadiationCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(solarRadiationCatchments) = paste0("solarRadiation_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						solarRadiationCatchments[i,j] = sum(solarRadiationCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(solarRadiationCatchments, "CDS_Copernicus_data/Catchment_areas_radia.csv", quote=F, row.names=F, sep=",")
	}
temperatureCatchments = read.csv("CDS_Copernicus_data/Catchment_areas_temp.csv")
precipitationCatchments = read.csv("CDS_Copernicus_data/Catchment_areas_prec.csv")
relativeHumidityCatchments = read.csv("CDS_Copernicus_data/Catchment_areas_rhum.csv")
solarRadiationCatchments = read.csv("CDS_Copernicus_data/Catchment_areas_radia.csv")
catchmentAreas@data = cbind(catchmentAreas@data, temperatureCatchments, relativeHumidityCatchments, solarRadiationCatchments)

	# 5.5. Extracting PM 2.5 emission value from gridded data

if (!file.exists("Rasters_de_irceline_be/Catchment_areas_PM_25.csv"))
	{
		template = raster("PM25_mean_2017.asc")
		grid = shapefile("Rasters_de_irceline_be/Rast_2020_daily_grids.shp")
		pm25 = read.csv("Rasters_de_irceline_be/PM25_2020_daily_grids.csv", head=T)
		days2 = gsub("-","",days1); days3 = ymd(days1)
		pm25CatchmentAreas = matrix(nrow=dim(catchmentAreas)[1], ncol=length(days1))
		colnames(pm25CatchmentAreas) = days1
		for (i in 1:length(days1))
			{
				index = which(colnames(pm25)==paste0("X",gsub("-","",days1[i])))
				if (length(index1) == 1)
					{
						temp1 = grid; temp1@data$vals = pm25[,index]
						temp2 = SpatialPolygons(temp1@polygons)
						rast = rasterize(temp2, template, field=temp1@data$vals)
						for (j in 1:dim(catchmentAreas@data)[1])
							{
								maxArea = 0; polIndex = 0
								for (k in 1:length(catchmentAreas@polygons[[j]]@Polygons))
									{
										if (maxArea < catchmentAreas@polygons[[j]]@Polygons[[k]]@area)
											{
												maxArea = catchmentAreas@polygons[[j]]@Polygons[[k]]@area; polIndex = k
											}
									}
								pol = catchmentAreas@polygons[[j]]@Polygons[[polIndex]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sf::st_as_sfc(sps); st_crs(pol) = communes1@proj4string; crs(rast) = crs(pol)
								pm25CatchmentAreas[j,i] = exact_extract(rast, pol, fun="mean")
							}
					}
			}
		colnames(pm25CatchmentAreas) = paste0("emissionPM25_",days1)
		write.table(pm25CatchmentAreas, "Rasters_de_irceline_be/Catchment_areas_PM_25.csv", quote=F, row.names=F, sep=",")
	}
pm25CatchmentAreas = read.csv("Rasters_de_irceline_be/Catchment_areas_PM_25.csv", head=T)
catchmentAreas@data = cbind(catchmentAreas@data, pm25CatchmentAreas)

if (writingFiles) write.csv(catchmentAreas@data, "Catchment_areas_2.csv", row.names=F, quote=F)

	# 5.6. Exploring the Relation between temporal variables and new hospitalisations

library(ggridges); library(ggplot2); library(viridis); library(hrbrthemes); library(tibble)

catchmentAreas@data = read.csv("Catchment_areas_2.csv", head=T); h = 1; movingAverageOf7days = TRUE
temporalVariables = c("proximusIndexAll", "temperature", "relativeHumidity", "solarRadiation", "emissionPM25")
if (showingPlots)
	{
		temporalVariable = temporalVariables[h]; temporalVariable = "proximusIndexAll"
		temporalVariable = temporalVariables[h]; temporalVariable = "temperature"
		temporalVariable = temporalVariables[h]; temporalVariable = "relativeHumidity"
		temporalVariable = temporalVariables[h]; temporalVariable = "solarRadiation"
		temporalVariable = temporalVariables[h]; temporalVariable = "emissionPM25"
		colnames(catchmentAreas@data) = gsub("\\.","-",colnames(catchmentAreas@data))
		temporalValues = catchmentAreas@data[,which(grepl(temporalVariable,colnames(catchmentAreas@data)))]
		newHospitaCases = catchmentAreas@data[,which(grepl("newHospitaliCases",colnames(catchmentAreas@data)))]
		dailyRatioHosps = catchmentAreas@data[,which(grepl("dailyRatioHosCases",colnames(catchmentAreas@data)))]
		Ds = c(1:30); R2s = matrix(nrow=dim(dailyRatioHosps)[1], ncol=length(Ds)); cors = matrix(nrow=dim(dailyRatioHosps)[1], ncol=length(Ds))
		for (i in 1:length(Ds))
			{
				for (j in 1:dim(dailyRatioHosps)[1])
					{
						days2 = as.character(ymd(days1)-Ds[i]); days3 = as.character(ymd(days1)-1)
						indices = which(days2%in%gsub(paste0(temporalVariable,"_"),"",colnames(temporalValues)))
						x1 = t(temporalValues[j,paste0(paste0(temporalVariable,"_"),days2[indices])]); x2 = x1
						y1 = t(newHospitaCases[j,paste0("newHospitaliCases_",days1[indices])]); y2 = y1
						y1 = t(dailyRatioHosps[j,paste0("dailyRatioHosCases_",days1[indices])]); y2 = y1
						if (movingAverageOf7days == TRUE)
							{
								for (k in 1:length(x2))
									{
										indices = seq(k-3,k+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(x2)))]
										if (!is.na(x2[k])) x2[k] = mean(as.numeric(x1[indices]), na.rm=T)
										if (!is.na(as.numeric(y2[k]))) y2[k] = mean(as.numeric(y1[indices]), na.rm=T)
									}
							}	else		{
								x2 = x1; y2 = y1
							}
						y = as.numeric(y2); x = as.numeric(x2)
						y[is.infinite(y)] = NA; x = x[!is.na(y)]; y = y[!is.na(y)]
						y = y[!is.na(x)]; x = x[!is.na(x)]; # plot(x,y)
						if (length(y) > 0)
							{
								lr = lm("y ~ x"); R2s[j,i] = summary(lr)$r.square; cors[j,i] = cor(x,y,method="spearman")
							}
					}
			}
		correlations = c(); days = c()
		for (i in 1:dim(cors)[2])
			{
				correlations = c(correlations, cors[,i])
				days = c(days, rep(i,dim(cors)[1]))
			}
		days = days[!is.na(correlations)]; correlations = correlations[!is.na(correlations)]
		df = as_tibble(data.frame(cbind(as.factor(days),correlations)))
				ggplot(df, aes(x=correlations, y=as.factor(days), fill=..x..)) +
		 	   geom_density_ridges_gradient(scale=3, rel_min_height=0.01) +
		scale_fill_viridis(name="Temp. [F]", option="C") + # theme_ipsum() +
		theme(legend.position="none", panel.spacing=unit(0.1,"lines"), strip.text.x=element_text(size=7))
	}
if (showingPlots)
	{
		pdf("Catchment_areas_NEW.pdf", width=9, height=6) # dev.new(width=9, height=6)
		par(mfrow=c(6,1), mar=c(0,2.5,0,0), oma=c(2,0,0,0), lwd=0.2, col="gray30")
		temporalVariables = c("proximusIndexAll", "temperature", "relativeHumidity", "solarRadiation", "emissionPM25"); movingAverageOf7days = TRUE
		temporalVariable_names = c("mobility index", "temperature", "relative humidity", "solar radiation", "PM 2.5 emission")
		selectedDays = c(16:length(days1)); newHospitaliCases = catchmentAreas@data[,which(grepl("newHospitaliCases",colnames(catchmentAreas@data)))]
		selected_dates = decimal_date(ymd("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-15"))
		selected_dates_name = c("","01-04","01-05","01-06","01-07","01-08","01-09","01-10","01-11","")
		xMin = decimal_date(ymd("2020-03-15")); xMax = decimal_date(ymd("2020-11-30"))
		allValues = catchmentAreas@data[,paste0("newHospitaliCases_",days1[selectedDays])]
		yMin = min(allValues,na.rm=T); yMax = max(allValues,na.rm=T)
		if (movingAverageOf7days == TRUE)
			{
				allYvalues = c()
				for (i in 1:dim(catchmentAreas@data)[1])
					{
						y1 = catchmentAreas@data[i,paste0("newHospitaliCases_",days1[selectedDays])]; y2 = y1
						for (j in 1:length(y2))
							{
								indices = seq(j-3,j+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(y1)))]
								if (!is.na(y2[j])) y2[j] = mean(as.numeric(y1[indices]), na.rm=T)
							}
						allYvalues = c(allYvalues, y2)
					}
				yMin = min(as.numeric(allYvalues),na.rm=T); yMax = max(as.numeric(allYvalues[!is.infinite(as.numeric(allYvalues))]),na.rm=T)
			}	else		{
				yMin = min(allValues,na.rm=T); yMax = max(allValues[!is.infinite(allValues)],na.rm=T)
			}
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				y1 = catchmentAreas@data[i,paste0("newHospitaliCases_",days1[selectedDays])]; y2 = y1
				if (movingAverageOf7days == TRUE)
					{
						for (j in 1:length(y2))
							{
								indices = seq(j-3,j+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(y1)))]
								values = y1[indices]; values = values[which(!is.infinite(as.numeric(values)))]
								if (!is.na(y2[j])) y2[j] = mean(as.numeric(values), na.rm=T)
							}
					}	else	{
						y2 = y1
					}
				if (i == 1)
					{
						plot(decimal_date(ymd(days1[selectedDays])), y2, lwd=0.1, col="red", type="l", axes=F, ann=F, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
					}	else	{
						lines(decimal_date(ymd(days1[selectedDays])), y2, lwd=0.1, col="red")
					}
			}
		axis(side=2, lwd.tick=0.2, cex.axis=0.9, lwd=0.2, tck=-0.05, col.axis="gray30", mgp=c(0,0.35,0))
		title(ylab="hospitalisations", cex.lab=1.1, mgp=c(1.4,0,0), col.lab="gray30")
		for (i in 1:length(temporalVariables))
			{
				temporalValues = catchmentAreas@data[,which(grepl(temporalVariables[i],colnames(catchmentAreas@data)))]
				yMin = min(temporalValues,na.rm=T); yMax = max(temporalValues,na.rm=T)
				if (movingAverageOf7days == TRUE)
					{
						allYvalues = c()
						for (j in 1:dim(catchmentAreas@data)[1])
							{
								y1 = temporalValues[j,paste0(temporalVariables[i],"_",days1[selectedDays])]; y2 = y1
								for (k in 1:length(y2))
									{
										indices = seq(k-3,k+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(y1)))]
										if (!is.na(y2[k])) y2[k] = mean(as.numeric(y1[indices]), na.rm=T)
									}
								allYvalues = c(allYvalues, y2)
							}
						yMin = min(as.numeric(allYvalues),na.rm=T); yMax = max(as.numeric(allYvalues),na.rm=T)
					}
				for (j in 1:dim(catchmentAreas@data)[1])
					{
						y1 = temporalValues[j,paste0(temporalVariables[i],"_",days1[selectedDays])]; y2 = y1
						if (movingAverageOf7days == TRUE)
							{
								for (k in 1:length(y2))
									{
										indices = seq(k-3,k+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(y1)))]
										if (!is.na(y2[k])) y2[k] = mean(as.numeric(y1[indices]), na.rm=T)
									}
							}	else	{
								y2 = y1
							}
						if (j == 1)
							{
								plot(decimal_date(ymd(days1[selectedDays])), y2, lwd=0.1, col="gray50", type="l", axes=F, ann=F, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
							}	else	{
								lines(decimal_date(ymd(days1[selectedDays])), y2, lwd=0.1, col="gray50")
							}
					}
				axis(side=2, lwd.tick=0.2, cex.axis=0.9, lwd=0.2, tck=-0.05, col.axis="gray30", mgp=c(0,0.35,0))
				title(ylab=temporalVariable_names[i], cex.lab=1.1, mgp=c(1.4,0,0), col.lab="gray30")
				if (i == length(temporalVariables))
					{
						axis(side=1, lwd.tick=0.2, cex.axis=0.9, lwd=0.2, tck=-0.05, col.axis="gray30", mgp=c(0,0.25,0), at=selected_dates, labels=selected_dates_name)
					}
			}
		dev.off()
	}

	# 5.7. Multivariate regression analyses considering different fixed lag times

D = 20; temporalVariables = c("proximusIndexAll", "temperature", "relativeHumidity", "solarRadiation", "emissionPM25")
betas = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(temporalVariables)); colnames(betas) = temporalVariables
pValues = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(temporalVariables)); colnames(pValues) = temporalVariables
for (i in 1:dim(catchmentAreas)[1])
	{
		days2 = as.character(ymd(days1)-D); days3 = as.character(ymd(days1)-1)
		dailyRatioHosps = catchmentAreas@data[,which(grepl("dailyRatioHosCases",colnames(catchmentAreas@data)))]
		responseV = t(dailyRatioHosps[i,paste0("dailyRatioHosCases_",days1)])
		predictors = matrix(nrow=length(days1), ncol=length(temporalVariables))
		for (j in 1:dim(predictors)[2])
			{
				for (k in 1:dim(predictors)[1])
					{
						if (paste0(paste0(temporalVariables[j],"_"),days2[k])%in%colnames(catchmentAreas@data))
							{
								predictors[k,j] = catchmentAreas@data[i,paste0(paste0(temporalVariables[j],"_"),days2[k])]
							}
					}		
			}
		indices = c()
		for (j in 1:dim(predictors)[1])
			{
				if ((sum(!is.na(predictors[j,])) == dim(predictors)[2])&(!is.na(responseV[j,1]))&(!is.infinite(responseV[j,1]))) indices = c(indices, j)
			}
		predictors = predictors[indices,]; responseV = responseV[indices,1]
		df = cbind(responseV, predictors); colnames(df) = c("dailyRatioHosCases",temporalVariables)
		for (j in 1:dim(df)[2]) df[,j] = (df[,j]-mean(df[,j]))/sd(df[,j])
		fm = paste0("dailyRatioHosCases ~ ",temporalVariables[1])
		for (j in 2:length(temporalVariables)) fm = paste0(fm," + ",temporalVariables[j])
		glm = glm(as.formula(fm), data=as.data.frame(df))
		for (j in 1:length(temporalVariables)) betas[i,j] = glm$coefficients[temporalVariables[j]]
		for (j in 1:length(temporalVariables)) pValues[i,j] = summary(glm)$coefficients[temporalVariables[j],"Pr(>|t|)"]
	}

