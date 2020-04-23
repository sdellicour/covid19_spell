library(ade4)
library(adephylo)
library(ape)
library(castor)
library(diagram)
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
library(rgeos)
library(seraphim)
library(spdep)
library(treeio)

# 1. Analyses of Sciensano data (generating tables)
# 2. Analyses of mobility data (aggregated mobile data)
# 3. Analyses at the province levels (hospitalisations)
# 4. Analyses at the commune levels (positive cases)
# 5. Analyses of hospital catchment areas
# 6. Phylogenetic and phylogeographic analyses

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
dates = unique(data_cases$DATE)[!is.na(unique(data_cases$DATE))]
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
		tab[i,"Tests"] = data_tests[index,"TESTS"]
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
		dates = unique(data_cases$DATE)[!is.na(unique(data_cases$DATE))]
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

# 2. Analyses of mobility data (aggregated mobile data)

communes = shapefile("Shapefile_communes/Shapefile_communes.shp")
fileNames = c("Mobility_data_mean_20200217_20200221","Mobility_data_mean_20200222_20200223",			  "Mobility_data_mean_20200316_20200320","Mobility_data_mean_20200321_20200322",			  "Mobility_data_mean_20200323_20200327","Mobility_data_mean_20200328_20200329")
periods = c("17/02 - 21/02/2020","16/03 - 20/03/2020","23/03 - 27/03/2020","22/02 - 23/02/2020","21/03 - 22/03/2020","28/03 - 29/03/2020",
			"17/02 - 21/02/2020","16/03 - 20/03/2020","23/03 - 27/03/2020","22/02 - 23/02/2020","21/03 - 22/03/2020","28/03 - 29/03/2020")
communes@data$period1In = rep(NA, dim(communes@data)[1]); communes@data$period1Out = rep(NA, dim(communes@data)[1])
communes@data$period2In = rep(NA, dim(communes@data)[1]); communes@data$period2Out = rep(NA, dim(communes@data)[1])
communes@data$period3In = rep(NA, dim(communes@data)[1]); communes@data$period3Out = rep(NA, dim(communes@data)[1])
communes@data$period4In = rep(NA, dim(communes@data)[1]); communes@data$period4Out = rep(NA, dim(communes@data)[1])
communes@data$period5In = rep(NA, dim(communes@data)[1]); communes@data$period5Out = rep(NA, dim(communes@data)[1])
communes@data$period6In = rep(NA, dim(communes@data)[1]); communes@data$period6Out = rep(NA, dim(communes@data)[1])
for (i in 1:length(fileNames))
	{
		data = read.csv(paste0("Data_mobility_Dahlberg/",fileNames[i],".csv"))
		for (k in 1:dim(communes@data)[1])
			{
				indices = which(data[,"trip_nis"]==communes@data[k,"NIS5"])
				if (length(indices) > 0)
					{
						communes@data[k,paste0("period",i,"In")] = sum(data[indices,"trips"]/data[indices,"subs"])
					}
				indices = which(data[,"home_nis"]==communes@data[k,"NIS5"])
				if (length(indices) > 0)
					{
						communes@data[k,paste0("period",i,"Out")] = sum(data[indices,"trips"]/data[indices,"subs"])
					}
			}
	}
if (showingPlots)
	{
		for (i in 1:length(variables))
			{
				communes@data[is.na(communes@data[,variables[i]]),variables[i]] = 0
			}
		variables = c("period1In","period3In","period5In","period2In","period4In","period6In",
					  "period1Out","period3Out","period5Out","period2Out","period4Out","period6Out")
		variableNames = c("# trips in / # subscribers","# trips in / # subscribers","# trips in / # subscribers",
						  "# trips in / # subscribers","# trips in / # subscribers","# trips in / # subscribers",
						  "# trips out / # subscribers","# trips out / # subscribers","# trips out / # subscribers",
						  "# trips out / # subscribers","# trips out / # subscribers","# trips out / # subscribers")
		communes_light = gSimplify(communes, 100); colourScales = list()
		for (i in 1:6) colourScales[[i]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101])
		for (i in 7:12) colourScales[[i]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlOrBr"))(131)[1:101])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				if (i == 1)
					{
						minV = min(communes@data[,variables[1:6]], na.rm=T); maxV = max(communes@data[,variables[1:6]], na.rm=T)
					}
				if (i == 7)
					{
						minV = min(communes@data[,variables[7:12]], na.rm=T); maxV = max(communes@data[,variables[7:12]], na.rm=T)
					}	
				values = communes@data[,variables[i]]
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]; legendRast = raster(as.matrix(c(minV,maxV)))		
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(communes_light, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames[i], cex=0.5, col="gray30", at=92000, line=-11.9)
				mtext(periods[i], cex=0.5, col="gray30", at=92000, line=-12.6)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}
	}
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
communes@data[which(communes@data[,"TX_PROV_DE"]=="NA"),"TX_PROV_DE"] = "Bruxelles"
communes@data[,"TX_PROV_DE"] = gsub("Anvers","Antwerpen",gsub("Limbourg","Limburg",communes@data[,"TX_PROV_DE"]))
communes@data[,"TX_PROV_DE"] = gsub("Province du ","",gsub("Province de ","",gsub("Province d\u0092","",communes@data[,"TX_PROV_DE"])))	
communes@data[,"TX_PROV_DE"] = gsub("Brabant flamand","Vlaams Brabant",gsub("Brabant wallon","Brabant Wallon",communes@data[,"TX_PROV_DE"]))
communes@data[,"TX_PROV_DE"] = gsub("Flandre occidentale","West-Vlaanderen",gsub("Flandre orientale","Oost-Vlaanderen",communes@data[,"TX_PROV_DE"]))	
provinces@data$period1In = rep(NA, dim(provinces@data)[1]); provinces@data$period1Out = rep(NA, dim(provinces@data)[1])
provinces@data$period2In = rep(NA, dim(provinces@data)[1]); provinces@data$period2Out = rep(NA, dim(provinces@data)[1])
provinces@data$period3In = rep(NA, dim(provinces@data)[1]); provinces@data$period3Out = rep(NA, dim(provinces@data)[1])
provinces@data$period4In = rep(NA, dim(provinces@data)[1]); provinces@data$period4Out = rep(NA, dim(provinces@data)[1])
provinces@data$period5In = rep(NA, dim(provinces@data)[1]); provinces@data$period5Out = rep(NA, dim(provinces@data)[1])
provinces@data$period6In = rep(NA, dim(provinces@data)[1]); provinces@data$period6Out = rep(NA, dim(provinces@data)[1])
for (i in 1:length(fileNames))
	{
		data = read.csv(paste0("Data_mobility_Dahlberg/",fileNames[i],".csv"))
		for (j in 1:dim(provinces@data)[1])
			{
				INs = 0; OUTs = 0
				indices1 = which(communes@data[,"TX_PROV_DE"]==provinces@data[j,"NAME_2"])
				for (k in 1:length(indices1))
					{
						indices2 = which(data[,"trip_nis"]==communes@data[indices1[k],"NIS5"])
						if (length(indices2) > 0)
							{
								for (l in 1:length(indices2))
									{
										if ((!is.na(data[indices2[l],"trips"]))&(!is.na(data[indices2[l],"subs"])))
											{
												INs = INs + (data[indices2[l],"trips"]/data[indices2[l],"subs"])
											}
									}
							}
						indices3 = which(data[,"home_nis"]==communes@data[indices1[k],"NIS5"])
						if (length(indices3) > 0)
							{
								for (l in 1:length(indices3))
									{
										if ((!is.na(data[indices3[l],"trips"]))&(!is.na(data[indices3[l],"subs"])))
											{
												OUTs = OUTs + (data[indices3[l],"trips"]/data[indices3[l],"subs"])
											}
									}
							}
					}
				provinces@data[j,paste0("period",i,"In")] = INs
				provinces@data[j,paste0("period",i,"Out")] = OUTs
			}
	}
df = provinces@data[,c("period1In","period2In","period3In","period4In","period5In","period6In",
				"period1Out","period2Out","period3Out","period4Out","period5Out","period6Out")]
row.names(df) = provinces@data[,"NAME_2"]; row.names(df) = gsub("Bruxelles","Brussels",gsub(" ","",row.names(df)))
if (writingFiles) write.csv(df, "Provinces_mobile_data_IN_OUT_indices.csv", quote=F)
if (showingPlots)
	{
		for (i in 1:length(variables))
			{
				provinces@data[is.na(communes@data[,variables[i]]),variables[i]] = 0
			}
		variables = c("period1In","period3In","period5In","period2In","period4In","period6In",
					  "period1Out","period3Out","period5Out","period2Out","period4Out","period6Out")
		variableNames = c("# trips in / # subscribers","# trips in / # subscribers","# trips in / # subscribers",
						  "# trips in / # subscribers","# trips in / # subscribers","# trips in / # subscribers",
						  "# trips out / # subscribers","# trips out / # subscribers","# trips out / # subscribers",
						  "# trips out / # subscribers","# trips out / # subscribers","# trips out / # subscribers")
		colourScales = list()
		for (i in 1:6) colourScales[[i]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101])
		for (i in 7:12) colourScales[[i]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlOrBr"))(131)[1:101])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				if (i == 1)
					{
						minV = min(provinces@data[,variables[1:6]], na.rm=T); maxV = max(provinces@data[,variables[1:6]], na.rm=T)
					}
				if (i == 7)
					{
						minV = min(provinces@data[,variables[7:12]], na.rm=T); maxV = max(provinces@data[,variables[7:12]], na.rm=T)
					}	
				values = provinces@data[,variables[i]]
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]; legendRast = raster(as.matrix(c(minV,maxV)))		
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(provinces, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames[i], cex=0.5, col="gray30", at=92000, line=-11.9)
				mtext(periods[i], cex=0.5, col="gray30", at=92000, line=-12.6)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}
	}

fileNames = c("Mobility_data_mean_20200217_20200221","Mobility_data_mean_20200316_20200320","Mobility_data_mean_20200323_20200327")
periods = c("17/02-21/02/2020","16/03-20/03/2020","23/03-27/03/2020")
for (i in 1:length(fileNames))
	{
		data = read.csv(paste0("Data_mobility_Dahlberg/",fileNames[i],".csv"))
		NIScodes = unique(data[,"home_nis"]); NIScodes = NIScodes[order(NIScodes)]
		nTrips = matrix(0, nrow=length(NIScodes), ncol=length(NIScodes))
		for (j in 2:dim(nTrips)[1])
			{
				for (k in 1:(j-1))
					{
						index = which((data[,"home_nis"]==NIScodes[j])&(data[,"trip_nis"]==NIScodes[k]))
						if (length(index) == 1) nTrips[j,k] = data[index,"trips"]
					}
			}
		locations = matrix(nrow=length(NIScodes), ncol=3)
		colnames(locations) = c("ind","x","y"); locations[,1] = NIScodes
		for (j in 1:dim(locations)[1])
			{
				index1 = which(communes@data[,"NIS5"]==locations[j,1])
				if (length(index1) == 1)
					{
						maxArea = 0; index2 = 0
						for (k in 1:length(communes@polygons[[index1]]@Polygons))
							{
								if (maxArea < communes@polygons[[index1]]@Polygons[[k]]@area)
									{
										maxArea = communes@polygons[[index1]]@Polygons[[k]]@area; index2 = k
									}
							}
						pol = communes@polygons[[index1]]@Polygons[[index2]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sps; proj4string(pol) = communes@proj4string
						locations[j,c("x","y")] = coordinates(pol)
					}
			}
		indices = which(!is.na(locations[,"x"]))
		locations = data.frame(locations[indices,])
		nTrips = as.matrix(nTrips[indices,indices])
		row.names(nTrips) = locations$ind; colnames(nTrips) = locations$ind
		hw = MAPI_EstimateHalfwidth(locations, crs=2154, beta=0.25)
		grid = MAPI_GridHexagonal(locations, crs=31370, hw=1000, buf=0, shift=F)
		mapi = MAPI_RunOnGrid(locations, metric=nTrips, grid, isMatrix=T, nbPermuts=1000)
		tails = MAPI_Tails(mapi, alpha=0.05)
		if (showingPlots) MAPI_Plot2(mapi, tails=tails)
		st_write(mapi, dsn="./Data_mobility_Dahlberg", layer=paste0(fileNames[i],"_1"), driver="ESRI Shapefile", update=T, delete_layer=T)
		st_write(tails, dsn="./Data_mobility_Dahlberg", layer=paste0(fileNames[i],"_2"), driver="ESRI Shapefile", update=T, delete_layer=T)		
	}
belgium = spTransform(raster::getData("GADM", country="BEL", level=0), crs(communes))
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
population = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(communes))
template = population; template[!is.na(template[])] = 0
for (i in 1:length(fileNames))
	{
		mapi = shapefile(paste0("Data_mobility_Dahlberg/",fileNames[i],"_1.shp"))
		rast = rasterize(mapi, template, field="avg_value")
		rast = raster::mask(crop(rast,belgium),belgium)
		writeRaster(rast, paste0("Data_mobility_Dahlberg/",fileNames[i],".tif"))
	}
if (showingPlots)
	{
		rasts = list()
		for (i in 1:length(fileNames))
			{
				rasts[[i]] = raster(paste0("Data_mobility_Dahlberg/",fileNames[i],".tif"))
				rasts[[i]][!is.na(rasts[[i]][])] = log(rasts[[i]][!is.na(rasts[[i]][])]+1)
				if (i == 1)
					{
						minV = min(rasts[[i]][], na.rm=T); maxV = max(rasts[[i]][], na.rm=T)
						# r = rasts[[i]]; r[!is.na(r)] = 0; contour = rasterToPolygons(r, dissolve=T)
					}
			}
		colourScale = colorRampPalette(brewer.pal(9,"BrBG"))(161)[31:131]; cols = list()
		dev.new(width=9,height=2.7); legendRast = raster(as.matrix(c(minV,maxV)))
		par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(rasts))
			{
				index1 = (((min(rasts[[i]][],na.rm=T)-minV)/(maxV-minV))*100)+1
				index2 = (((max(rasts[[i]][],na.rm=T)-minV)/(maxV-minV))*100)+1
				cols[[i]] = colourScale[index1:index2]
				plot(rasts[[i]], col=cols[[i]], ann=F, axes=F, box=F, legend=F)
				mtext(paste0("Number of trips (log)"), cex=0.54, col="gray30", at=55000, line=-13.3)
				mtext(paste0(periods[i]), cex=0.54, col="gray30", at=55000, line=-14.2)
				plot(legendRast, legend.only=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0,0.47,0.10,0.12),
	 				 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
	 				 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.2,0)))
			}
	}

# 3. Analyses at the province levels (hospitalisations)

periods = c("16/03-22/03/2020","23/03-29/03/2020","30/03-05/04/2020")
selectedDays1 = ymd(c("2020-03-22","2020-03-29","2020-04-05"))
firstDay = ymd("2020-01-30"); D = 7 # time interval
selectedDays2 = as.numeric(selectedDays1-firstDay)
provinces = raster::getData("GADM", country="BEL", level=2)
provinces@data$NAME_3 = c("Brussels","Antwerpen","Limburg","OostVlaanderen","VlaamsBrabant",
				 "WestVlaanderen","BrabantWallon","Hainaut","Li\xe8ge","Luxembourg","Namur")
data = read.csv("Data_Sciensano_0604/COVID19BE_HOSP.csv")
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
				if (length(index1) == 0)	
					{
						QTD1 = 0
					}	else	{
						QTD1 = temp[index1,"TOTAL_IN"]
					}
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

if (showingPlots)
	{
		DTmax = ceiling(max(provinces@data[,c("DT1","DT2","DT3")])); DTmax = 30
		colourScale = colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101]; cols = list()
		dev.new(width=3.2,height=7); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(selectedDays2))
			{
				values = provinces@data[,paste0("DT",i)]
				values[which(values>DTmax)] = DTmax
				cols[[i]] = colourScale[1+((values/DTmax)*100)]
				plot(provinces, border="gray30", col=cols[[i]])
				mtext(paste0("Doubling time - ",periods[i]), cex=0.54, col="gray30", at=3.55, line=-14.2)
				plot(legendRast, legend.only=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
	 				 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
	 				 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.2,0), at=seq(0,DTmax,5), labels=c("0","5","10","15","20","25","30")))
			}
		DT_values = c(0,1,2,4,8,16,32); DTmax = max(DT_values); cols = list()
		colourScale1 = colorRampPalette(brewer.pal(9,"YlGn"))(12)[1:length(DT_values)]
		colourScale1 = colorRampPalette(brewer.pal(11,"RdYlGn"))(12)[3:(length(DT_values)+2)]
		colourScale2 = c("gray90",colourScale1)
		dev.new(width=9,height=2.7); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(2,2,2,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(selectedDays2))
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
tab = provinces@data[,c("DT1","DT2","DT3")]; row.names(tab) = provinces@data$NAME_2; colnames(tab) = periods
if (writingFiles) write.csv(round(tab,2), "Doubling_times_provinces.csv", quote=F)

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

if (showingPlots)
	{
		cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#d3d3d3")
		DTmax = ceiling(max(c(max(doublingTHosCases))))
		dev.new(width=7,height=5); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(1,1), mar=c(2.9,3.1,1,1), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		xLabels = c(paste0(c(22:31),"-03"),paste0("0",c(1:5),"-04")); dates = c(1:length(xLabels))
		for (i in 1:dim(doublingTHosCases)[1])
			{
				if (i == 1)
					{
						plot(dates,doublingTHosCases[i,], col=cols[i], lwd=1, ylim=c(1.8,33), axes=F, ann=F, type="l")
					}	else	{
						lines(dates,doublingTHosCases[i,], col=cols[i], lwd=1)
					}
				points(dates,doublingTHosCases[i,], col=cols[i], cex=0.7, pch=16)
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.013, col.axis="gray30", mgp=c(0,0.05,0), at=c(-1,dates), labels=c(-1,xLabels))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-2,32,2))
		title(ylab="doubling time hospitalisation (time window = 7 days)", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		title(xlab="day", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
		legend(1, 33, provinces@data$NAME_2, col=cols, text.col="gray30", pch=16, pt.cex=1.2, box.lty=0, cex=0.7, y.intersp=1.3)
	}

# 4. Analyses at the commune levels (positive cases)

communes = shapefile("Shapefile_communes/Shapefile_communes.shp")
communes_light = gSimplify(communes, 100)
equivalence = read.csv("Shapefile_communes/Postal_codes_vs_NIS.csv", header=T)

	# 4.1. Computing positive cases doubling times for two time periods

periods = c("18-26/03/2020","27/03-04/04/2020")
selectedDays1 = ymd(c("2020-03-26","2020-04-04"))
firstDay = ymd("2020-01-30"); D = 9 # time interval
selectedDays2 = as.numeric(selectedDays1-firstDay)
data = read.csv("Google_Drive_N_Hens/Data_Hospit_1_05-04.csv", sep=";")
daysSinceTheFirstCase = seq(1,200,1)
cumulatedCases_list = list()
for (i in 1:dim(communes@data)[1])
	{
		cumulatedCases = rep(0, length(daysSinceTheFirstCase))
		NIS = communes@data[i,"NIS5"]
		indices1 = which(equivalence[,"code_INS"]==NIS)
		if (length(indices1) != 0)
			{
				indices2 = c()
				for (j in 1:length(indices1))
					{
						postalCode = equivalence[indices1[j],"code_Postal"]
						indices2 = c(indices2, which(data[,"Postcode"]==postalCode))
					}
				if (length(indices2) != 0)
					{
						temp = data[indices2,]
						temp = temp[which(temp[,"dateusedforstatistics"]!=""),]
						dates = temp[,"dateusedforstatistics"]
						dates = gsub("Jan","-01-",dates)
						dates = gsub("Feb","-02-",dates)
						dates = gsub("Mar","-03-",dates)
						dates = gsub("Apr","-04-",dates)
						days = as.numeric(dmy(dates)-firstDay)
						for (j in 1:length(daysSinceTheFirstCase))
							{
								cumulatedCases[j] = sum(days <= j) 
							}
					}
			}
		cumulatedCases_list[[i]] = cumulatedCases
	}
communes@data$cases18March = rep(NA, dim(communes@data)[1])
communes@data$cases26March = rep(NA, dim(communes@data)[1])
communes@data$cases27March = rep(NA, dim(communes@data)[1])
communes@data$cases04April = rep(NA, dim(communes@data)[1])
communes@data$DT1 = rep(NA,dim(communes@data)[1])
communes@data$DT2 = rep(NA,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		temp = cumulatedCases_list[[i]]
		for (j in 1:length(selectedDays2))
			{
				QTD1 = temp[selectedDays2[j]-D]
				QTD2 = temp[selectedDays2[j]]
				if (j == 1)
					{
						communes@data[i,"cases18March"] = QTD1
						communes@data[i,"cases26March"] = QTD2
					}
				if (j == 2)
					{
						communes@data[i,"cases27March"] = QTD1
						communes@data[i,"cases04April"] = QTD2
					}
				if ((QTD1 >= 5)&(QTD2 >= 5)&(QTD1 != QTD2))
					{
						DT = (D*log(2))/(log(QTD2/QTD1))
					}	else	{
						DT = NA
					}
				communes@data[i,paste0("DT",j)] = DT
			}
	}
if (showingPlots)
	{
		DTmax = ceiling(max(communes@data[,c("DT1","DT2")], na.rm=T)); DTmax = 32
		colourScale1 = colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101]; cols = list()
		colourScale2 = c("gray90",colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101])
		dev.new(width=7,height=3); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(1,1,1,1), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(selectedDays2))
			{
				values = communes@data[,paste0("DT",i)]
				values[is.na(values)] = 0; values[values[]>30] = 30
				cols[[i]] = colourScale2[1+((values/DTmax)*100)]
				plot(communes_light, border="gray30", col=cols[[i]], lwd=0.1)
				mtext(paste0("Doubling time - ",periods[i]), cex=0.6, col="gray30", at=90000, line=-11.3)
				plot(legendRast, legend.only=T, col=colourScale1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.55, lwd=0,
			 		 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.0,0)))
			}
		DT_values = c(0,1,2,4,8,16,32); DTmax = max(DT_values)
		colourScale1 = colorRampPalette(brewer.pal(9,"YlGn"))(12)[1:length(DT_values)]; colourScale2 = c("gray90",colourScale1)
		dev.new(width=7,height=3); par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(1,1,1,1), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(selectedDays2))
			{
				values1 = communes@data[,paste0("DT",i)]; values2 = communes@data[,paste0("DT",i)]
				values2[values1[]>32] = 7; values2[values1[]<=32] = 6; values2[values1[]<=16] = 5
				values2[values1[]<=8] = 4; values2[values1[]<=4] = 3; values2[values1[]<=2] = 2
				values2[values1[]<=1] = 1; values2[is.na(values1)] = 0; cols[[i]] = colourScale2[values2+1]
				plot(communes_light, border="gray30", col=cols[[i]], lwd=0.1)
				mtext(paste0("Doubling time"), cex=0.65, col="gray30", at=70000, line=-10.3)
				mtext(paste0(periods[i]), cex=0.65, col="gray30", at=70000, line=-11)
			}
		legend(290000, 250000, c("0 - 1 day","1 - 2 days","2 - 4 days","4 - 8 days","8 - 16 days","16 - 32 days",">32 days"),
			   col=colourScale1, text.col="gray30", pch=16, pt.cex=1.0, box.lty=0, cex=0.65, y.intersp=1.1)
	}

	# 4.2. Extracting and assigning covariate values to each commune

communes@data$xCentroid = rep(0,dim(communes@data)[1])
communes@data$yCentroid = rep(0,dim(communes@data)[1])
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
		centroidCoordinates = coordinates(pol)
		communes@data[i,"xCentroid"] = centroidCoordinates[1,1]
		communes@data[i,"yCentroid"] = centroidCoordinates[1,2]
	}
data = read.csv("Data_Sciensano_1704/COVID19BE_CASES_MUNI_CUM.csv")
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
pm25 = raster("PM25_anmean_2017.asc")
parks = raster("Urban_parcs_UNamur.tif")
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
		communes@data[i,"pm10"] = exact_extract(pm10, sf::st_as_sfc(pol_light), fun='mean')
		communes@data[i,"pm25"] = exact_extract(pm10, sf::st_as_sfc(pol_light), fun='mean')
		urban_park_areas = exact_extract(pm10, sf::st_as_sfc(pol_light), fun='sum')
		communes@data[i,"ratioParksPopulation"] = urban_park_areas/communes@data[i,"population"]
	}
if (!file.exists("CorineLandCover.asc"))
	{
		clc = crop(raster("CorineLandCover18.tif"), extent(3500000,4500000,2500000,4000000))
		clc_crs = clc@crs; communes_clc = spTransform(communes, clc_crs)
		clc = mask(crop(clc, communes_clc), communes_clc)
		clc@crs = clc_crs; writeRaster(clc, "CorineLandCover.asc")
	}
if (!file.exists("Communes_urbanP.csv"))
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
		write.csv(propUrbanArea, "Shapefile_communes/Proportions_urbanA.csv", row.names=F, quote=F)
	}
communes@data$propUrbanArea = read.csv("Shapefile_communes/Proportions_urbanA.csv")[,3]

variables = c("popDensityLog","medianIncome","sectorP","sectorS","sectorT",
			  "medianAge","moreThan65","pm10","pm25","propUrbanArea"); dfs = list()
df1 = communes@data[,c("NIS5","xCentroid","yCentroid","DT1",variables)]; df1 = df1[!is.na(df1[,"DT1"]),]; dfs[[1]] = df1
df2 = communes@data[,c("NIS5","xCentroid","yCentroid","DT2",variables)]; df2 = df2[!is.na(df2[,"DT2"]),]; dfs[[2]] = df2
if (writingFiles == TRUE)
	{
		df = communes@data[,c("NIS5","xCentroid","yCentroid","cases18March","cases26March",
							  "cases27March","cases04April","DT1","DT2",variables)]
		write.csv(df, "Covariate_values_commune.csv", row.names=F, quote=F)		
	}

	# 4.3. Plotting the doubling time estimates and each covariate

if (showingPlots)
	{
		variables = c("incidences","DT1","DT2","medianIncome","propUrbanArea","popDensityLog",
					  "sectorP","sectorS","sectorT","moreThan65","pm10","pm25")
		variableNames = c("# cases per 1000 persons","Doubling times 1° period","Doubling times 2° period",
						  "Median declared income (€)","Urban area proportion","Population density (log)",
						  "% in primary sector","% in secundary sector","% in tertiary sector",
						  ">= 65 years (proportion)","PM 1.0 emission","PM 2.5 emission")
		communes_light = gSimplify(communes, 100); colourScales = list()
		colourScales[[1]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101])
		colourScales[[2]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101])
		colourScales[[3]] = c("#E5E5E5",colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101])
		colourScales[[4]] = c(colorRampPalette(brewer.pal(9,"RdPu"))(151)[1:101])
		colourScales[[5]] = c(colorRampPalette(brewer.pal(9,"Purples"))(151)[1:101])
		colourScales[[6]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(151)[1:101])
		colourScales[[7]] = c(colorRampPalette(brewer.pal(9,"Greens"))(151)[1:101])
		colourScales[[8]] = c(colorRampPalette(brewer.pal(9,"Oranges"))(151)[1:101])
		colourScales[[9]] = c(colorRampPalette(brewer.pal(9,"Blues"))(151)[1:101])
		colourScales[[10]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
		colourScales[[11]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(151)[1:101])
		colourScales[[12]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(151)[1:101])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (i in 1:length(variables))
			{
				values = communes@data[,variables[i]]
				if ((i == 2)|(i == 3))
					{
						values[is.na(values)] = 0; values[values[]>30] = 30
					}
				minV = min(values); maxV = max(values)
				legendCols = colourScales[[i]][1:length(colourScales[[i]])]; legendRast = raster(as.matrix(c(minV,maxV)))		
				cols = colourScales[[i]][(((values-minV)/(maxV-minV))*100)+1]
				plot(communes_light, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames[i], cex=0.54, col="gray30", at=92000, line=-12.4)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}
	}

	# 4.4. Performing and plotting the first axes of an exploratory PCA

if (showingPlots)
	{
		pca = dudi.pca(df1[,variables], scannf=F, nf=length(variables)); lis = pca$li[,1:2]; cos = pca$co
		colourScale = colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101]
		DTmax = max(df1[,"DT1"]); legendRast = raster(as.matrix(seq(0,DTmax,1)))
		cols = colourScale[1+((df1[,"DT1"]/DTmax)*100)]
		dev.new(width=6, height=6); par(mar=c(3,3,1.5,1.5), lwd=0.2, col="gray30")
		plot(lis, col="gray50", cex=0.3, pch=16, ann=F, axes=F, xlim=c(-5.5,7.5), ylim=c(-5.5,4.0))
		points(lis, col="gray30", cex=0.75, pch=1, lwd=0.3); points(lis, col=cols, cex=0.70, pch=16); 
		s.corcircle(2*cos, xax=1, yax=2, box=F, sub="", csub=0.7, clabel=0.7, possub="topleft", grid=F, cgrid=1, full=F, add.plot=T)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-9,9,1))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-7,9,1))
		title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
		mtext(paste0("Doubling time - 18-26/03/2020"), cex=0.75, col="gray30", at=2.5, line=-1)
		plot(legendRast, legend.only=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.73,0.90,0.91),
			 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.55, lwd=0,
			 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.0,0)))
	}

	# 4.5. Assessing spatial autocorrelation with the Moran's I test

for (i in 1:length(dfs))
	{
		geoDists = as.matrix(dist(dfs[[i]][,c("xCentroid","yCentroid")]))
		weights = 1/geoDists; diag(weights) = 0
		print(Moran.I(dfs[[i]][,paste0("DT",i)], weights))
	}

	# 4.6. Univariate (LR) followed by multivariate regression (GLM) analyses

selectedVariables = list()
for (i in 1:length(dfs))
	{
		buffer = c()
		for (j in 1:length(variables))
			{
				formula = paste0("DT",i," ~ ",variables[j])
				lr = glm(formula, data=dfs[[i]])
				pValue = summary(lr)$coefficients[2,4]
				if (pValue < 0.05)
					{
						buffer = c(buffer, variables[j])
					}
			}
		selectedVariables[[i]] = buffer
	}
for (i in 1:length(dfs))
	{
		df_z = dfs[[i]]
		for (j in 1:dim(dfs[[i]])[2])
			{
				df_z[,j] = zTransformation(df_z[,j])
			}
		formula = paste0("DT",i," ~ ",selectedVariables[[i]][1])
		if (length(selectedVariables[[i]]) > 1)
			{
				for (j in 2:length(selectedVariables[[i]]))
					{
						formula = paste0(formula," + ",selectedVariables[[i]][j])
					}
			}
		glm = glm(formula, data=df_z); print(summary(glm))
	}

	# 4.7. GAM (generalised additive model) analyses

gams = list()
zTransformations = FALSE
for (i in 1:length(dfs))
	{
		df = dfs[[i]]; colnames(df) = gsub(paste0("DT",i),"DT",colnames(df))
		if (zTransformations == TRUE)
			{
				for (j in 1:dim(dfs[[i]])[2]) df[,j] = zTransformation(df[,j])
			}
		gam = gam(DT ~ s(popDensityLog) + s(propUrbanArea)+ s(medianIncome) + s(sectorP) + s(sectorS) + s(sectorT)
					 + s(medianAge) + s(moreThan65) + s(pm10) + s(pm25) + s(xCentroid,yCentroid), data=df, method="REML")
		print(summary(gam)); gams[[i]] = gam
		if (showingPlots)
			{
				dev.new(); plot(gam, pages=1)
			}
	}
if (showingPlots)
	{
		gam = gams[[1]]; responseCurves = list()
		curves = plot(gam, pages=1); dev.off()
		selectedVariables = c("popDensityLog","propUrbanArea")
		variableNames = c("population density (log)","urban area proportion")
		dev.new(width=6.5,height=3)
		par(mfrow=c(1,2), mar=c(3,3,1,1), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, col="gray30", bty="o")
		for (i in 1:length(selectedVariables))
			{
				index = which(colnames(dfs[[1]])==selectedVariables[i])-4
				lower_l = curves[[index]]$fit-curves[[index]]$se
				upper_l = curves[[index]]$fit+curves[[index]]$se
				yLim = c(min(c(lower_l,upper_l)),max(c(lower_l,upper_l)))
				xx_l = c(curves[[index]]$x,rev(curves[[index]]$x)); yy_l = c(lower_l,rev(upper_l))
				plot(curves[[index]]$x, curves[[index]]$fit, ylim=yLim, ann=F, axes=F, type="l", col="gray30", lwd=1.0)
				polygon(xx_l, yy_l, col=rgb(100,100,100,100,maxColorValue=255), border=0)
				axis(side=1, pos=yLim[1], lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				axis(side=2, pos=min(curves[[index]]$x), lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				title(xlab=variableNames[i], cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
				title(ylab=paste0("s(",variableNames[i],")"), cex.lab=0.7, mgp=c(0.6,0,0), col.lab="gray30")
			}
		dev.new(width=6.5,height=3)
		par(mfrow=c(1,2), mar=c(3,3,1,1), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, col="gray30", bty="o")
		for (i in 1:length(selectedVariables))
			{				
				df = cbind(dfs[[1]][5:dim(dfs[[1]])[2]],dfs[[1]][,c("xCentroid","yCentroid")])
				for (j in 1:dim(df)[2])
					{
						if (colnames(df)[j] == selectedVariables[i])
							{
								df[,j] = seq(min(df[,j]),max(df[,j]),(max(df[,j])-min(df[,j]))/(dim(df)[1]-1))
							}	else	{
								df[,j] = median(df[,j])
							}
					}
				plot(df[,selectedVariables[i]], predict(gam, df), ann=F, axes=F, type="l", col="gray30", lwd=1.0)
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025,
					 col.tick="gray30", col.axis="gray30", col="gray30")
				title(xlab=variableNames[i], cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
				title(ylab=paste0("response"), cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
			}
	}

	# 4.8. Distributing community health workers

population = raster("WorldPop_pop_raster.tif")
communes@data$population = rep(0,dim(communes@data)[1])
communes@data$xCentroid = rep(0,dim(communes@data)[1])
communes@data$yCentroid = rep(0,dim(communes@data)[1])
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
		centroidCoordinates = coordinates(pol)
		communes@data[i,"population"] = exact_extract(population , sf::st_as_sfc(pol), fun='sum')
		communes@data[i,"xCentroid"] = centroidCoordinates[1,1]
		communes@data[i,"yCentroid"] = centroidCoordinates[1,2]
	}
communes@data$popDensityLog = log(communes@data$population/(communes@data$Shape_Area/10^6))
geoDistances = matrix(nrow=dim(communes@data)[1], ncol=dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		for (j in 1:dim(communes@data)[1])
			{
				d1 = communes@data[i,"xCentroid"]-communes@data[j,"xCentroid"]
				d2 = communes@data[i,"yCentroid"]-communes@data[j,"yCentroid"]
				geoDistances[i,j] = sqrt(((d1)^2)+((d2)^2))
			}
	}
data = read.csv("Google_Drive_N_Hens/Commune_cases_13-04.csv", head=T)
communes@data$casesLastFiveDays = rep(0,dim(communes@data)[1])
communes@data$cumulatedCases = rep(0,dim(communes@data)[1])
lastFiveDays = c("07-04-2020","08-04-2020","09-04-2020","10-04-2020","11-04-2020")
data[,"dateusedforstatistics"] = gsub("Jan","-01-",data[,"dateusedforstatistics"])
data[,"dateusedforstatistics"] = gsub("Feb","-02-",data[,"dateusedforstatistics"])
data[,"dateusedforstatistics"] = gsub("Mar","-03-",data[,"dateusedforstatistics"])
data[,"dateusedforstatistics"] = gsub("Apr","-04-",data[,"dateusedforstatistics"])
for (i in 1:dim(communes@data)[1])
	{
		nCases = 0
		for (j in 1:length(lastFiveDays))
			{
				indices = which((data[,"NIS5"]==communes@data[i,"NIS5"])&(data[,"dateusedforstatistics"]==lastFiveDays[j]))
				if (length(indices) > 0)
					{
						for (k in 1:length(indices)) nCases = nCases + 1
					}
			}
		communes@data[i,"casesLastFiveDays"] = nCases
		indices = which(data[,"NIS5"]==communes@data[i,"NIS5"])
		communes@data[i,"cumulatedCases"] = length(indices)
	}
communes@data$slidingWindow1 = rep(NA,dim(communes@data)[1])
communes@data$slidingWindow2 = rep(NA,dim(communes@data)[1])
communes@data$slidingWindow3 = rep(NA,dim(communes@data)[1])
for (i in 1:dim(communes@data)[1])
	{
		distances1 = geoDistances[i,]; distances2 = distances1[order(distances1)]
		distances3 = distances2[2:6]; indices = which(distances1%in%distances3)
		communes@data[i,"slidingWindow1"] = mean(c(communes@data[i,"casesLastFiveDays"],mean(communes@data[indices,"casesLastFiveDays"])))
		distances3 = distances2[2:11]; indices = which(distances1%in%distances3)
		communes@data[i,"slidingWindow2"] = mean(c(communes@data[i,"casesLastFiveDays"],mean(communes@data[indices,"casesLastFiveDays"])))
	}
communes@data$healthWorker1 = round(2000*(communes@data$slidingWindow1/(sum(communes@data$slidingWindow1))))
communes@data$healthWorker2 = round(2000*(communes@data$slidingWindow2/(sum(communes@data$slidingWindow2))))
communes@data$healthWorker1[communes@data$healthWorker1==0] = 1
communes@data$healthWorker2[communes@data$healthWorker2==0] = 1
provinces = unique(communes@data[,"TX_PROV_DE"])
provincesData = matrix(nrow=length(provinces), ncol=2)
colnames(provincesData) = c("population","casesLastFiveDays")
row.names(provincesData) = provinces
for (i in 1:length(provinces))
	{
		indices = which(communes@data[,"TX_PROV_DE"]==provinces[i])
		provincesData[i,1] = sum(communes@data[indices,"population"])
		provincesData[i,2] = sum(communes@data[indices,"casesLastFiveDays"])
	}
for (i in 1:dim(communes@data)[1])
	{
		index = which(provinces==communes@data[i,"TX_PROV_DE"])
		proportion = communes@data[i,"population"]/provincesData[index,"population"]
		communes@data[i,"slidingWindow3"] = round(provincesData[index,"casesLastFiveDays"]*proportion)
	}
communes@data$healthWorker1 = round(2000*(communes@data$slidingWindow1/(sum(communes@data$slidingWindow1))))
communes@data$healthWorker2 = round(2000*(communes@data$slidingWindow2/(sum(communes@data$slidingWindow2))))
communes@data$healthWorker3 = round(2000*(communes@data$slidingWindow3/(sum(communes@data$slidingWindow3))))
communes@data$healthWorker1[communes@data$healthWorker1==0] = 1
communes@data$healthWorker2[communes@data$healthWorker2==0] = 1
communes@data$healthWorker3[communes@data$healthWorker3==0] = 1
if (showingPlots)
	{
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (h in 1:3)
			{
				variableNames = c("# cases 07-11/04/2020",paste0("# cases (smoothing ",h,")"),"Community health workers")
				# variables = c("cumulatedCases","cumulatedCases","popDensityLog"); colourScales = list()
				variables = c("casesLastFiveDays",paste0("slidingWindow",h),paste0("healthWorker",h)); colourScales = list()
				colourScales[[1]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101]
				colourScales[[2]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101]
				colourScales[[3]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(121)[1:101]
				for (i in 1:length(variables))
					{
						if (i < 3)
							{
								minV = min(communes@data[,c(variables[1],variables[2])])
								maxV = max(communes@data[,c(variables[1],variables[2])])
								legendRast = raster(as.matrix(c(0,maxV))); # maxV = 100
							}	else	{
								minV = min(communes@data[,c(variables[3])])
								maxV = max(communes@data[,c(variables[3])])
								legendRast = raster(as.matrix(c(0,maxV))); # maxV = 20
							}
						values = communes@data[,variables[i]]; values[values[]>maxV] = maxV
						cols = colourScales[[i]][(((values -minV)/(maxV-minV))*100)+1]
						plot(communes_light, border="gray30", lwd=0.2, col=cols)
						mtext(paste0(variableNames[i]), cex=0.55, col="gray30", at=93000, line=-12.4)
						plot(legendRast, legend.only=T, col=colourScales[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 				 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
			 				 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.2,0)))
					}
			}
	}
if (writingFiles)
	{
		df = communes@data[,c("NIS5","population","casesLastFiveDays","slidingWindow1","slidingWindow2","slidingWindow3","healthWorker1","healthWorker2","healthWorker3")]
		colnames(df) = c("NIS","population","cases_07-11-20","smoothing_5_communes","smoothing_10_communes","smoothing_by_pronvinces",
						 "health_workers_smoothing_5_communes","health_workers_smoothing_10_communes","health_workers_smoothing_by_provinces")
		df$population = round(df$population); write.csv(df, "Community_health_workers.csv", row.names=F, quote=F)
	}

# 5. Analyses of hospital catchment areas

communes1 = shapefile("Shapefile_communes/Shapefile_NIS5_codes.shp")
communes2 = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
communes2 = spTransform(communes2, proj4string(communes1))
equivalence = read.csv("Shapefile_communes/Postal_codes_vs_NIS.csv", header=T)
catchmentAreas = shapefile("Hosp_catchment_areas/Hospital_catchment_areas_080420.shp")
catchmentAreas@data[,"X_ID"] = gsub(" ","",catchmentAreas@data[,"X_ID"])
catchmentAreas@data$area = as.numeric(catchmentAreas@data$area)

	# 5.1. Establishing the link between catchment areas and communes

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
				pol1 = sps; proj4string(pol1) = crs(catchmentAreas_pop)
				population = population + exact_extract(population_WP, sf::st_as_sfc(pol1), fun='sum')
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

	# 5.2. Extracting and assigning covariate values to each commune

data = read.csv("Data_Sciensano_1904/COVID19BE_CASES_MUNI_CUM.csv")
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

	# 5.3. 	Extracting and assigning covariate values to each catchment area

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
				pol1 = sps; proj4string(pol1) = crs(catchmentAreas_pop)
				population = population + exact_extract(population_WP, sf::st_as_sfc(pol1), fun='sum')
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
pm10 = raster("PM10_anmean_2017.asc")
pm25 = raster("PM25_anmean_2017.asc")
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
		pol = sps; proj4string(pol) = communes1@proj4string
		catchmentAreas@data[i,"pm10"] = exact_extract(pm10, sf::st_as_sfc(pol), fun='mean')
		catchmentAreas@data[i,"pm25"] = exact_extract(pm25, sf::st_as_sfc(pol), fun='mean')
		catchmentAreas@data[i,"ratioParkPop"] = exact_extract(pm25, sf::st_as_sfc(pol), fun='sum')/catchmentAreas@data[i,"population"]
	}
if (!file.exists("Hosp_catchment_areas/Proportion_of_urban_areas_080420.csv"))
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
		write.csv(propUrbanArea, "Hosp_catchment_areas/Proportion_of_urban_areas_080420.csv", row.names=F, quote=F)
	}
catchmentAreas@data$propUrbanArea = read.csv("Hosp_catchment_areas/Proportion_of_urban_areas_080420.csv")[,3]

	# 5.4. 	Extracting mobility data from Dalberg (mobile phone data)

if (!file.exists("Data_mobility_Dahlberg/Mobility_data_commune_matrix_20200416.csv"))
	{
		data = read.csv("Data_mobility_Dahlberg/Mobility_data_trips_per_capita_20200416.csv")
		days1 = c(paste0("2020-03-0",c(1:9)),paste0("2020-03-",c(10:31)),paste0("2020-04-0",c(1:9)),paste0("2020-04-",c(10:12)))
		days2 = gsub("-","",days1); days3 = ymd(days1)
		mobilityIndexCommunes = matrix(nrow=dim(communes2@data)[1], ncol=length(days1))
		row.names(mobilityIndexCommunes) = communes2@data[,"nouveau_PO"]; colnames(mobilityIndexCommunes) = days1
		for (i in 1:dim(communes2@data)[1])
			{
				for (j in 1:length(days2))
					{
						indices = which((data[,"postalcode"]==communes2@data[i,"nouveau_PO"])&(data[,"date"]==days2[j]))
						if (length(indices) > 0) mobilityIndexCommunes[i,j] = mean(data[indices,"Proximus_Trips_Per_Capita"])
					}
			}
		mobilityIndexCatchments = matrix(nrow=dim(catchmentAreas@data)[1], ncol=length(days1))
		colnames(mobilityIndexCatchments) = paste0("mobilityIndex_",days1)
		for (i in 1:dim(catchmentAreas@data)[1])
			{
				for (j in 1:length(days2))
					{
						mobilityIndexCatchments[i,j] = sum(mobilityIndexCommunes[,j]*populations2[,i],na.rm=T)/catchmentAreas@data[i,"population"]
					}
			}
		write.table(mobilityIndexCatchments, "Data_mobility_Dahlberg/Mobility_data_catchement_areas_20200416.csv", quote=F, row.names=F, sep=",")
	}
mobilityIndexCatchments = read.csv("Data_mobility_Dahlberg/Mobility_data_catchement_areas_20200416.csv")
colnames(mobilityIndexCatchments) = gsub("\\.","-",colnames(mobilityIndexCatchments))
catchmentAreas@data = cbind(catchmentAreas@data, mobilityIndexCatchments)

	# 5.5. Cumputing doubling times for hospitalisations and ICU

data = read.csv("Google_Drive_N_Hens/Data_Hospit_2_20-04.csv", head=T, sep=";")
days1 = c(paste0("2020-03-0",c(1:9)),paste0("2020-03-",c(10:31)),paste0("2020-04-0",c(1:9)),paste0("2020-04-",c(10:20)))
firstDay = ymd("2020-01-30"); days2 = ymd(days1); days3 = as.numeric(days2-firstDay); D = 7 # time interval
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
		temp1 = data[lines,c("days","NewPatientsNotReferred","Confirmed.patients.in.hospital","Confirmed.patients.in.ICU")]
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
						newHospitaliCases[i,j] = temp[index2,"NewPatientsNotReferred"]
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

	# 5.6. Saving and plotting the variables assigned to each area

if (writingFiles)
	{
		df = catchmentAreas@data
		indices = which(grepl("cumulated",colnames(df)))
		for (i in 1:length(indices)) df[which(is.na(df[,indices[i]])),indices[i]] = 0
		indices = which(grepl("incidence",colnames(df)))
		for (i in 1:length(indices)) df[which(is.na(df[,indices[i]])),indices[i]] = 0
		colnames(df) = gsub("cumulatedH","h",colnames(df)); colnames(df) = gsub("cumulatedI","I",colnames(df))
		write.csv(df, "All_figures_&_outputs/Covariables_catchment_areas_200420.csv", row.names=F, quote=F)
	}
if (showingPlots)
	{
		periods = c("16/03-22/03/2020","23/03-29/03/2020","30/03-05/04/2020")
		selectedDays1 = ymd(c("2020-03-22","2020-03-29","2020-04-05"))
		variables = c("incidenceHosCases_2020-03-22","incidenceHosCases_2020-03-29","incidenceHosCases_2020-04-05",
					  "doublingTHosCases_2020-03-22","doublingTHosCases_2020-03-29","doublingTHosCases_2020-04-05",
					  "incidenceICUCases_2020-03-22","incidenceICUCases_2020-03-29","incidenceICUCases_2020-04-05",
					  "doublingTICUCases_2020-03-22","doublingTICUCases_2020-03-29","doublingTICUCases_2020-04-05")
		variableNames1 = c("Hosp. incidence 22/03/20","Hosp. incidence 29/03/20","Hosp. incidence 05/04/20",
						  "Hosp. doubling time","Hosp. doubling time","Hosp. doubling time",
						  "ICU incidence 22/03/20","ICU incidence 29/03/20","ICU incidence 05/04/20",
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
		variables = c("growthRHosCases_2020-03-22","growthRHosCases_2020-03-29","growthRHosCases_2020-04-05",
					  "doublingTHosCases_2020-03-22","doublingTHosCases_2020-03-29","doublingTHosCases_2020-04-05",
					  "growthRICUCases_2020-03-22","growthRICUCases_2020-03-29","growthRICUCases_2020-04-05",
					  "doublingTICUCases_2020-03-22","doublingTICUCases_2020-03-29","doublingTICUCases_2020-04-05")
		variableNames1 = c("Hosp. growth rate 22/03/20","Hosp. growth rate 29/03/20","Hosp. growth rate 05/04/20",
						  "Hosp. doubling time","Hosp. doubling time","Hosp. doubling time",
						  "ICU growthR 22/03/20","ICU growthR 29/03/20","ICU growthR 05/04/20",
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
		variables = c("mobilityIndex_2020-03-08","mobilityIndex_2020-03-15","mobilityIndex_2020-03-22",
					  "ratioParkPop","medianIncome","popDensityLog","sectorP","sectorS","sectorT","moreThan65","bedsInMRs","pm10")
		variableNames = c("Mobility index 22/03/20","Mobility index 29/03/20","Mobility index 05/04/20",
						  "Ratio park area/population","Median declared income (€)","Population density (log)",
						  "% in primary sector","% in secundary sector","% in tertiary sector",
						  ">= 65 years (proportion)","# beds in rest houses","PM 1.0 emission")
		colourScales = list()
		colourScales[[1]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[2]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[3]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101])
		colourScales[[4]] = c(colorRampPalette(brewer.pal(9,"YlGn"))(151)[1:101])
		colourScales[[5]] = c(colorRampPalette(brewer.pal(9,"RdPu"))(151)[1:101])
		colourScales[[6]] = c(colorRampPalette(brewer.pal(9,"BuPu"))(151)[1:101])
		colourScales[[7]] = c(colorRampPalette(brewer.pal(9,"Greens"))(151)[1:101])
		colourScales[[8]] = c(colorRampPalette(brewer.pal(9,"Oranges"))(151)[1:101])
		colourScales[[9]] = c(colorRampPalette(brewer.pal(9,"Blues"))(151)[1:101])
		colourScales[[10]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
		colourScales[[11]] = c(colorRampPalette(brewer.pal(9,"PuBuGn"))(151)[1:101])
		colourScales[[12]] = c(colorRampPalette(brewer.pal(9,"YlOrBr"))(151)[1:101])
		dev.new(width=7,height=8); par(mfrow=c(4,3), mar=c(0,0,0,0), oma=c(2,2,1,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
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
				plot(catchmentAreas, border="gray30", col=cols, lwd=0.1)
				mtext(variableNames[i], cex=0.54, col="gray30", at=92000, line=-12.4)
				plot(legendRast, legend.only=T, col=legendCols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.5,0.10,0.12),
			 		 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.7, lwd=0,
					 lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.13,0)))
			}
	}
	
	# 5.7. Relation between mobility and hospitalisation

reloadingFile = FALSE
if (reloadingFile)
	{
		catchmentAreas = shapefile("Hosp_catchment_areas/Hospital_catchment_areas_080420.shp")
		catchmentAreas@data = read.csv("All_figures_&_outputs/Covariables_catchment_areas_090420.csv")
		colnames(catchmentAreas@data) = gsub("\\.","-",colnames(catchmentAreas@data))
		colnames(catchmentAreas@data) = gsub("hosCases_","cumulatedHosCases_",colnames(catchmentAreas@data))
	}
mobilityIndices = catchmentAreas@data[,which(grepl("mobilityIndex",colnames(catchmentAreas@data)))]
newHospitaCases = catchmentAreas@data[,which(grepl("newHospitaliCases",colnames(catchmentAreas@data)))]
growthRateHosps = catchmentAreas@data[,which(grepl("growthRHosCases",colnames(catchmentAreas@data)))]
days1 = c(paste0("2020-03-",c(15:31)),paste0("2020-04-0",c(1:9)),paste0("2020-04-",c(10:12)))
Ds = c(1:20); R2s = matrix(nrow=dim(growthRateHosps)[1], ncol=length(Ds)); cors = matrix(nrow=dim(growthRateHosps)[1], ncol=length(Ds))
for (i in 1:length(Ds))
	{
		for (j in 1:dim(growthRateHosps)[1])
			{
				days2 = as.character(ymd(days1)-Ds[i]); days3 = as.character(ymd(days1)-1)
				indices = which(days2%in%gsub("mobilityIndex_","",colnames(mobilityIndices)))
				x = t(mobilityIndices[j,paste0("mobilityIndex_",days2[indices])])
				y = t(newHospitaCases[j,paste0("newHospitaliCases_",days1[indices])])
				y = t(growthRateHosps[j,paste0("growthRHosCases_",days1[indices])])
				y[is.infinite(y)] = NA; x = x[!is.na(y)]; y = y[!is.na(y)]; # plot(x,y) }
				if (length(y) > 0)
					{
						lr = lm("y ~ x"); R2s[j,i] = summary(lr)$r.square; cors[j,i] = cor(x,y)
					}
			}
	}
library(ggridges); library(ggplot2); library(viridis); library(hrbrthemes)
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
	   scale_fill_viridis(name="Temp. [F]", option="C") + theme_ipsum() +
	   theme(legend.position="none", panel.spacing=unit(0.1,"lines"), strip.text.x=element_text(size=7))

newHospitalisations = c(); days = c()
for (i in 1:length(days1))
	{
		newHospitalisations = c(newHospitalisations, newHospitaCases[,paste0("newHospitaliCases_",days1[i])])
		days = c(days, rep(i,dim(newHospitaCases)[1]))
	}
df = as_tibble(data.frame(cbind(as.factor(days),newHospitalisations)))
ggplot(df, aes(x=newHospitalisations, y=as.factor(days), fill=..x..)) +
	   geom_density_ridges_gradient(scale=3, rel_min_height=0.01) +
	   scale_fill_viridis(name="Temp. [F]", option="C") + theme_ipsum() +
	   theme(legend.position="none", panel.spacing=unit(0.1,"lines"), strip.text.x=element_text(size=7))

mobilityIndices = catchmentAreas@data[,which(grepl("mobilityIndex",colnames(catchmentAreas@data)))]
lagAndAveragedM = catchmentAreas@data[,which(grepl("incidenceHosCases",colnames(catchmentAreas@data)))]
lagAndAveragedM[,] = NA; colnames(lagAndAveragedM) = gsub("incidenceHosCases","lagAndAveragedM",colnames(lagAndAveragedM))
days1 = ymd(gsub("lagAndAveragedM_","",colnames(lagAndAveragedM)))
index = which(as.character(days1)=="2020-03-14")
for (i in index:length(days1))
	{
		days2 = as.character(days1[i]-c(6:15))
		for (j in 1:dim(lagAndAveragedM)[1])
			{
				indices = which(days2%in%gsub("mobilityIndex_","",colnames(mobilityIndices)))
				lagAndAveragedM[j,i] = mean(as.numeric(mobilityIndices[j,paste0("mobilityIndex_",days2[indices])]), na.rm=T)
			}
	}
lagAndAveragedM = lagAndAveragedM[,index:dim(lagAndAveragedM)[2]]
catchmentAreas@data = cbind(catchmentAreas@data, lagAndAveragedM)

	# 5.8. Performing and plotting the first axes of an exploratory PCA

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

	# 5.9. Assessing spatial autocorrelation with the Moran's I test

variableNames = c("xCentroid","yCentroid","popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
				  "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea",
				  "lagAndAveragedM_2020-03-14","lagAndAveragedM_2020-03-21","lagAndAveragedM_2020-03-28","lagAndAveragedM_2020-04-05",
				  "incidenceHosCases_2020-03-13","incidenceHosCases_2020-03-20","incidenceHosCases_2020-03-27","incidenceHosCases_2020-04-04",
				  "incidenceHosCases_2020-03-14","incidenceHosCases_2020-03-21","incidenceHosCases_2020-03-28","incidenceHosCases_2020-04-05",
				  "growthRHosCases_2020-03-14","growthRHosCases_2020-03-21","growthRHosCases_2020-03-28","growthRHosCases_2020-04-05")
df = catchmentAreas@data[,variableNames]; colnames(df) = gsub("-","",colnames(df))
responseVariables = c("incidenceHosCases_20200314","incidenceHosCases_20200321","incidenceHosCases_20200328","incidenceHosCases_20200405")
responseVariables = c("growthRHosCases_20200321","growthRHosCases_20200328","growthRHosCases_20200405")
for (i in 1:length(responseVariables))
	{
		values = df[,responseVariables[i]]
		indices = which((!is.na(values))&(!is.infinite(values)))
		geoDists = as.matrix(dist(df[indices,c("xCentroid","yCentroid")]))
		weights = 1/geoDists; diag(weights) = 0
		print(Moran.I(values[indices],weights)$p.value)
	}

	# 5.10. Univariate (LR) followed by multivariate regression (GLM) analyses

selectedVariables = list()
for (i in 1:length(responseVariables))
	{
		if (grepl("incidence",responseVariables[i]))
			{
				predictors = c("popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
			   				   "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea",
			   				   "lagAndAveragedM_20200314","lagAndAveragedM_20200321","lagAndAveragedM_20200328","lagAndAveragedM_20200405")
				day = as.character(ymd(unlist(strsplit(responseVariables[i],"_"))[2])-1)
				predictors = c(predictors, paste0("incidenceHosCases_",gsub("-","",day)))
			}
		if (grepl("growthR",responseVariables[i]))
			{
				predictors = c("popDensity","popDensityLog","medianAge","moreThan65","maisonsDeRepos","bedsInMRs",
			  				   "medianIncome","sectorP","sectorS","sectorT","cases","pm10","pm25","ratioParkPop","propUrbanArea",
							   "lagAndAveragedM_20200314","lagAndAveragedM_20200321","lagAndAveragedM_20200328","lagAndAveragedM_20200405")
			}
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

	# 5.11. GAM (generalised additive model) analyses

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
		if (grepl("incidence",responseVariables[i]))
			{
				if (i == 1) gam = gam(responseVariable ~ s(sectorT) + s(propUrbanArea) +
						  			  s(incidenceHosCases_20200313) + s(lagAndAveragedM_20200314), data=tmp, method="REML")
				if (i == 2) gam = gam(responseVariable ~ s(sectorT) + s(propUrbanArea) +
						  			  s(incidenceHosCases_20200320) + s(lagAndAveragedM_20200321), data=tmp, method="REML")
				if (i == 3) gam = gam(responseVariable ~ s(sectorT) + s(propUrbanArea) +
						  			  s(incidenceHosCases_20200327) + s(lagAndAveragedM_20200328), data=tmp, method="REML")
				if (i == 4) gam = gam(responseVariable ~ s(sectorT) + s(propUrbanArea) +
						  			  s(incidenceHosCases_20200404) + s(lagAndAveragedM_20200405), data=tmp, method="REML")
			}
		if (grepl("growthR",responseVariables[i]))
			{
				if (i == 1) gam = gam(responseVariable ~ s(propUrbanArea) + s(lagAndAveragedM_20200321) +
						  			  s(xCentroid,yCentroid), data=tmp, method="REML")
				if (i == 2) gam = gam(responseVariable ~ s(propUrbanArea) + s(lagAndAveragedM_20200328) +
									  s(xCentroid,yCentroid), data=tmp, method="REML")
				if (i == 3) gam = gam(responseVariable ~ s(propUrbanArea) + s(lagAndAveragedM_20200405) +
									  s(xCentroid,yCentroid), data=tmp, method="REML")
			}
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
		selectedVariables = c("propUrbanArea")
		variableNames = c("proportion of urban area")
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
				tmp = df[,c("sectorT","propUrbanArea","lagAndAveragedM_20200321","lagAndAveragedM_20200328","lagAndAveragedM_20200405")]
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

	# 5.12. Classic correlation analyses

variableNames = c("incidenceHosCases_2020-03-14","incidenceHosCases_2020-03-21","incidenceHosCases_2020-03-28","incidenceHosCases_2020-04-05",
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
	
# 6. Phylogenetic and phylogeographic analyses

	# 6.1. Preparing the input files for the discrete phylogeographic analyses 

tree = read.tree("Phylogenetic_analyses/Nextstrain_200420.tree")
data = read.csv("Phylogenetic_analyses/Nextstrain_200420.csv", sep=";")
	# N.B.: the ".tsv" file has first to be opened in Excel and then exported as comma ".csv"
labs = unique(data[which(data[,"Country"]=="Belgium"),"Submitting.Lab"])
if (showingPlots)
	{
		dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
	}
txt = c(); tab1 = c(); tab2 = c(); tab3 = c()
selectedCountries = c("Belgium","Italy","Spain","France","England","Germany","Holland","Luxembourg")
for (i in 1:length(tree$tip.label))
	{
		index = which(data[,1]==tree$tip.label[i])
		date = as.character(data[index,"Collection.Data"])
		if (date != "")
			{
				txt = c(txt, paste0(">",tree$tip.label[i]),"NNNN")
				location = unlist(strsplit(tree$tip.label[i],"\\/"))[1]
				tab1 = rbind(tab1, cbind(tree$tip.label[i],location,date))
				region = gsub(" ","",as.character(data[index,"Region"]))
				if (!location%in%selectedCountries) location = region
				tab2 = rbind(tab2, cbind(tree$tip.label[i],location,date))
				if (location != "Belgium") location = "other"
				tab3 = rbind(tab3, cbind(tree$tip.label[i],location,date))
			}
	}
write(txt, "Phylogenetic_analyses/Nextstrain_200420.fasta")
colnames(tab1) = c("trait","location","collection_date")
colnames(tab2) = c("trait","location","collection_date")
colnames(tab3) = c("trait","location","collection_date")
write.table(tab1, "Phylogenetic_analyses/Nextstrain_200420_1.txt", row.names=F, quote=F, sep="\t")
write.table(tab2, "Phylogenetic_analyses/Nextstrain_200420_2.txt", row.names=F, quote=F, sep="\t")
write.table(tab3, "Phylogenetic_analyses/Nextstrain_200420_3.txt", row.names=F, quote=F, sep="\t")

	# 6.2. Analysing the outputs of the preliminary discrete phylogeographic analysis 

computingHPDInterval = FALSE # N.B.: long analysis
if (computingHPDInterval)
	{
		trees = readAnnotatedNexus("Phylogenetic_analyses/Nextstrain_2004_3.trees")
		belgianBranches_list = rep(NA,length(trees))
		belgianIntroductions_list = rep(NA,length(trees))
		belgianTipBranches_list = rep(NA,length(trees))
		for (i in 1:length(trees))
			{
				belgianBranches = 0; belgianIntroductions = 0; belgianTipBranches = 0
				for (j in 1:dim(trees[[i]]$edge)[1])
					{
						if (tree$annotations[[j]]$location == "Belgium")
							{
								belgianBranches = belgianBranches + 1
								index = which(trees[[i]]$edge[,2]==trees[[i]]$edge[j,1])
								if (trees[[i]]$annotations[[index]]$location != "Belgium")
									{
										belgianIntroductions = belgianIntroductions + 1
									}
								if (!trees[[i]]$edge[j,2]%in%trees[[i]]$edge[,1])
									{
										belgianTipBranches = belgianTipBranches + 1
									}
							}
					}
				belgianBranches_list[i] = belgianBranches
				belgianIntroductions_list[i] = belgianIntroductions
				belgianTipBranches_list[i] = belgianTipBranches
			}
		quantiles = quantile(belgianIntroductions_list,probs=c(0.025,0.975))
		cat("A minimum number of ",median(belgianIntroductions_list)," lineage introductions (95% HPD interval = [",quantiles[1],"-",quantiles[2],"])",
			" identified from the global phylogenetic analysis of ",belgianTipBranches," SARS-CoV-2 sampled in Belgium (07-04-2020)",sep="")
		# A minimum number of 125 lineage introductions (95% HPD interval = [112-139]) identified from the global phylogenetic analysis of 253 SARS-CoV-2 sampled in Belgium 
	}
if (showingPlots)
	{
		tree = readAnnotatedNexus("Phylogenetic_analyses/Nextstrain_2004_3_MCC.tree")
		samplingDates = decimal_date(dmy(gsub("\\/","-",tab1[,"collection_date"]))); mostRecentSamplingYear = max(samplingDates)
		selectedDates = decimal_date(ymd(c("2019-12-19","2020-01-01","2020-01-15","2020-02-01","2020-02-15","2020-03-01","2020-03-15","2020-04-01","2020-04-15")))
		rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingYear-rootHeight
		selectedLabels = c("15-19","01-01","15-01","01-02","15-02","01-03","15-03","01-04","15-04")
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.3, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.3, col="gray30", lwd=0.5)
					}
				if (tree$annotations[[i]]$location == "Belgium")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Belgium")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		cols = rep("gray50",dim(tree$edge)[1]); lwds = rep(0.05,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		dev.new(width=11, height=8)
		plot(tree, show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$annotations[[i]]$location == "Belgium")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Belgium")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}	else	{
								if (!tree$edge[i,2]%in%tree$edge[,1])
									{	}
							}
					}
			}
		axis(lwd=0.2, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.65, mgp=c(0,0.1,-0.9), lwd.tick=0.2, 
			 col.lab="gray30", col="gray30", tck=-0.01, side=1)
	}

tree = readAnnotatedNexus("Phylogenetic_analyses/Nextstrain_2004_3_MCC.tree")
data1 = read.csv("Sequences_metadata/SARS-CoV-2_KULeuven_080420.csv")
data2 = read.csv("Sequences_metadata/SARS-CoV-2_ULiegeSeq_220420.csv")
communes = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
belgium = spTransform(raster::getData("GADM", country="BEL", level=0), crs(communes))
rast = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(communes)); rast[] = log(rast[])
belgianBranches = c(); belgianIntroductions = c()
belgianTipBranches = c(); sampledSequences = c()
for (i in 1:dim(tree$edge)[1])
	{
		if (tree$annotations[[i]]$location == "Belgium")
			{
				belgianBranches = c(belgianBranches,i)
				index = which(tree$edge[,2]==tree$edge[i,1])
				if (tree$annotations[[index]]$location != "Belgium")
					{
						belgianIntroductions = c(belgianIntroductions, i)
					}
				if (!tree$edge[i,2]%in%tree$edge[,1])
					{
						belgianTipBranches = c(belgianTipBranches, i)
						sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]])
					}
			}
	}
for (i in 1:length(belgianIntroductions))
	{
		if (i == 1) clusters1 = list()
		if (tree$edge[belgianIntroductions[i],2]%in%tree$edge[,1])
			{
				subtree = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				clusters1[[i]] = subtree$tip.label
			}	else		{
				clusters1[[i]] = tree$tip.label[tree$edge[belgianIntroductions[i],2]]	
			}
	}
if (!file.exists(paste0("Phylogenetic_analyses/Sampling_Belgium.csv")))
	{
		samplingData = matrix(nrow=length(sampledSequences), ncol=5)
		colnames(samplingData) = c("sequenceID","collectionDate","postCode","longitude","latitude")
		samplingData[,"sequenceID"] = sampledSequences
		for (i in 1:dim(samplingData)[1])
			{
				index = which(data[,"Strain"]==samplingData[i,"sequenceID"])
				date = dmy(gsub("\\/","-",data[index,"Collection.Data"]))
				samplingData[i,"collectionDate"] = decimal_date(date)
				ID = unlist(strsplit(samplingData[i,"sequenceID"],"\\/"))[2]
				index1 = which(grepl(ID,data1[,"sequence.name"]))
				if (length(index1) == 1)
					{
						samplingData[i,"postCode"] = data1[index1,"Post.code"]
						indices = which(communes@data[,"nouveau_PO"]==data1[index1,"Post.code"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
				index1 = which(grepl(ID,data2[,"gisaid.virus.name"]))
				if (length(index1) == 1)
					{
						samplingData[i,"postCode"] = as.character(data2[index1,"Postal.code"])
						indices = which(communes@data[,"nouveau_PO"]==data2[index1,"Postal.code"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
			}
		print(samplingData[which(is.na(samplingData[,"postCode"])),"sequenceID"])
		write.csv(samplingData, "Phylogenetic_analyses/Sampling_Belgium.csv", quote=F, row.names=F)
	}	
samplingData = read.csv("Phylogenetic_analyses/Sampling_Belgium.csv", head=T)
for (i in 1:length(belgianIntroductions))
	{
		tab = c()
		if (i == 1)
			{
				clusters2 = list(); centroids = list()
			}
		for (j in 1:length(clusters1[[i]]))
			{
				index = which(samplingData[,"sequenceID"]==clusters1[[i]][j])
				if (length(index) == 1)
					{
						line = cbind(as.numeric(samplingData[index,"collectionDate"]),as.numeric(samplingData[index,"longitude"]),as.numeric(samplingData[index,"latitude"]))
						row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
					}
			}
		colnames(tab) = c("collectionDate","longitude","latitude"); clusters2[[i]] = tab
		centroids[[i]] = cbind(mean(tab[!is.na(tab[,1]),1]), mean(tab[!is.na(tab[,2]),2]))
	}
if (showingPlots)
	{
		dev.new(width=7, height=6); par(oma=c(0,0,0,0), mar=c(0.5,4,1,0), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"BuPu"))(161)[1:101])
		cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
		plot(rast, col=cols, axes=F, ann=F, box=F, legend=F)
		plot(belgium, border="gray40", lwd=0.4, add=T)
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						if (length(clusters1[[i]]) > 1)
							{
								for (j in 1:dim(clusters2[[i]])[1])
									{
										if (!is.na(clusters2[[i]][j,1]))
											{
												segments(centroids[[i]][,1],centroids[[i]][,2],clusters2[[i]][j,1],clusters2[[i]][j,2], lwd=0.5, col="gray30")	
											}
									}
							}
					}
			}
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						for (j in 1:dim(clusters2[[i]])[1])
							{
								points(clusters2[[i]][,1], clusters2[[i]][,2], pch=16, cex=0.8, col="chartreuse3")
								points(clusters2[[i]][,1], clusters2[[i]][,2], pch=1, cex=0.8, col="gray30", lwd=0.2)
							}
					}
			}
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						if (length(clusters1[[i]]) > 1)
							{
								points(centroids[[i]][,1], centroids[[i]][,2], pch=16, cex=0.6, col="red")
								points(centroids[[i]][,1], centroids[[i]][,2], pch=1, cex=0.6, col="gray30", lwd=0.2)
							}
					}
			}
		legendRast = raster(as.matrix(c(min(rast[],na.rm=T),max(rast[],na.rm=T))))
		mtext("Human population (log-transformed)", col="gray30", cex=0.7, line=-23, at=78000)
		plot(legendRast, legend.only=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.141,0.409,0.18,0.19),
			 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.55, lwd=0,
			 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,-0.05,0)))
	}

phylogeneticDistances = as.matrix(distTips(tree, tips=samplingData[,"sequenceID"], method="patristic"))
geographicDistances = matrix(0, nrow=dim(samplingData)[1], ncol=dim(samplingData)[1])
for (i in 2:dim(geographicDistances)[1])
	{
		for (j in 1:(i-1))
			{
				x1 = as.numeric(samplingData[i,"longitude"]); y1 = as.numeric(samplingData[i,"latitude"])
				x2 = as.numeric(samplingData[j,"longitude"]); y2 = as.numeric(samplingData[j,"latitude"])
				geographicDistances[i,j] = (sqrt(((x2-x1)^2)+((y2-y1)^2)))/1000
			}
	}
ibds = matrix(nrow=dim(samplingData)[1], ncol=2)
geos = matrix(nrow=dim(samplingData)[1], ncol=2)
for (i in 1:dim(samplingData)[1])
	{
		for (j in 1:length(clusters1))
			{
				if (samplingData[i,"sequenceID"]%in%clusters1[[j]])
					{
						if (length(clusters1[[j]]) > 1)
							{
								seqIDs = clusters1[[j]][clusters1[[j]][]!=samplingData[i,"sequenceID"]]
								indices1 = which((samplingData[,"sequenceID"]%in%seqIDs)&(samplingData[,"collectionDate"]<samplingData[i,"collectionDate"]))
								if (length(indices1) > 0)
									{
										phylogeneticDis = phylogeneticDistances[i,indices1]
										geographicDis = geographicDistances[i,indices1]
										indices2 = which(phylogeneticDis==min(phylogeneticDis)); vS1 = c(); vS2 = c()
										for (k in 1:length(indices2))
											{
												vS1 = c(vS1, geographicDis[indices2[k]]/phylogeneticDis[indices2[k]])
												vS2 = c(vS2, geographicDis[indices2[k]])
											}
										ibds[i,1] = samplingData[i,"collectionDate"]; ibds[i,2] = mean(vS1)
										geos[i,1] = samplingData[i,"collectionDate"]; geos[i,2] = mean(vS2)
									}
							}
					}
			}
	}
ibds = ibds[which(!is.na(ibds[,2])),]; buffer = cbind(as.numeric(ibds[,1]), as.numeric(ibds[,2])); ibds = buffer[order(buffer[,1]),]
geos = geos[which(!is.na(geos[,2])),]; buffer = cbind(as.numeric(geos[,1]), as.numeric(geos[,2])); geos = buffer[order(buffer[,1]),]

	# 6.4. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)

coordinates = matrix(nrow=dim(tab3)[1], ncol=2)
colnames(coordinates) = c("latitude","longitude")
tab4 = cbind(tab3, coordinates)
for (i in 1:dim(tab4)[1])
	{
		index = which(samplingData[,"sequenceID"]==tab4[i,"trait"])
		if (length(index) > 0)
			{
				tab4[i,"latitude"] = samplingData[index,"latitude"]
				tab4[i,"longitude"] = samplingData[index,"longitude"]
			}
	}
write.table(tab4, "Phylogenetic_analyses/Nextstrain_200420_4.txt", row.names=F, quote=F, sep="\t")

analyses = c()
template = scan("Phylogenetic_analyses/Nextstrain_RRW_temp.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3))
			{
				xml = gsub("TEMPLATE", paste0("Clade_",i), template); analyses = c(analyses, paste0("Clade_",i))
				tre = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				tips_to_drop = tre$tip.label[which(!tre$tip.label%in%row.names(clusters2[[i]]))]
				if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
				write.tree(tre, paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".tre"))
				tre = scan(paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".tre"), what="", sep="\n", quiet=T)
				sink(file=paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".xml"))
				for (j in 1:length(xml))
					{
						cat(xml[j]); cat("\n")
						if (xml[j]=="\t<taxa id=\"taxa\">")
							{
								for (k in 1:dim(clusters2[[i]])[1])
									{
										cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[i]])[k],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",clusters2[[i]][k,"collectionDate"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"latitude\">\n")
										cat(paste0("\t\t\t\t",clusters2[[i]][k,"latitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"longitude\">\n")
										cat(paste0("\t\t\t\t",clusters2[[i]][k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"coordinates\">\n")
										cat(paste0("\t\t\t\t",clusters2[[i]][k,"latitude"]," ",clusters2[[i]][k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
							}
						if (xml[j]=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")
							{
								for (k in 1:dim(clusters2[[i]])[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[i]])[k],"\"/>","\n"))
										cat("\t\t\tNNNN\n")
										cat("\t\t</sequence>\n")
									}
							}
						if (xml[j]=="\t<newick id=\"startingTree\">")
							{
								cat(paste0("\t\t",tre,"\n"))
							}
					}
				sink(NULL)
			}
	}

	# 6.5. Running BEAST and building the maximum clade consensus (MCC) tree

source("MCC_tree_extraction.r")
sink(file=paste0("Phylogenetic_analyses/Phylogeographic_runs/Analyses.sh"))
for (i in 1:length(analyses))
	{
		cat(paste0("java -jar beast_1104.jar -overwrite ",analyses[i],".xml\n"))
	}
sink(NULL)
wd = getwd()
setwd(paste0(wd,"/Phylogenetic_analyses/Phylogeographic_runs/"))
system("bash Analyses.sh", ignore.stdout=T, ignore.stderr=F)
for (i in 1:length(analyses))
	{
		system(paste0("BEAST_1104/bin/treeannotator -burninTrees 101 -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
	}

	# 6.6. Extracting the spatio-temporal information embedded in the MCC and posterior trees

for (i in 1:length(analyses))
	{
		index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
		mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
		mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree"))
		mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
	}
nberOfTreesToSample = 1000; burnIn = 101; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
for (i in 1:length(analyses))
	{
		index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
		mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
		allTrees = scan(file=paste0(analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F); localTreesDirectory = paste0(analyses[i],"_ext")
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}
for (i in 1:length(analyses))
	{
		tab = read.csv(paste0(analyses[i],".csv"), head=T)
		if (i == 1)
			{
				all = tab
			}	else	{
				maxNodeID = max(all[,c("node1","node2")])
				tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
				all = rbind(all, tab)
			}
	}
write.csv(all, "All_clades.csv", row.names=F, quote=F)
dir.create(file.path("All_clades_ext1"), showWarnings=F)
dir.create(file.path("All_clades_ext2"), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfTreesToSample)
	{
		for (j in 1:length(analyses))
			{
				tab = read.csv(paste0(analyses[j],"_ext/TreeExtractions_",i,".csv"), head=T)
				if (j == 1)
					{
						all = tab
					}	else	{
						maxNodeID = max(all[,c("node1","node2")])
						tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
						all = rbind(all, tab)
					}
			}
		write.csv(all, paste0("All_clades_ext1/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		temp1 = all[,c("startLon","startLat")]; temp2 = all[,c("endLon","endLat")]
		coordinates(temp1) = ~ startLon + startLat; crs(temp1) = crs(communes)
		coordinates(temp2) = ~ endLon + endLat; crs(temp2) = crs(communes)
		temp1 = spTransform(temp1, CRS("+init=epsg:4326"))@coords
		temp2 = spTransform(temp2, CRS("+init=epsg:4326"))@coords
		all[,c("startLon","startLat")] = temp1; all[,c("endLon","endLat")] = temp2
		write.csv(all, paste0("All_clades_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
setwd(wd)

	# 6.7. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)

localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext1")
percentage = 80; prob = percentage/100; precision = 1/(365/3.5)
mcc = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T); startDatum = min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
if (showingPlots)
	{
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[21:121])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		firstTimePeriod = TRUE; secondTimePeriod = FALSE
		firstTimePeriod = FALSE; secondTimePeriod = TRUE
		firstTimePeriod = FALSE; secondTimePeriod = FALSE
		dev.new(width=7.3, height=6); par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		plot(belgium, border="gray30", col="gray95", lwd=0.4); croppingPolygons = TRUE
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(polygons))
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(belgium)
						if (croppingPolygons == TRUE) pol = crop(pol, belgium)
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		selectedBranches = 1:dim(mcc)[1]
		if (firstTimePeriod) selectedBranches = which(mcc[,"endYear"]<decimal_date(ymd("2020-03-18")))
		if (secondTimePeriod) selectedBranches = which(mcc[,"startYear"]>decimal_date(ymd("2020-03-18")))
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (i in 1:length(clusters2))
			{
				if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) == 2)
							{
								indices = which(!is.na(clusters2[[i]][,"longitude"]))
								if (firstTimePeriod) indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]<decimal_date(ymd("2020-03-18"))))
								if (secondTimePeriod) indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]>decimal_date(ymd("2020-03-18"))))
								if (length(indices) == 2)
									{
										curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
													cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
													arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray30", arr.col=NA, arr.pos=F,
													curve=0.1, dr=NA, endhead=F)
									}
							}
					}
			}
		for (i in 1:length(clusters2))
			{
				if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
					{
						for (j in 1:dim(clusters2[[i]])[1])
							{
								if (!is.na(clusters2[[i]][j,"longitude"]))
									{
										plotTheNode = TRUE
										if ((firstTimePeriod==TRUE)&(clusters2[[i]][j,"collectionDate"]>decimal_date(ymd("2020-03-18")))) plotTheNode = FALSE
										if ((secondTimePeriod==TRUE)&(clusters2[[i]][j,"collectionDate"]<decimal_date(ymd("2020-03-18")))) plotTheNode = FALSE
										if (plotTheNode)
											{
												index = (((clusters2[[i]][j,"collectionDate"]-minYear)/(maxYear-minYear))*100)+1
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=16, col=colourScale[index], cex=0.8)
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=1, col="gray30", cex=0.8, lwd=0.4)
											}
									}
							}
					}
			}
		for (i in rev(selectedBranches))
			{
				if (!mcc[i,"node1"]%in%mcc[,"node2"])
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=0.8)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=0.8, lwd=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.8, lwd=0.4)
			}
		selectedDates = decimal_date(ymd(c("2020-03-03","2020-03-18","2020-04-03")))
		selectedLabels = c("03-03-2020","18-03-2020","03-04-2020")
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.1,0.5,0.100,0.112),
			 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.00,0),
		     at=selectedDates, labels=selectedLabels))
	}

	# 6.8. Estimating and plotting dispersal statistics associated with lineages
	
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext2")
timSlices = 50; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 5; slidingWindow = 1/(365/7)
outputName = "Phylogenetic_analyses/All_dispersal_statistics/Nextstrain_2004"
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
if (showingPlots)
	{
		dev.new(width=5, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
		tab1 = read.table("Phylogenetic_analyses/All_dispersal_statistics/Nextstrain_2004_median_weighted_branch_dispersal_velocity.txt", header=T)
		tab2 = read.table("Phylogenetic_analyses/All_dispersal_statistics/Nextstrain_2004_95%HPD_weighted_branch_dispersal_velocity.txt", header=T)
		tab1[,2] = tab1[,2]/366; tab2[,2:3] = tab2[,2:3]/366 # to have the lineage dispersal velocity in km/day
		mcc = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T)
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,3000/366), xlim=c(minYear,maxYear), col=NA)
		slicedTimes = tab1[,1]; branchDispersalVelocityMeanValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
		lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=col1)
		ats = c(minYear, decimal_date(ymd("2020-03-14")), decimal_date(ymd("2020-03-18")), maxYear)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=ats)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
	}

	# 6.9. Analysing the outputs of the overall discrete phylogeographic analysis

tree = readAnnotatedNexus("Phylogenetic_analyses/Nextstrain_2004_2_MCC.tree")
rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingYear-rootHeight
locations = c("Belgium","Africa","Asia","England","Europe","France","Germany","Italy","Luxembourg","NorthAmerica","Oceania","SouthAmerica","Spain")
locationNames = c("Belgium","Africa","Asia","England","Europe","France","Germany","Italy","Luxembourg","North America","Oceania","South America","Spain")
colours = c("gray30",colorRampPalette(brewer.pal(11,"Spectral"))(length(locations)-1)); ancestralLocations = c()
if (showingPlots)
	{
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						location = tree$annotations[[i]]$location
						cols[i] = colours[which(locations==location)]
						if (location == "Belgium") lwds[i] = 0.5
					}
			}
		dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((tree$edge[i,2]%in%tree$edge[,1])&(tree$annotations[[i]]$location!="Belgium"))
					{
						indices = which(tree$edge[,1]==tree$edge[i,2])
						tipLocations = c(tree$annotations[[indices[1]]]$location, tree$annotations[[indices[2]]]$location)
						if (sum(tipLocations=="Belgium") > 0)
							{
								location = tree$annotations[[i]]$location
								col = colours[which(locations==location)]
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col=col)
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
								ancestralLocations = c(ancestralLocations, location)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		dev.new(width=3.5, height=4); par(oma=c(0,0,0,0)); plot.new()
		legend(0.2, 1.09, locationNames, col=colours, text.col="gray30", pch=16, pt.cex=1.2, box.lty=0, cex=0.7, y.intersp=1.3)
	}

