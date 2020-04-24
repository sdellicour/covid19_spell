library(diagram)
library(lubridate)
library(seraphim)
library(treeio)

# 1. Preparing the input files for the discrete phylogeographic analyses 
# 2. Analysing the outputs of the preliminary discrete phylogeographic analysis 
# 3. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)
# 4. Running BEAST and building the maximum clade consensus (MCC) tree
# 5. Extracting spatio-temporal information embedded in MCC and posterior trees
# 6. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)
# 7. Estimating and plotting dispersal statistics associated with lineages
# 8. Analysing the outputs of the overall discrete phylogeographic analysis

writingFiles = FALSE
showingPlots = FALSE

# 1. Preparing the input files for the discrete phylogeographic analyses 

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
colnames(tab1) = c("trait","location","collection_date")
colnames(tab2) = c("trait","location","collection_date")
colnames(tab3) = c("trait","location","collection_date")
if (writingFiles) write.table(tab1, "Phylogenetic_analyses/Nextstrain_200420_1.txt", row.names=F, quote=F, sep="\t")
if (writingFiles) write.table(tab2, "Phylogenetic_analyses/Nextstrain_200420_2.txt", row.names=F, quote=F, sep="\t")
if (writingFiles) write.table(tab3, "Phylogenetic_analyses/Nextstrain_200420_3.txt", row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, "Phylogenetic_analyses/Nextstrain_200420.fasta")

# 2. Analysing the outputs of the preliminary discrete phylogeographic analysis 

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
			" identified from the global phylogenetic analysis of ",belgianTipBranches," SARS-CoV-2 sampled in Belgium (20-04-2020)",sep="")
		# A minimum number of 165 lineage introductions (95% HPD interval = [155-177]) identified from the global phylogenetic analysis of 391 SARS-CoV-2 sampled in Belgium 
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
pop = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(communes)); pop[] = log(pop[])
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
clusterSizes = rep(NA, length(clusters1))
collectionDates = c()
for (i in 1:length(clusters1))
	{
		clusterSizes[i] = length(clusters1[[i]])
		collectionDates = c(collectionDates, clusters2[[i]][,"collectionDate"])
	}
if (showingPlots)
	{
		dev.new(width=3.3, height=8); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,1,1), lwd=0.2, col="gray30")
		hist(clusterSizes, breaks=35, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		hist(collectionDates, breaks=50, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=decimal_date(ymd(c("2020-02-03","2020-02-18","2020-03-03","2020-03-18","2020-04-03"))),
			 labels=c("03-02-2020","18-02-2020","03-03-2020","18-03-2020","03-04-2020"))
	}
if (showingPlots)
	{
		dev.new(width=7, height=6); par(oma=c(0,0,0,0), mar=c(0.5,4,1,0), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"BuPu"))(161)[1:101])
		cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
		plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
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
tryringIBDlikeAnalysis = FALSE
if (tryringIBDlikeAnalysis)
	{
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
	}

# 3. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)

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
if (writingFiles) write.table(tab4, "Phylogenetic_analyses/Nextstrain_200420_4.txt", row.names=F, quote=F, sep="\t")

analyses = c()
template = scan("Phylogenetic_analyses/Nextstrain_RRW_temp.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3))
			{
				analyses = c(analyses, paste0("Clade_",i))
				if (!file.exists(paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".xml")))
					{
						xml = gsub("TEMPLATE", paste0("Clade_",i), template)
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
	}

# 4. Running BEAST and building the maximum clade consensus (MCC) tree

source("MCC_tree_extraction.r")
sink(file=paste0("Phylogenetic_analyses/Phylogeographic_runs/Analyses.sh"))
for (i in 1:length(analyses))
	{
		cat(paste0("java -jar beast_1104.jar -overwrite ",analyses[i],".xml\n"))
	}
sink(NULL)
runningNewAnalyses = FALSE
wd = getwd(); setwd(paste0(wd,"/Phylogenetic_analyses/Phylogeographic_runs/"))
if (runningNewAnalyses)
	{
		system("bash Analyses.sh", ignore.stdout=T, ignore.stderr=F)
		for (i in 1:length(analyses))
			{
				system(paste0("BEAST_1104/bin/treeannotator -burninTrees 101 -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

# 5. Extracting spatio-temporal information embedded in MCC and posterior trees

wd = getwd(); setwd(paste0(wd,"/Phylogenetic_analyses/Phylogeographic_runs/"))
for (i in 1:length(analyses))
	{
		if (!file.exists(paste0(analyses[i],".csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree"))
				mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
			}
	}
nberOfTreesToSample = 1000; burnIn = 101; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				allTrees = scan(file=paste0(analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
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

# 6. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)

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
		cexNode = 0.8; LWD = 1.0
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		firstTimePeriod = TRUE; secondTimePeriod = FALSE
		firstTimePeriod = FALSE; secondTimePeriod = TRUE
		firstTimePeriod = FALSE; secondTimePeriod = FALSE
		if ((firstTimePeriod == TRUE)|(secondTimePeriod == TRUE)) { cexNode = 1.1; LWD = 2.0 }
		dev.new(width=7.3, height=6); par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
		cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
		if ((firstTimePeriod == TRUE)|(secondTimePeriod == TRUE))
			{
				plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(provinces, border="white", col=NA, add=T, lwd=LWD)
			}	else		{
				plot(provinces, border=NA, col="gray95", lwd=LWD)
			}
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
		plot(provinces, border="white", col=NA, add=T, lwd=LWD)
		plot(belgium, border="gray30", col=NA, add=T, lwd=0.4); croppingPolygons = TRUE
		selectedBranches = 1:dim(mcc)[1]
		if (firstTimePeriod) selectedBranches = which(mcc[,"endYear"]<decimal_date(ymd("2020-03-18")))
		if (secondTimePeriod) selectedBranches = which(mcc[,"startYear"]>decimal_date(ymd("2020-03-18")))
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
							{
								if (sum(!is.na(clusters2[[i]][,"longitude"])) == 2)
									{
										indices = which(!is.na(clusters2[[i]][,"longitude"]))
										if (firstTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]<decimal_date(ymd("2020-03-18"))))
											}
										if (secondTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]>decimal_date(ymd("2020-03-18"))))
											}
										if (length(indices) == 2)
											{
												curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
															cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
															arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
															curve=0.1, dr=NA, endhead=F)
											}
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
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
											}
									}
							}
					}
			}
		for (i in rev(selectedBranches))
			{
				if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode,)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
			}
		selectedDates = decimal_date(ymd(c("2020-03-03","2020-03-18","2020-04-03")))
		selectedLabels = c("03-03-2020","18-03-2020","03-04-2020")
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.1,0.5,0.100,0.112),
			 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.00,0),
		     at=selectedDates, labels=selectedLabels))
	}
if (showingPlots)
	{
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0
		for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
			{
				if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
					{
						maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
					}
			}
		pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pop_liege = raster::mask(crop(pop,pol),pol)
		vS1 = raster::extract(pop_liege, mcc[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, mcc[,c("endLon","endLat")])
		sub = mcc[which((!is.na(vS1))&(!is.na(vS2))),]		
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[21:121])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((sub[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((sub[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		dev.new(width=12, height=3); par(mfrow=c(1,3), oma=c(0,2,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		for (i in 1:2)
			{
				cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
				cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
				plot(pop_liege, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(pol, border="gray30", col=NA, add=T, lwd=0.4)
				if (i == 1) selectedBranches = which(sub[,"endYear"]<decimal_date(dmy("18-03-2020")))
				if (i == 2) selectedBranches = which(sub[,"startYear"]>=decimal_date(dmy("18-03-2020")))
				for (j in selectedBranches)
					{
						curvedarrow(cbind(sub[j,"startLon"],sub[j,"startLat"]), cbind(sub[j,"endLon"],sub[j,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (j in rev(selectedBranches))
					{
						if (!sub[j,"node1"]%in%sub[selectedBranches,"node2"])
							{
								points(sub[j,"startLon"], sub[j,"startLat"], pch=16, col=startYears_colours[j], cex=cexNode)
								points(sub[j,"startLon"], sub[j,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
							}
						points(sub[j,"endLon"], sub[j,"endLat"], pch=16, col=endYears_colours[j], cex=cexNode,)
						points(sub[j,"endLon"], sub[j,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
			}
		plot(pop_liege, col=cols, axes=F, ann=F, box=F, legend=F)
		plot(pol, border="gray30", col=NA, add=T, lwd=0.4)
		legendRast = raster(as.matrix(c(min(rast[],na.rm=T),max(rast[],na.rm=T))))
		mtext("Human population (log-transformed)", col="gray30", cex=0.7, line=-19.5, at=731000)
		plot(legendRast, legend.only=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.1,0.58,0.11,0.13),
			 alpha=1, horizontal=T, legend.args=list(text="", cex=1.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.9, lwd=0,
			 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.3,0)))
	}

# 7. Estimating and plotting dispersal statistics associated with lineages
	
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext2")
timSlices = 50; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 5; slidingWindow = 1/(365/7)
nberOfExtractionFiles = 1000; outputName = "Phylogenetic_analyses/All_dispersal_statistics/50_time_slices_7days_sliding_window/Nextstrain_2004"
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
if (showingPlots)
	{
		dev.new(width=5, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
		directory = "Phylogenetic_analyses/All_dispersal_statistics/50_time_slices_7days_sliding_window/"
		tab1 = read.table(paste0(directory,"Nextstrain_2004_median_weighted_branch_dispersal_velocity.txt"), header=T)
		tab2 = read.table(paste0(directory,"Nextstrain_2004_95%HPD_weighted_branch_dispersal_velocity.txt"), header=T)
		tab1[,2] = tab1[,2]/366; tab2[,2:3] = tab2[,2:3]/366 # to have the lineage dispersal velocity in km/day
		mcc = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T)
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,3000/366), xlim=c(minYear,maxYear), col=NA)
		slicedTimes = tab1[,1]; branchDispersalVelocityMeanValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
		lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=col1)
		ats = c(minYear, decimal_date(ymd("2020-03-01")), decimal_date(ymd("2020-03-14")), decimal_date(ymd("2020-03-28")), maxYear)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("","01-03","14-03","28-03",""))
	}

mcc = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T)
geoDistances = rep(NA, dim(mcc)[1]); branchDurations = rep(NA, dim(mcc)[1])
colours = c("#FAA521","#4676BB"); cols = rep(NA, dim(mcc)[1])
for (i in 1:dim(mcc)[1])
	{
		geoDistances[i] = sqrt(((mcc[i,"startLon"]-mcc[i,"endLon"])^2)+((mcc[i,"startLat"]-mcc[i,"endLat"])^2))
		dt = mcc[i,"endYear"]-mcc[i,"startYear"]; branchDurations[i] = dt*366
		if (mcc[i,"endYear"] < decimal_date(ymd("2020-03-18"))) cols[i] = colours[1]
		if (mcc[i,"startYear"] > decimal_date(ymd("2020-03-18"))) cols[i] = colours[2]
	}
if (showingPlots)
	{
		dev.new(width=5, height=4); par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(3,3.5,1.5,2))
		cols1 = cols; cols2 = cols1; cols2[which(!is.na(cols2))] = paste0(cols2[which(!is.na(cols2))],"50")
		plot(log(geoDistances), log(branchDurations), col=cols2, pch=16, cex=0.8, xlim=c(3.1,11.7), ylim=c(-5.9,3.2), axes=F, ann=F)
		points(log(geoDistances), log(branchDurations), col=cols1, pch=1, cex=0.8)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, at=seq(2,12,1), col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, at=seq(-7,4,1), col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="phylogenetic branch durations (days, log-transformed)", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		title(xlab="geographic distance (km, log-transformed)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
	}

# 8. Analysing the outputs of the overall discrete phylogeographic analysis

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

