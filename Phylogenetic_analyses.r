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
# 7. Comparing intra- and inter-province/municipality lineage migration events
# 8. Estimating and plotting dispersal statistics associated with lineages

analysis = "Nextstrain_200420"; removeSuspiciousHomoplasies = FALSE # 1st analysis
analysis = "Nextstrain_200420"; removeSuspiciousHomoplasies = TRUE  # 2nd analysis
analysis = "Nextstrain_200420"; removeSuspiciousHomoplasies = FALSE # 3rd analysis
analysis = "TreeTime_260420_1"; removeSuspiciousHomoplasies = FALSE # 4th analysis
analysis = "TreeTime_260420_2"; removeSuspiciousHomoplasies = FALSE # 5th analysis
analysis = "TreeTime_100620"; removeSuspiciousHomoplasies = FALSE # 6th analysis

	# 1st analysis: analysis based on the Nextstrain tree of the 20-04-20
	# 2nd analysis: analysis based on the Nextstrain tree of the 20-04-20 but for which several branches are dropped, i.e. branches corresponding 
	#		to 21 sequences with suspicious homoplasies: C3130T, A4050C, T8022G, T13402G, C27046T, T28785G (see also the csv file "Suspicious_homopl")
	# 3rd analysis: analysis based on the Nextstrain tree of the 20-04-20, but which was based on an alignment without these 21 sequences
	# 4th analysis: analysis based on a new TreeTime tree, which was based on an GISAID alignment of the 26-04-20 but without these 21 sequences
	# 5th analysis: analysis based on a new TreeTime tree, which was based on an GISAID alignment of the 26-04-20 but with the following sites masked
	# 		because carrying suspicious homoplasies listed by De Maio and colleagues: 11083, 21575, 11074, 6990, 29353, 16887, 10323, 3145, 4050, 13408,
	# 		13402, 14408, 3130, 8022, 15324 - and a second-tier list of sites that were also masked: 6255, 21137, 14724, 26681, 28077, 241, 
	# 		335, 2094, 3037, 8782, 9223, 10741, 11704, 14786, 14805, 17247, 19684, 24034, 24378, 25563, 26461, 27384, 28826, 28854, 29700
	# 6th analysis: analysis based on a new TreeTime tree, which was based on all the Belgian sequences available on GISAID on the 10-06-20 plus non-Belgian
	#		sequences available on GISAID and included in Nextstrain (which used a subsampling per country to avoid including too many sequences)
	
data1 = read.csv("Sequences_metadata/SARS-CoV-2_KULeuven_100620.csv", sep=";")
data2 = read.csv("Sequences_metadata/SARS-CoV-2_ULiegeSeq_020620.csv", sep=";")
writingFiles = FALSE; showingPlots = FALSE

# 1. Preparing the input files for the discrete phylogeographic analyses 

tree = read.tree(paste0("Phylogenetic_analyses/",analysis,".tre"))
if (grepl("TreeTime",analysis) == TRUE)
	{
		seqIDs = tree$tip.label; countries = rep(NA, length(seqIDs)); collectionDates = rep(NA, length(seqIDs))
		for (i in 1:length(seqIDs))
			{
				if (grepl("hCoV-19",seqIDs[i]))
					{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[2]
					}	else	{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[1]
					}
				if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 3)
					{
						collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))]
					}
				if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 4)
					{
						collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))-1]
					}
			}
		tab = cbind(seqIDs,countries,collectionDates); colnames(tab) = c("Strain","Country","Collection Data")
		write.csv(tab, paste0("Phylogenetic_analyses/",analysis,".csv"), row.names=F, quote=F)
		data = read.csv(paste0("Phylogenetic_analyses/",analysis,".csv"))
	}	else	{
		data = read.csv(paste0("Phylogenetic_analyses/",analysis,".csv"), sep=";")
		# N.B.: for Nextstrain, the ".tsv" file has first to be opened in Excel and then exported as comma ".csv"
	}
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
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("ULG-",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
	}
if (removeSuspiciousHomoplasies)
	{
		suspicious_homoplasies = read.csv("Phylogenetic_analyses/Suspicious_homopl.csv")[,"sequence"]
		suspicious_homoplasies = as.character(unique(suspicious_homoplasies[!is.na(suspicious_homoplasies)]))
		labs = unique(data[which(data[,"Country"]=="Belgium"),"Submitting.Lab"])
		sedIDs = gsub("Belgium\\/","",gsub("\\/2020","",data[which(data[,"Country"]=="Belgium"),"Strain"]))
		sedIDs_KUL = sedIDs[(!grepl("ULG-",sedIDs))&(!grepl("UGent-",sedIDs))]
		temp = data1[which(data1[,"GisAID"]=="OK"),]
		temp[,"sequence.name"] = gsub("SARS2-CoV\\/Belgium\\/Human\\/","",gsub("\\/2020","",temp[,"sequence.name"]))
		temp[,"sequence.name"] = gsub("SARS2-Cov\\/Belgium\\/Human\\/","",gsub("\\/2020","",temp[,"sequence.name"]))
		indices = which(!temp[,"sequence.name"]%in%sedIDs_KUL)
		discardedKULsequences = temp[indices,"sequence.name"]
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
						if ((!tree$edge[i,2]%in%tree$edge[,1]) & (tree$tip.label[tree$edge[i,2]]%in%suspicious_homoplasies))
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="red")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
				add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
			}
		tree = drop.tip(tree, suspicious_homoplasies)
		if (writingFiles) write.nexus(tree, file="Newick_tree_for_XML_file.tree")
	}
txt = c(); tab = c()
for (i in 1:length(tree$tip.label))
	{
		index = which(data[,1]==tree$tip.label[i])
		date = as.character(data[index,"Collection.Data"])
		if (date != "")
			{
				txt = c(txt, paste0(">",tree$tip.label[i]),"NNNN")
				if (grepl("Nextstrain",analysis))
					{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[1]
					}
				if (grepl("TreeTime",analysis))
					{
						if (grepl("hCoV-19",tree$tip.label[i]))
							{
								location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
							}	else	{
								location = unlist(strsplit(tree$tip.label[i],"\\/"))[1]
							}
					}
				if (location != "Belgium") location = "other"
				tab = rbind(tab, cbind(tree$tip.label[i],location,date))
			}
	}
colnames(tab) = c("trait","location","collection_date")
if (writingFiles) write.table(tab, paste0("Phylogenetic_analyses/",analysis,".txt"), row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, paste0("Phylogenetic_analyses/",analysis,".fasta"))

# 2. Analysing the outputs of the preliminary discrete phylogeographic analysis 

burnIn = 101; computingHPDInterval = FALSE # N.B.: long analysis
if (computingHPDInterval)
	{
		trees = readAnnotatedNexus("Phylogenetic_analyses/TreeTime_100620.trees")
		belgianBranches_list = rep(NA,length(trees))
		belgianIntroductions_list = rep(NA,length(trees))
		belgianTipBranches_list = rep(NA,length(trees))
		for (i in burnIn:length(trees))
			{
				belgianBranches = 0; belgianIntroductions = 0; belgianTipBranches = 0
				for (j in 1:dim(trees[[i]]$edge)[1])
					{
						if (trees[[i]]$annotations[[j]]$location == "Belgium")
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
		quantiles = quantile(belgianIntroductions_list[!is.na(belgianIntroductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(belgianIntroductions_list[!is.na(belgianIntroductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",belgianTipBranches," SARS-CoV-2 sampled in Belgium (20-04-2020)",sep="")
		# Results for the 1° analysis based on the Nextstrain tree of the 20-04-20 (without removing tip branches associated with suspicious homoplasie):
			# a minimum number of 166 lineage introductions (95% HPD interval = [161-171]) identified from the global phylogenetic analysis of 391 SARS-CoV-2 sampled in Belgium
		# Results for the 2° analysis based on the Nextstrain tree of the 20-04-20 (when dropping tip branches associated with suspicious homoplasie):
			# a minimum number of 157 lineage introductions (95% HPD interval = [151-162]) identified from the global phylogenetic analysis of 370 SARS-CoV-2 sampled in Belgium
		# Results for the 3° analysis based on the Nextstrain tree of the 20-04-20 (based on a lignment without sequences associated with suspicious homoplasie):
			# a minimum number of 144 lineage introductions (95% HPD interval = [138-150]) identified from the global phylogenetic analysis of 370 SARS-CoV-2 sampled in Belgium
		# Results for the 6° analysis based on the Nextstrain tree of the 10-06-20 (based on a lignment without sequences associated with suspicious homoplasie):
			# a minimum number of XXX lineage introductions (95% HPD interval = [XXX-XXX]) identified from the global phylogenetic analysis of XXX SARS-CoV-2 sampled in Belgium
	}
tree = readAnnotatedNexus("Phylogenetic_analyses/TreeTime_100620.tree")
if (showingPlots)
	{
		samplingDates = decimal_date(ymd(gsub("\\/","-",tab[,"collection_date"]))); mostRecentSamplingYear = max(samplingDates)
		selectedDates = decimal_date(ymd(c("2019-11-01","2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01")))
		rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingYear-rootHeight
		selectedLabels = c("01-11-2019","01-12-2019","01-01-2020","01-02-2020","01-03-2020","01-04-2020","01-05-2020","01-06-2020")
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
				clusters1[[i]] = gsub("'","",subtree$tip.label)
			}	else		{
				clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[belgianIntroductions[i],2]])
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
				date = ymd(gsub("\\/","-",data[index,"Collection.Data"]))
				samplingData[i,"collectionDate"] = decimal_date(date)
				ID = unlist(strsplit(samplingData[i,"sequenceID"],"\\/"))[3]
				index1 = which(grepl(gsub("Rega","rega",ID),data1[,"sequence.name"]))
				if (length(index1) == 1)
					{
						samplingData[i,"postCode"] = data1[index1,"ZIP"]
						indices = which(communes@data[,"nouveau_PO"]==data1[index1,"ZIP"])
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
		centroids[[i]] = cbind(mean(tab[!is.na(tab[,"longitude"]),"longitude"]), mean(tab[!is.na(tab[,"latitude"]),"latitude"]))
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
		collectionDates_filetered = collectionDates#[which(collectionDates>decimal_date(ymd("2020-03-01")))]
		dev.new(width=3.3, height=8); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,1,1), lwd=0.2, col="gray30")
		hist(clusterSizes, breaks=40, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		hist(collectionDates_filetered, breaks=65, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=decimal_date(ymd(c("2020-02-01","2020-03-01","2020-04-01","2020-05-01"))),
			 labels=c("01-02-2020","01-03-2020","01-04-2020","01-05-2020"))
	}
if (showingPlots)
	{
		dev.new(width=7.3, height=6); par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
		plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
		plot(provinces, border="white", col=NA, add=T, lwd=1.0)
		plot(belgium, border="gray30", col=NA, add=T, lwd=0.4)
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
												segments(centroids[[i]][,1],centroids[[i]][,2],clusters2[[i]][j,"longitude"],clusters2[[i]][j,"latitude"], lwd=0.5, col="gray30")	
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
								points(clusters2[[i]][,"longitude"], clusters2[[i]][,"latitude"], pch=16, cex=0.8, col="chartreuse3")
								points(clusters2[[i]][,"longitude"], clusters2[[i]][,"latitude"], pch=1, cex=0.8, col="gray30", lwd=0.2)
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
		legendRast = raster(as.matrix(c(min(pop[],na.rm=T),max(pop[],na.rm=T))))
		mtext("Human population (log-transformed)", col="gray30", cex=0.7, line=-23, at=601000)
		plot(legendRast, legend.only=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.141,0.409,0.18,0.19),
			 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.55, lwd=0,
			 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,-0.05,0)))
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
if (writingFiles) write.table(tab4, paste0("Phylogenetic_analyses/",analysis,"_4.txt"), row.names=F, quote=F, sep="\t")

analyses = c()
template = scan("Phylogenetic_analyses/RRW_XMLtemplate.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3))
			{
				analyses = c(analyses, paste0("Clade_",i)); cluster = clusters2[[i]][which(!is.na(clusters2[[i]][,"longitude"])),]
				xml = gsub("TEMPLATE", paste0("Clade_",i), template)
				tre = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				tips_to_drop = tre$tip.label[which(!tre$tip.label%in%row.names(cluster))]
				if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
				write.tree(tre, paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".tre"))
				tre = scan(paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".tre"), what="", sep="\n", quiet=T)
				sink(file=paste0("Phylogenetic_analyses/Phylogeographic_runs/Clade_",i,".xml"))
				for (j in 1:length(xml))
					{
						cat(xml[j]); cat("\n")
						if (xml[j]=="\t<taxa id=\"taxa\">")
							{
								for (k in 1:dim(cluster)[1])
									{
										cat(paste0("\t\t<taxon id=\"",row.names(cluster)[k],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",cluster[k,"collectionDate"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"latitude\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"latitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"longitude\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"coordinates\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"latitude"]," ",cluster[k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
							}
						if (xml[j]=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")
							{
								for (k in 1:dim(cluster)[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",row.names(cluster)[k],"\"/>","\n"))
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
		burnIns = rep(1001, length(analyses)); burnIns[which(analyses%in%paste0("Clade_",c(5,96,123)))] = 5001
		for (i in 1:length(analyses))
			{
system(paste0("BEAST_1104/bin/treeannotator -burninTrees ",burnIns[i]," -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

# 5. Extracting spatio-temporal information embedded in MCC and posterior trees

prov_WGS = raster::getData("GADM", country="BEL", level=2)
polIndex1 = which(prov_WGS@data[,"NAME_2"]=="Liège")
maxArea = 0; polIndex2 = 0
for (i in 1:length(prov_WGS@polygons[[polIndex1]]@Polygons))
	{
		if (maxArea < prov_WGS@polygons[[polIndex1]]@Polygons[[i]]@area)
			{
				maxArea = prov_WGS@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
			}
	}
pol = prov_WGS@polygons[[polIndex1]]@Polygons[[polIndex2]]
p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
pol = sps; proj4string(pol) = crs(raster("WorldPop_pop_raster.tif"))
pop_liege_WGS = raster::mask(crop(raster("WorldPop_pop_raster.tif"),pol),pol)
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
nberOfTreesToSample = 1000; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2]); burnIn = burnIns[i]
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
dir.create(file.path("Bf_180320_ext2"), showWarnings=F)
dir.create(file.path("Af_180320_ext2"), showWarnings=F)
dir.create(file.path("Bf_180320_Liege"), showWarnings=F)
dir.create(file.path("Af_180320_Liege"), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles)
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
		bf1 = all[which(all[,"endYear"]<decimal_date(dmy("18-03-2020"))),]
		af1 = all[which(all[,"endYear"]>=decimal_date(dmy("18-03-2020"))),]
		write.csv(bf1, paste0("Bf_180320_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af1, paste0("Af_180320_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		vS1 = raster::extract(pop_liege_WGS, bf1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, bf1[,c("endLon","endLat")])
		bf2 = bf1[which((!is.na(vS1))&(!is.na(vS2))),]
		vS1 = raster::extract(pop_liege_WGS, af1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, af1[,c("endLon","endLat")])
		af2 = af1[which((!is.na(vS1))&(!is.na(vS2))),]
		write.csv(bf2, paste0("Bf_180320_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af2, paste0("Af_180320_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
setwd(wd)

# 6. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)

localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext1"); nberOfExtractionFiles = 1000
percentage = 80; prob = percentage/100; precision = 1/(365/7); croppingPolygons = TRUE
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
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
			}
		selectedDates = decimal_date(ymd(c("2020-02-03","2020-03-03","2020-03-14","2020-03-18","2020-04-03","2020-05-03")))
		selectedLabels = c("03-02-2020","03-03-2020","","18-03-2020","03-04-2020","03-05-2020")
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.1,0.5,0.100,0.112),
			 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.00,0),
		     at=selectedDates, labels=selectedLabels))
	}
if (showingPlots)
	{
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0; cexNode = 0.85
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
		dev.new(width=11, height=4); par(mfrow=c(1,2), oma=c(0,2,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
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
						points(sub[j,"endLon"], sub[j,"endLat"], pch=16, col=endYears_colours[j], cex=cexNode)
						points(sub[j,"endLon"], sub[j,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
			}
	}

# 7. Comparing intra- and inter-province/municipality lineage migration events

communes = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
provinces_sf = sf::st_as_sf(provinces); communes_sf = sf::st_as_sf(communes)
nberOfExtractionFiles = 1000; localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext1")
for (i in 0:nberOfExtractionFiles)
	{
		if (i == 0) tab = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T)
		if (i >= 1) tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
		startEndProvinces = matrix(nrow=dim(tab)[1], ncol=2); colnames(startEndProvinces) = c("startProvince","endProvince")
		pts = sf::st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(provinces))
		indices1 = sf::st_intersects(pts, provinces_sf); indices2 = rep(NA, length(indices1))
		for (j in 1:length(indices1))
			{
				if (length(indices1[[j]]) > 0) indices2[j] = indices1[[j]]
			}
		startEndProvinces[,"startProvince"] = gsub(" ","",provinces@data[indices2,"NAME_2"])
		pts = sf::st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(provinces))
		indices1 = sf::st_intersects(pts, provinces_sf); indices2 = rep(NA, length(indices1))
		for (j in 1:length(indices1))
			{
				if (length(indices1[[j]]) > 0) indices2[j] = indices1[[j]]
			}
		startEndProvinces[,"endProvince"] = gsub(" ","",provinces@data[indices2,"NAME_2"])
		tab = cbind(tab, startEndProvinces)
		startEndCommunes = matrix(nrow=dim(tab)[1], ncol=2); colnames(startEndCommunes) = c("startCommune","endCommune")
		pts = sf::st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(communes))
		indices1 = sf::st_intersects(pts, communes_sf); indices2 = rep(NA, length(indices1))
		for (j in 1:length(indices1))
			{
				if (length(indices1[[j]]) > 0) indices2[j] = indices1[[j]]
			}
		startEndCommunes[,"startCommune"] = gsub(" ","",communes@data[indices2,"nouveau_PO"])
		pts = sf::st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(communes))
		indices1 = sf::st_intersects(pts, communes_sf); indices2 = rep(NA, length(indices1))
		for (j in 1:length(indices1))
			{
				if (length(indices1[[j]]) > 0) indices2[j] = indices1[[j]]
			}
		startEndCommunes[,"endCommune"] = gsub(" ","",communes@data[indices2,"nouveau_PO"])
		tab = cbind(tab, startEndCommunes)
		if (i == 0) write.csv(tab, "Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", row.names=F, quote=F)
		if (i >= 1) write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
timeSlices = 100; slidingWindow = 14/366
minStartYear = 2021; maxEndYear = 0
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
		if (minStartYear > min(tab[,"startYear"])) minStartYear = min(tab[,"startYear"])
		if (maxEndYear < max(tab[,"endYear"])) maxEndYear = max(tab[,"endYear"])
	}
timeInterval = (maxEndYear-minStartYear)/timeSlices
startEndTimes = matrix(nrow=timeSlices, ncol=3)
for (i in 1:timeSlices)
	{
		time = minStartYear+((i-1)*timeInterval)+(timeInterval/2)
		startTime = time - (slidingWindow/2)
		endTime = time + (slidingWindow/2)
		startEndTimes[i,1:3] = cbind(time, startTime, endTime)
	}
directory = "Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_all_branches_ext1/"
if (!file.exist(paste0(directory,"All_within_province_transition_events.csv")))
	{
		withinProvinceTransitions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)
		amongProvincesTransitions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)		
		amongProvincesProportions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)
		withinCommuneTransitions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)
		amongCommunesTransitions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)		
		amongCommunesProportions = matrix(0, nrow=timeSlices, ncol=nberOfExtractionFiles)
		for (t in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",t,".csv"), head=T)
				tab = tab[with(tab,order(endYear, startYear)),]; tab = tab[order(tab[,"endYear"]),]
				for (i in 1:timeSlices)
					{
						time = startEndTimes[i,1]
						startTime = startEndTimes[i,2]
						endTime = startEndTimes[i,3]
						for (j in 1:dim(tab)[1])
							{
								branchInTimeInterval = FALSE
								if ((tab[j,"startYear"]<startTime)&(tab[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
								if ((tab[j,"startYear"]>startTime)&(tab[j,"startYear"]<endTime)) branchInTimeInterval = TRUE
								if ((tab[j,"endYear"]>startTime)&(tab[j,"endYear"]<endTime)) branchInTimeInterval = TRUE
								if (branchInTimeInterval == TRUE)
									{
										if ((!is.na(tab[j,"startProvince"]))&(!is.na(tab[j,"endProvince"])))
											{
								if (as.character(tab[j,"startProvince"]) == as.character(tab[j,"endProvince"])) withinProvinceTransitions[i,t] = withinProvinceTransitions[i,t]+1
								if (as.character(tab[j,"startProvince"]) != as.character(tab[j,"endProvince"])) amongProvincesTransitions[i,t] = amongProvincesTransitions[i,t]+1
											}
										if ((!is.na(tab[j,"startCommune"]))&(!is.na(tab[j,"endCommune"])))
											{
								if (as.character(tab[j,"startCommune"]) == as.character(tab[j,"endCommune"])) withinCommuneTransitions[i,t] = withinCommuneTransitions[i,t]+1
								if (as.character(tab[j,"startCommune"]) != as.character(tab[j,"endCommune"])) amongCommunesTransitions[i,t] = amongCommunesTransitions[i,t]+1
											}
									}
							}
						amongProvincesProportions[i,t] = amongProvincesTransitions[i,t]/(withinProvinceTransitions[i,t]+amongProvincesTransitions[i,t])
						amongCommunesProportions[i,t] = amongCommunesTransitions[i,t]/(withinCommuneTransitions[i,t]+amongCommunesTransitions[i,t])
					}
			}
		row.names(withinProvinceTransitions) = startEndTimes[,1]; row.names(amongProvincesTransitions) = startEndTimes[,1]
		write.csv(withinProvinceTransitions, paste0(directory,analysis,"_all_within_province_transition_events.csv"), quote=F)
		write.csv(amongProvincesTransitions, paste0(directory,analysis,"_all_among_provinces_transition_events.csv"), quote=F)
		row.names(withinCommuneTransitions) = startEndTimes[,1]; row.names(amongCommunesTransitions) = startEndTimes[,1]
		write.csv(withinCommuneTransitions, paste0(directory,analysis,"_all_within_commune_transition_events.csv"), quote=F)
		write.csv(amongCommunesTransitions, paste0(directory,analysis,"_all_among_communes_transition_events.csv"), quote=F)
		row.names(amongProvincesProportions) = startEndTimes[,1]; row.names(amongCommunesProportions) = startEndTimes[,1]
		write.csv(amongProvincesProportions, paste0(directory,analysis,"_all_among_provinces_proportion_events.csv"), quote=F)
		write.csv(amongCommunesProportions, paste0(directory,analysis,"_all_among_communes_proportion_events.csv"), quote=F)
	}	else	{
		withinProvinceTransitions = as.matrix(read.csv(paste0(directory,analysis,"_all_within_province_transition_events.csv"), header=T))
		amongProvincesTransitions = as.matrix(read.csv(paste0(directory,analysis,"_all_among_provinces_transition_events.csv"), header=T))
		amongProvincesProportions = as.matrix(read.csv(paste0(directory,analysis,"_all_among_provinces_proportion_events.csv"), header=T))
		withinCommuneTransitions = as.matrix(read.csv(paste0(directory,analysis,"_all_within_commune_transition_events.csv"), header=T))
		amongCommunesTransitions = as.matrix(read.csv(paste0(directory,analysis,"_all_among_communes_transition_events.csv"), header=T))
		amongCommunesProportions = as.matrix(read.csv(paste0(directory,analysis,"_all_among_communes_proportion_events.csv"), header=T))
	}
withinProvinceTransitions_median_95HPD = matrix(nrow=dim(withinProvinceTransitions)[1], ncol=3)
amongProvincesTransitions_median_95HPD = matrix(nrow=dim(amongProvincesTransitions)[1], ncol=3)
amongProvincesProportions_median_95HPD = matrix(nrow=dim(amongProvincesProportions)[1], ncol=3)
withinCommuneTransitions_median_95HPD = matrix(nrow=dim(withinCommuneTransitions)[1], ncol=3)
amongCommunesTransitions_median_95HPD = matrix(nrow=dim(amongCommunesTransitions)[1], ncol=3)
amongCommunesProportions_median_95HPD = matrix(nrow=dim(amongCommunesProportions)[1], ncol=3)
for (i in 1:dim(withinProvinceTransitions)[1])
	{
		withinProvinceTransitions_median_95HPD[i,1] = median(withinProvinceTransitions[i,])
		withinProvinceTransitions_median_95HPD[i,2] = quantile(withinProvinceTransitions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		withinProvinceTransitions_median_95HPD[i,3] = quantile(withinProvinceTransitions[i,], probs=c(0.025,0.975), na.rm=T)[2]
		amongProvincesTransitions_median_95HPD[i,1] = median(amongProvincesTransitions[i,])
		amongProvincesTransitions_median_95HPD[i,2] = quantile(amongProvincesTransitions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		amongProvincesTransitions_median_95HPD[i,3] = quantile(amongProvincesTransitions[i,], probs=c(0.025,0.975), na.rm=T)[2]
		amongProvincesProportions_median_95HPD[i,1] = median(amongProvincesProportions[i,])
		amongProvincesProportions_median_95HPD[i,2] = quantile(amongProvincesProportions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		amongProvincesProportions_median_95HPD[i,3] = quantile(amongProvincesProportions[i,], probs=c(0.025,0.975), na.rm=T)[2]
		withinCommuneTransitions_median_95HPD[i,1] = median(withinCommuneTransitions[i,])
		withinCommuneTransitions_median_95HPD[i,2] = quantile(withinCommuneTransitions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		withinCommuneTransitions_median_95HPD[i,3] = quantile(withinCommuneTransitions[i,], probs=c(0.025,0.975), na.rm=T)[2]
		amongCommunesTransitions_median_95HPD[i,1] = median(amongCommunesTransitions[i,])
		amongCommunesTransitions_median_95HPD[i,2] = quantile(amongCommunesTransitions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		amongCommunesTransitions_median_95HPD[i,3] = quantile(amongCommunesTransitions[i,], probs=c(0.025,0.975), na.rm=T)[2]
		amongCommunesProportions_median_95HPD[i,1] = median(amongCommunesProportions[i,])
		amongCommunesProportions_median_95HPD[i,2] = quantile(amongCommunesProportions[i,], probs=c(0.025,0.975), na.rm=T)[1]
		amongCommunesProportions_median_95HPD[i,3] = quantile(amongCommunesProportions[i,], probs=c(0.025,0.975), na.rm=T)[2]
	}
if (writingFiles)
	{
		tab = cbind(startEndTimes[,1], withinProvinceTransitions_median_95HPD, amongProvincesTransitions_median_95HPD)
		colnames(tab) = c("time", "within_province_transitions_median", "within_province_transitions_lower95HPD", "within_province_transitions_higher95HPD",
						  "among_provinces_transitions_median", "among_provinces_transitions_lower95HPD", "among_provinces_transitions_higher95HPD")
		write.csv(tab, paste0(directory,analysis,"_within_vs_among_provinces_transitions.csv"), row.names=F, quote=F)
		tab = cbind(startEndTimes[,1], withinCommuneTransitions_median_95HPD, amongCommunesTransitions_median_95HPD)
		colnames(tab) = c("time", "within_commune_transitions_median", "within_commune_transitions_lower95HPD", "within_commune_transitions_higher95HPD",
						  "among_communes_transitions_median", "among_communes_transitions_lower95HPD", "among_communes_transitions_higher95HPD")
		write.csv(tab, paste0(directory,analysis,"_within_vs_among_communes_transitions.csv"), row.names=F, quote=F)
		tab = cbind(startEndTimes[,1], amongProvincesProportions_median_95HPD, amongCommunesProportions_median_95HPD)
		colnames(tab) = c("time", "among_provinces_proportion_median", "among_provinces_proportion_lower95HPD", "among_provinces_proportion_higher95HPD",
						  "among_communes_proportion_median", "among_communes_proportion_lower95HPD", "among_communes_proportion_higher95HPD")
		write.csv(tab, paste0(directory,analysis,"_among_provinces_communes_proportions.csv"), row.names=F, quote=F)
	}

# 8. Estimating and plotting dispersal statistics associated with lineages

nberOfExtractionFiles = 1000; timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 5; slidingWindow = 1/(365/14)
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/All_clades_ext2")
dir.create(file.path("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_all_branches_ext2/"), showWarnings=F)
outputName = paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_all_branches_ext2/",analysis)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/Bf_180320_ext2")
outputName = paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_ext2/",analysis)
dir.create(file.path("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_ext2/"), showWarnings=F)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/Af_180320_ext2")
outputName = paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_ext2/",analysis); dir.create(file.path(outputName), showWarnings=F)
dir.create(file.path("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_ext2/"), showWarnings=F)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/Bf_180320_Liege")
outputName = paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_Liege/",analysis); dir.create(file.path(outputName), showWarnings=F)
dir.create(file.path("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_Liege/"), showWarnings=F)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
localTreesDirectory = paste0("Phylogenetic_analyses/Phylogeographic_runs/Af_180320_Liege")
outputName = paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_Liege/",analysis); dir.create(file.path(outputName), showWarnings=F)
dir.create(file.path("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_Liege/"), showWarnings=F)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 

tab = read.table(paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_ext2/",analysis,"_estimated_dispersal_statistics.txt"), header=T)
vS = round(tab[,"weighted_branch_dispersal_velocity"]/366,1); cat(median(vS)," km/day (95% HPD interval = [",quantile(vS,0.025),"-",quantile(vS,0.975),"])",sep="")
tab = read.table(paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_ext2/",analysis,"_estimated_dispersal_statistics.txt"), header=T)
vS = round(tab[,"weighted_branch_dispersal_velocity"]/366,1); cat(median(vS)," km/day (95% HPD interval = [",quantile(vS,0.025),"-",quantile(vS,0.975),"])",sep="")
tab = read.table(paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Bf_180320_Liege/",analysis,"_estimated_dispersal_statistics.txt"), header=T)
vS = round(tab[,"weighted_branch_dispersal_velocity"]/366,1); cat(median(vS)," km/day (95% HPD interval = [",quantile(vS,0.025),"-",quantile(vS,0.975),"])",sep="")
tab = read.table(paste0("Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_Af_180320_Liege/",analysis,"_estimated_dispersal_statistics.txt"), header=T)
vS = round(tab[,"weighted_branch_dispersal_velocity"]/366,1); cat(median(vS)," km/day (95% HPD interval = [",quantile(vS,0.025),"-",quantile(vS,0.975),"])",sep="")

if (showingPlots)
	{
		directory = "Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_all_branches_ext2/"
		mcc = read.csv("Phylogenetic_analyses/Phylogeographic_runs/All_clades.csv", head=T)
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		tab1 = read.table(paste0(directory,analysis,"_median_weighted_branch_dispersal_velocity.txt"), header=T)
		tab2 = read.table(paste0(directory,analysis,"_95%HPD_weighted_branch_dispersal_velocity.txt"), header=T)
		tab1[,2] = tab1[,2]/366; tab2[,2:3] = tab2[,2:3]/366 # to have the lineage dispersal velocity in km/day
		col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
		dev.new(width=6, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		plot(tab1[,1], tab1[,3], type="l", axes=F, ann=F, ylim=c(0,700), xlim=c(decimal_date(ymd("2020-02-15")),maxYear), col=NA)
		lines(slicedTimes, numberOfBranchesMedianValue, lwd=1, col="red")
		ats = c(decimal_date(ymd("2020-02-15")), decimal_date(ymd("2020-03-03")), decimal_date(ymd("2020-03-14")),
				decimal_date(ymd("2020-03-18")), decimal_date(ymd("2020-04-03")), decimal_date(ymd("2020-05-03")), maxYear)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("","03-03","14-03","18-03","03-04","03-05",""))
		dev.new(width=6, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,12), xlim=c(decimal_date(ymd("2020-02-15")),maxYear), col=NA)
		slicedTimes = tab1[,1]; numberOfBranchesMedianValue = tab1[,3]
		branchDispersalVelocityMedianValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
		lines(slicedTimes, branchDispersalVelocityMedianValue, lwd=1, col=col1)
		ats = c(decimal_date(ymd("2020-02-15")), decimal_date(ymd("2020-03-03")), decimal_date(ymd("2020-03-14")),
				decimal_date(ymd("2020-03-18")), decimal_date(ymd("2020-04-03")), decimal_date(ymd("2020-05-03")), maxYear)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("","03-03","14-03","18-03","03-04","03-05",""))
		directory = "Phylogenetic_analyses/All_dispersal_statistics/Dispersal_statistics_all_branches_ext1/"
		tab = read.csv(paste0(directory,analysis,"_among_provinces_communes_proportions.csv"), header=T)
		col1a = "#FAA521"; col1b = "#FAA52150"; col2a = "#4676BB"; col2b = "#4676BB50"
		dev.new(width=6, height=2.5); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		plot(tab[,1], tab[,2], type="l", axes=F, ann=F, ylim=c(0,0.40), xlim=c(decimal_date(ymd("2020-02-15")),maxYear), col=NA)
		slicedTimes = tab[,1]; lower_l_1 = tab[,3]; upper_l_1 = tab[,4]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2b, border=0)
		lines(slicedTimes, tab[,2], lwd=1, col=col2a)
		ats = c(decimal_date(ymd("2020-02-15")), decimal_date(ymd("2020-03-03")), decimal_date(ymd("2020-03-14")),
				decimal_date(ymd("2020-03-18")), decimal_date(ymd("2020-04-03")), decimal_date(ymd("2020-05-03")), maxYear)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("","03-03","14-03","18-03","03-04","03-05",""))
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
		plot(log(geoDistances), log(branchDurations), col=cols2, pch=16, cex=0.8, xlim=c(1.5,11.5), ylim=c(-7,3.5), axes=F, ann=F)
		points(log(geoDistances), log(branchDurations), col=cols1, pch=1, cex=0.8)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, at=seq(0,13,1), col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, at=seq(-8,5,1), col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="phylogenetic branch durations (days, log-transformed)", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		title(xlab="geographic distance (km, log-transformed)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
	}

