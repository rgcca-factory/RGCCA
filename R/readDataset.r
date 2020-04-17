readDataset=function(bloc=c("C","I","L","M","T"))
{	# lit le fichier de donnees dans un repertoire qui contient les noms du bloc en csv
	listMat=list()
	for(i in 1:length(bloc))
	{
		listMat[[i]]=read.table(paste(bloc[i],".csv",sep=""),sep=";",row.names=1)
	}
	names(listMat)=bloc
	return(listMat)
}
