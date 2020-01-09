function () 
{
	library(pheatmap)

	options(stringsAsFactors=FALSE)

	mtab4818=function(){
		data=read.csv("E-MTAB-4818-query-results_mod.fpkms.tsv",sep="\t",header=T)
		ix=grep("Solyc01g094190|Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,3:dim(data_f)[2]]
		data_f=na.omit(data_f)
		for(x in 1:dim(data_f)[2]){
			my_mean=mean(data_f[,x])
			my_sd=mean(data_f[,x])
			data_f[,x]=(data_f[,x]-my_mean)/my_sd
		}
		print(data_f)
		par(mfrow=c(2,2))
		pheatmap(data_f,display_numbers=TRUE)
		pheatmap(cor(data_f),display_numbers=TRUE)
		pheatmap(cor(t(data_f)),display_numbers=TRUE)
	}

	mtab4818()
}
