function (sel_study) 
{
	library(pheatmap)

	options(stringsAsFactors=FALSE)

	do_norm=function(data_f){
		for(x in 1:dim(data_f)[2]){
			my_mean=mean(data_f[,x])
			my_sd=mean(data_f[,x])
			data_f[,x]=(data_f[,x]-my_mean)/my_sd
		}
		return(data_f)
	}

	mtab4818=function(){
		data=read.csv("E-MTAB-4818-query-results_mod.fpkms.tsv",sep="\t",header=T)
		ix=grep("Solyc01g094190|Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,3:dim(data_f)[2]]
		data_f=na.omit(data_f)
		data_f=do_norm(data_f)
		pheatmap(data_f,display_numbers=TRUE)
		pheatmap(cor(data_f),display_numbers=TRUE)
		pheatmap(cor(t(data_f)),display_numbers=TRUE)
		return(data_f)
	}

	s12864=function(){
		data=read.csv("s12864.txt",sep="\t",dec=",",header=T)
		ix=grep("Solyc01g094190|Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,3:dim(data_f)[2]]
		data_f=na.omit(data_f)
		data_f=do_norm(data_f)
		print(data_f)
		pheatmap(data_f,display_numbers=TRUE)
		pheatmap(cor(data_f),display_numbers=TRUE)
		pheatmap(cor(t(data_f)),display_numbers=TRUE)
		return(data_f)
	}



	s12870=function(sel_picture){
		data=read.csv("study_1.txt",sep="\t",dec=",",row.names=1)
		ix=grep("Solyc01g094190|Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",rownames(data))
		data_f=data[ix,2:7]
		data_f=do_norm(data_f)
		if(sel_picture==1){
			pheatmap(data_f,display_numbers=TRUE)
		}
		else if(sel_picture==2){
			pheatmap(cor(data_f),display_numbers=TRUE)
		}
		else if(sel_picture==3){
			pheatmap(cor(t(data_f)),display_numbers=TRUE)
		}
		return(data_f)
	}
	if(sel_study=="mtab4818"){
		h1=mtab4818()
	}else if(sel_study=="s12870"){
		h1=s12870(1)
	}else if(sel_study=="s12864"){
		h1=s12864()
	}else if(sel_study=="both"){
		h1=mtab4818()
		h2=s12870(1)

		for(x in 1:dim(h1)[1]){
			rn=rownames(h1)[x]
			rn=strsplit(rn,"\\.")[[1]]
			rownames(h1)[x]=rn
		}

		for(x in 1:dim(h2)[1]){
			rn=rownames(h2)[x]
			rn=strsplit(rn,"\\.")[[1]]
			rownames(h2)[x]=rn
		}

		print(h1)
		print(h2)

		uid=unique(rownames(h1),rownames(h2))
		uid=unique(uid)
		print(uid)

		h1=h1[uid,]
		h2=h2[uid,]

		print(dim(h1))
		print(dim(h2))

		df=data.frame(h1)
		df=cbind(df,h2)
		print(df)

		pheatmap(df)
	}
}
