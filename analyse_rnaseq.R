library(pheatmap)
options(stringsAsFactors=FALSE)

analyse_rnaseq=function (sel_study) 
{
	#
	# normalization of the data
	#
	do_norm=function(data_f){
		for(x in 1:dim(data_f)[2]){
			my_mean=mean(data_f[,x])
			my_sd=mean(data_f[,x])
			data_f[,x]=(data_f[,x]-my_mean)/my_sd
		}
		return(data_f)
	}

	#
	#	tomato fruit at four developmental stages
	#
	mtab4818=function(){
		data=read.csv("E-MTAB-4818-query-results_mod.fpkms.tsv",sep="\t",header=T)
		ix=grep("Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,3:dim(data_f)[2]]
		data_f=na.omit(data_f)
		data_f=do_norm(data_f)
		print("tomato fruit at four developmental stages")
		print(colnames(data))
		return(data_f)
	}

	#
	#	root transcriptome in chilling-sensitive tomato (S. lycopersium) c.t. more cold-tolerant (S.moneymaker)
	#
	mtab4855=function(){
		data=read.csv("E-MTAB-4855.txt",dec=",",sep="\t",header=T)
		print(head(data))
		ix=grep("Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,2:dim(data_f)[2]]
		data_f=na.omit(data_f)
		data_f=do_norm(data_f)
		print("root transcriptome in chilling-sensitive tomato (S. lycopersium) c.t. more cold-tolerant (S.moneymaker)")
		print(colnames(data))
		return(data_f)
	}

	#
	#	role of red light in resistance to P. syringae
	#
	s12864=function(){
		data=read.csv("s12864.txt",sep="\t",dec=",",header=T)
		ix=grep("Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",data[,1])
		data_f=data[ix,]
		rownames(data_f)=data_f[,1]
		data_f=data_f[,3:dim(data_f)[2]]
		data_f=na.omit(data_f)
		data_f=do_norm(data_f)
		print("role of red light in resistance to P. syringae")
		print(colnames(data))
		return(data_f)
	}


	#
	#	low temperature that differ in freezing S. lycopersium / S. habrochaites
	#
	s12870=function(sel_picture){
		data=read.csv("study_1.txt",sep="\t",dec=",",row.names=1)
		ix=grep("Solyc01g101060|Solyc09g008280|Solyc10g083970|Solyc12g099000",rownames(data))
		data_f=data[ix,2:7]
		data_f=do_norm(data_f)
		print("low temperature that differ in freezing S. lycopersium / S. habrochaites")
		print(colnames(data))
		return(data_f)
	}
	if(sel_study=="mtab4818"){
		h1=mtab4818()
	}else if(sel_study=="s12870"){
		h1=s12870(1)
	}else if(sel_study=="s12864"){
		h1=s12864()
	}else if(sel_study=="mtab4855"){
		h1=mtab4855()
	}else if(sel_study=="p1"){
	
		h1=mtab4818()
		h2=s12870(1)
		h4=mtab4855()

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

		for(x in 1:dim(h4)[1]){
			rn=rownames(h4)[x]
			rn=strsplit(rn,"\\.")[[1]][1]
			rownames(h4)[x]=rn
		}

		uid=unique(c(rownames(h1),rownames(h2),rownames(h4)))
		uid=unique(uid)

		h1=h1[uid,]
		h2=h2[uid,]

		df=data.frame(h1)
		df=cbind(df,h2)
		df=cbind(df,h4)

		pheatmap(df,cluster_cols=FALSE)

	}else if(sel_study=="p2"){
	
		h3=s12864()

		for(x in 1:dim(h3)[1]){
			rn=rownames(h3)[x]
			rn=strsplit(rn,"gene:")[[1]][2]
			rownames(h3)[x]=rn
		}

		uid=unique(rownames(h3))
		uid=unique(uid)

		h3=h3[uid,]

		df=data.frame(h3)

		pheatmap(df,cluster_cols=FALSE)
	}
	
}

analyse_rnaseq("p1")

