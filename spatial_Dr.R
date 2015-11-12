spatial_Dr <- function(data,sp_list,x_col,y_col,x.range,y.range,nsample=F){
	Dws <- unlist(lapply(sp_list,D_within,data=data,x_col=x_col,y_col=y_col,x.range=x.range,y.range=y.range,nsample=nsample))
	Dws	
}

D_within <- function(data,sp,x_col,y_col,x.range,y.range,nsample=F){
	data_sp <- data[data$sp==sp,]
	if(nsample){
		data_sp <- data_sp[sample(1:dim(data_sp)[1],nsample),]
	}
	#MPD(data_sp[,c(x_col,y_col)])
	if(sum(as.numeric(dist(data_sp[,c(x_col,y_col)]))==0)){
		print(paste(as.character(sp),"duplicate",sum(as.numeric(dist(data_sp[,c(x_col,y_col)]))==0)))
	}
	#sp_pattern <- ppp(data_sp[,x_col], data_sp[,y_col], x.range,y.range) ##range flexible for each species
	sp_pattern <- ppp(data_sp[,x_col], data_sp[,y_col], range(data_sp[,x_col]),range(data_sp[,y_col])) ##range flexible for each species
	clarkevans(sp_pattern,"donnelly")
}

MPD <- function(xy_mat){
	PDs <- as.numeric(dist(xy_mat))
	sd(PDs)/mean(PDs)
}

sp <- sps[1]
data_sp <- data[data$sp==sp,]
sp_pattern <- ppp(data_sp[,3], data_sp[,4], c(0,1000),c(0,900))
clarkevans(sp_pattern,"cdf")