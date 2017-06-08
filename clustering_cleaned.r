
crossCor <- function(x,y,lagMax){
    tmp<-ccf(x,y,lag.max=lagMax,plot=FALSE,na.action=na.pass)
    return(c(max(abs(tmp$acf)),tmp$lag[which.max(abs(tmp$acf))]))
}


dissimilarity_matrix <- function(X,lagMax=3,normalize=FALSE,ignore_threshold=1,out_file=paste('correlation_distance_',lagMax,'.nc',sep='')){
    # X should have the following shape: (ID,time) NOT (lon,lat,time)

    ID_len=dim(X)[1]
    time_Len=dim(X)[2]
    
    # normalize
    if (normalize){
        for (q in 1:ID_len){
            if (length(which(!is.na(X[q,])))>time_Len/2){
                X[q,]=X[q,]/sd(X[q,],na.rm=TRUE)
            }
            else{X[q,]=NA}
        }
        print(proc.time())
    }

    # to many nas -> flat 
    X[is.na(X)]=0

    # cross correlation
    distance_mat=array(NA,c(ID_len,ID_len))
    lag_choice_mat=array(NA,c(ID_len,ID_len))
    for (i in 1:ID_len){
        # das ist mit taeglichen Daten echt langsam, deswegen wollte ich mir das anzeigen lassen
        print(i)
        print(proc.time())
        # hier könnte man das recht einfach optimieren, weil einiges doppelt ausgerechnet wird
        for (j in 1:ID_len){
            if (length(which(!is.na(X[i,])))>ignore_threshold & length(which(!is.na(X[j,])))>ignore_threshold){
                tmp=crossCor(X[i,],X[j,],lagMax=lagMax)
                distance_mat[i,j]=1-abs(tmp[1])
                lag_choice_mat[i,j]=tmp[2]
            }
        }
    }

    # write results to netcdf
    nc_out<-create.nc(out_file)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "distance matrix for time series")
    att.put.nc(nc_out, "NC_GLOBAL", "max lags in both directions", "NC_CHAR", paste(lagMax))
    
    dim.def.nc(nc_out,"ID",dimlength=ID_len, unlim=FALSE)   
    dim.def.nc(nc_out,"ID2",dimlength=ID_len, unlim=FALSE)   

    var.def.nc(nc_out,"distance_mat","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "distance_mat", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "distance_mat", "dim_explanation", "NC_CHAR", "ID-ID")
    att.put.nc(nc_out, "distance_mat", "explanation", "NC_CHAR", "distences between IDs 1 time lag allowed")

    var.def.nc(nc_out,"lag_choice_mat","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "lag_choice_mat", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "lag_choice_mat", "dim_explanation", "NC_CHAR", "ID-ID")
    att.put.nc(nc_out, "lag_choice_mat", "explanation", "NC_CHAR", "which lag is closest")

    var.put.nc(nc_out,"distance_mat",distance_mat)      
    var.put.nc(nc_out,"lag_choice_mat",lag_choice_mat)      
 
    close.nc(nc_out) 

    return(distance_mat)    
}


cluster_evaluation <- function(distance_mat,method="ward.D2",untilGroup=11,plot_file=paste(method,'_clustering.pdf',sep=''),out_file=paste(method,'_clustering.nc',sep='')){
    distance_mat[is.na(distance_mat)]=1
    ID_len=dim(distance_mat)[1]


    # to many nas -> flat 
    ID_select_not_flat<-which(!is.na(distance_mat[1,]))

    attribution=array(NA,c(untilGroup,ID_len))

    if (method!="kmeans"){
        pdf(file=plot_file,width=10,height=3)
        par(mar=c(0,5,0,0))
        tmp<<-hclust(as.dist(distance_mat[ID_select_not_flat,ID_select_not_flat]),method=method)
        for (nGroup in 1:untilGroup){
            print(nGroup)
            attribution[nGroup,ID_select_not_flat]=cutree(tmp,nGroup)
            plot(tmp,main='',xlab='')  
            if (nGroup>=2){
                rect.hclust(tmp, k = nGroup, border = "red")  
            }        
        }
        # clustering height
        tmp_length<-length(tmp$height)
        clustering_height=tmp$height[(tmp_length):(tmp_length-untilGroup+1)]
    }
    graphics.off()

    if (method=="kmeans"){
        attribution[1,ID_select_not_flat]=array(1,length(ID_select_not_flat))
        for (nGroup in 2:untilGroup){
            tmp<<-pam(x=as.dist(distance_mat[ID_select_not_flat,ID_select_not_flat]),k=nGroup,diss=TRUE)
            attribution[nGroup,ID_select_not_flat]=tmp$clustering
        }
    }

    # within sum of dissimilarities
    within_sum_dissimilarity<-array(0,untilGroup)
    for (k in 1:untilGroup){
        for (G in unique(attribution[k,!is.na(attribution[k,])])){
            inGroup<-which(attribution[k,]==G)
            within_sum_dissimilarity[k]<-within_sum_dissimilarity[k]+sum(distance_mat[inGroup,inGroup]^2,na.rm=TRUE)
        }
    }

    # write results to netcdf
    nc_out<-create.nc(out_file)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "different cuttrees with indices etc")
    att.put.nc(nc_out, "NC_GLOBAL", "method", "NC_CHAR", method)
    
    dim.def.nc(nc_out,"ID",dimlength=ID_len, unlim=FALSE)   
    dim.def.nc(nc_out,"clusters",dimlength=untilGroup, unlim=FALSE)   

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(1,0))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "nClusters-ID")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "attributions")

    var.def.nc(nc_out,"clustering_height","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "clustering_height", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "clustering_height", "dim_explanation", "NC_CHAR", "nClusters-indices")

    var.def.nc(nc_out,"within_sum_dissimilarity","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "within_sum_dissimilarity", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "within_sum_dissimilarity", "dim_explanation", "NC_CHAR", "nClusters-indices")

    var.put.nc(nc_out,"attribution",attribution)      
    var.put.nc(nc_out,"clustering_height",clustering_height)      
    var.put.nc(nc_out,"within_sum_dissimilarity",within_sum_dissimilarity)      
 
    close.nc(nc_out)         
}

# load libraries
library(RNetCDF)
library(cluster)

# load data
nc  = open.nc('/home/peter/Dokumente/pik/backuped/data/_TMean/HadGHCND_TMean_data3D.day1-365.1950-2014.nc')
tas = var.get.nc(nc,"tas")
print(dim(tas))
# reshape
tas<-array(tas,dim=c(1319,365*65))

# compute dissimilarity matrix (die Funktion ist recht langsam)
distance_matrix=dissimilarity_matrix(tas[1:100,1:1000],lagMax=1,normalize=TRUE)

# start clustering
cluster_evaluation(distance_matrix)

