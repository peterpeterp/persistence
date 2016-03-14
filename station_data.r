point_comparison <- function(a,b,center){
    if (a[1] - center[1] >= 0 && b[1] - center[1] < 0){return(TRUE)}
    if (a[1] - center[1] < 0 && b[1] - center[1] >= 0){return(FALSE)}
    if (a[1] - center[1] == 0 && b[1] - center[1] == 0){
        if (a[2] - center[2] >= 0 || b[2] - center[2] >= 0){return(a[2] > b[2])}
        else{return(b[2] > a[2])}      
    }

    #compute the cross product of vectors (center -> a) x (center -> b)
    det = (a[1] - center[1]) * (b[2] - center[2]) - (b[1] - center[1]) * (a[2] - center[2])
    if (det < 0){return(TRUE)}
    if (det > 0){return(FALSE)}

    #points a and b are on the same line from the center
    #check which point is closer to the center
    d1 = (a[1] - center[1]) * (a[1] - center[1]) + (a[2] - center[2]) * (a[2] - center[2])
    d2 = (b[1] - center[1]) * (b[1] - center[1]) + (b[2] - center[2]) * (b[2] - center[2])
    return(d1>d2)
}

sorting <- function(points,center){
    for (i in 1:dim(points)[1]){
        for (j in 1:dim(points)[1]){
            for (k in j:dim(points)[1]){
                if (point_comparison(points[j,],points[k,],center)){
                    lower<-points[k,]
                    higher<-points[j,]
                    points[j,]<-lower
                    points[k,]<-higher
                }
            }
        }
    }
    return(points)
}

region_border_to_file <- function(region_name="ward23"){
    # loads attribution file and visualizes regions on map
    attribution<-read.table(paste("../data/_TMean/ID_regions/",region_name,".txt",sep=""))[,1]
    mids<-read.table(paste("../data/_TMean/ID_regions/",region_name,"_mids.txt",sep=""))    

    regNumb<-dim(mids)[1]
    polygons<-array(NA,c(regNumb,200,2))
    for (G in 1:regNumb){
        inside<-which(attribution==G)
        points<-array(NA,c(length(inside)*10,2))
        index<-0
        for (q in inside){
            print(index)
            if (!((dat$lon[q]+3.75) %in% dat$lon[inside[which(dat$lat[inside]==dat$lat[q])]])){
                points[index<-index+1,]=c(dat$lon[q]+1.875,dat$lat[q]-1.25)
                points[index<-index+1,]=c(dat$lon[q]+1.875,dat$lat[q]+1.25)
            }
            if (!((dat$lon[q]-3.75) %in% dat$lon[inside[which(dat$lat[inside]==dat$lat[q])]])){
                points[index<-index+1,]=c(dat$lon[q]-1.875,dat$lat[q]-1.25)
                points[index<-index+1,]=c(dat$lon[q]-1.875,dat$lat[q]+1.25)             
            }
            if (!((dat$lat[q]+2.5) %in% dat$lat[inside[which(dat$lon[inside]==dat$lon[q])]])){
                points[index<-index+1,]=c(dat$lon[q]-1.875,dat$lat[q]+1.25)
                points[index<-index+1,]=c(dat$lon[q]+1.875,dat$lat[q]+1.25)                 
            }
            if (!((dat$lat[q]-2.5) %in% dat$lat[inside[which(dat$lon[inside]==dat$lon[q])]])){
                points[index<-index+1,]=c(dat$lon[q]-1.875,dat$lat[q]-1.25)
                points[index<-index+1,]=c(dat$lon[q]+1.875,dat$lat[q]-1.25)                 
            }
        }
        polygons[G,1:index,]=sorting(points[1:index,],rowMeans(points[1:index,]))
        polygon(x=polygons[G,,1],y=polygons[G,,2],border="red",col=rgb(0.5,0.5,1,0.5)) 
        polygon(x=points[1:index,1],y=points[1:index,2],border="green",col=rgb(0.5,1,0.5,0.5)) 

        poligon=cbind(x=points[1:index,1],y=points[1:index,2])
        points=cbind(x=dat$lon,y=dat$lat)
        inside=pnt.in.poly(points,poligon)$pip
        points(dat$lon[inside],dat$lat[inside],col="blue")
    }
}