library(RNetCDF)
library(gtools)

markov_chain_estimation <- function(dataset="_TMean",trendID="91_5",trend_style="_mean",additional_style="",yearPeriod=c(1950,2014),ID_select=1:1319,ID_length=1319,seasons=1:5,order=3){

    print(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind_seasonal_median",".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")

    season_names=c("MAM","JJA","SON","DJF","4seasons")
    seasonStart=c(60,152,244,335,1)
    seasonStop=c(151,243,334,424,365)

    # different memory and outcome possibilities
    eventPossibilities=permutations(n=2,r=order,v=c(-1,1),repeats.allowed=TRUE)
    firstOutcome=(0:(2^(order-1)-1)*2+1)
    secondOutcome=(0:(2^(order-1)-1)*2+2)
    combinations=2^order

    # array for end results
    eventResult=array(NA,c(length(season_names),ID_length,4,combinations))
    
    for (sea in seasons){
        # take days in one season
        if (sea<4){
            indSea=ind
            indSea[,1:(seasonStart[sea]-1),]=0
            indSea[,(seasonStop[sea]+1):365,]=0
        }
        if (sea==4){
            indSea=ind
            indSea[,(seasonStop[sea]-365+1):(seasonStart[sea]-1),]=0
        }
        if (sea==5){
            indSea=ind
        }
        # take chosen years
        indSea=indSea[,,(yearPeriod[1]-1949):(yearPeriod[2]-1949)]

        percentage=0
        cat(paste("\n",season_names[sea],"\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage=percentage+5
            }
            indLoc=as.vector(indSea[q,,])
            len=length(indLoc)

            for (combi in 1:combinations){
                resultCond=array(TRUE,len-order)
                for (o in 1:order){
                    memoryCond=(indLoc==eventPossibilities[combi,o])[(order-o+1):(len-o)]
                    resultCond=(resultCond==TRUE & memoryCond==TRUE)
                }

                eventResult[sea,q,1,combi]=length(which(resultCond))

            }
            eventResult[sea,q,2,firstOutcome]=eventResult[sea,q,1,firstOutcome]+eventResult[sea,q,1,secondOutcome]
            eventResult[sea,q,2,secondOutcome]=eventResult[sea,q,1,firstOutcome]+eventResult[sea,q,1,secondOutcome]
            eventResult[sea,q,3,]=eventResult[sea,q,1,]/eventResult[sea,q,2,]

            #print(eventResult[sea,q,2,])
            #print(eventResult[sea,q,3,])

        }

    }

    print(paste("../data/",trendID,"/",dataset,additional_style,"/markov/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,dataset,"_",yearPeriod[1],"-",yearPeriod[2],"_markov_order",order,".nc",sep=""))
    nc_out <- create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/markov/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,dataset,"_",yearPeriod[1],"-",yearPeriod[2],"_markov_order",order,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR","gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "markov chain probabilities")
    att.put.nc(nc_out, "NC_GLOBAL", "markov process order", "NC_CHAR", paste(order-1))
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", paste(yearPeriod[1],"-",yearPeriod[2]))

    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)

    dim.def.nc(nc_out,"outs",dimlength=4,unlim=FALSE)

    dim.def.nc(nc_out,"combinations",dimlength=combinations,unlim=FALSE)
    dim.def.nc(nc_out,"memoryDepth",dimlength=order,unlim=FALSE)

    var.def.nc(nc_out,"eventPossibilities","NC_DOUBLE",c(3,4))
    att.put.nc(nc_out, "eventPossibilities", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "eventPossibilities", "dim_explanation", "NC_CHAR", "combinations-memoryDepth")
    att.put.nc(nc_out, "eventPossibilities", "explanation", "NC_CHAR", "chain of events")

    var.def.nc(nc_out,"eventResult","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "eventResult", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "eventResult", "dim_explanation", "NC_CHAR", "season-ID-outStyle-combinations")
    att.put.nc(nc_out, "eventResult", "explanation", "NC_CHAR", "outStyle: counts of result, counts events with that memory, probability of getting result with memory")
        
    var.put.nc(nc_out,"eventPossibilities",eventPossibilities)      
    var.put.nc(nc_out,"eventResult",eventResult)      
 
    close.nc(nc_out) 


}


markov_chain_estimation(order=4)