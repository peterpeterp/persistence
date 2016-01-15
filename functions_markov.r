library(RNetCDF)

markov_chain_estimation <- function(dataset="_TMean",trendID="91_5",trend_style="_mean",additional_style="",yearPeriod=c(1950,2014),ID_select=c(488),seasons=1:5,orderTot=3){

    print(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind_seasonal_median",".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")


    seasonStart=c(60,152,244,335,1)
    seasonStop=c(151,243,334,424,365)

    for (sea in seasons){
        # take days in one season
        if (sea!=4){
            indSea=ind[,seasonStart[sea]:seasonStop[sea],]
        }
        if (sea==4){
            indSea=c(ind[,seasonStart[sea]:365,],ind[,1:(seasonStop[sea]-365),])
        }
        # take chosen years
        indSea=indSea[,,(yearPeriod[1]-1949):(yearPeriod[2]-1949)]

        markov=array(NA,c(length(ID_select),orderTot*orderTot))
        for (q in ID_select){
            indLoc=as.vector(indSea[q,,])
            print(length(indLoc))


            for (state in c(-1,1)){
                len=length(indLoc)
                for (order in 1:orderTot){
                    resultCond=array(TRUE,len-order)
                    print("------------")
                    print(order)
                    for (o in 2:order){
                        memoryCond=(diff(indLoc,o)==0)[(order-o+1):(len-o)]
                        resultCond=(resultCond==TRUE & memoryCond==TRUE)
                    }
                    outcomeCond=(indLoc==state)[(order+1):len]
                    resultCond=(resultCond==TRUE & outcomeCond==TRUE)

                    stayCond=(diff(indLoc,1)==0)[(order-1+1):(len-1)]
                    stayCond=(resultCond==TRUE & stayCond==TRUE)

                    changeCond=(resultCond==TRUE & stayCond!=TRUE)

                    print(length(which(resultCond)))

                    stays=length(which(stayCond))
                    changes=length(which(changeCond))
                    print(paste(stays,changes))
                    print(stays/(stays+changes))
                }
            }

        }

    }

}

markov_chain_estimation(seasons=5,orderTot=6)