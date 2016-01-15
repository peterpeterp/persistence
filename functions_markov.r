library(RNetCDF)
library(gtools)

markov_chain_estimation <- function(dataset="_TMean",trendID="91_5",trend_style="_mean",additional_style="",yearPeriod=c(1950,2014),ID_select=c(488),ID_length=1319,seasons=1:5,orderTot=3){

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

        combinations=2^orderTot
        eventCounts=array(NA,c(ID_length,orderTot,combinations))
        memory=array(NA,c(orderTot,combinations,orderTot))
        for (order in 1:orderTot){memory[order,1:2^order,1:order]=permutations(n=2,r=order,v=c(-1,1),repeats.allowed=TRUE)}



        for (q in ID_select){
            indLoc=as.vector(indSea[q,,])
            #print(length(indLoc))

            len=length(indLoc)
            for (order in 1:orderTot){
                for (combi in 1:combinations){
                    resultCond=array(TRUE,len-order)
                    #print(memory[order,combi,])
                    for (o in 1:order){
                        #memoryCond=(diff(indLoc,o)==0)[(order-o+1):(len-o)]
                        #print(paste(order,combi,o))
                        
                        memoryCond=(indLoc==memory[order,combi,o])[(order-o+1):(len-o)]
                        #print(paste(memory[order,combi,o]))
                        #print(length(which(memoryCond)))
                        resultCond=(resultCond==TRUE & memoryCond==TRUE)
                        #print(length(which(resultCond)))
                    }
                    #outcomeCond=(indLoc==state)[(order+1):len]
                    #stayProb=length(which(resultCond==TRUE & outcomeCond==TRUE))

                    #resultProb=length(which(resultCond))

                    #outcomeCond=(indLoc!=state)[(order+1):len]
                    #changeProb=length(which(resultCond==TRUE & outcomeCond==TRUE))

                    #print(resultProb)
                    #print(paste(stayProb,changeProb))
                    #print(stayProb/resultProb)
                    
                    #print(length(which(resultCond)))
                    eventCounts[q,order,combi]=length(which(resultCond))
                }
                
            }
            for (order in 1:orderTot){
                print(memory[order,1:(2^order),1:order])
                print(eventCounts[q,order,1:(2^order)])

                print(0:(2^(order-1)-1)*2)

                print(markov[q,order,(0:(2^(order-1)-1)*2+1)])
                print(markov[q,order,(0:(2^(order-1)-1)*2+2)])
                

                eventCounts[order,]=eventCounts[q,order,(0:(2^(order-1)-1)*2+1)]/(eventCounts[q,order,(0:(2^(order-1)-1)*2+1)]+eventCounts[q,order,(0:(2^(order-1)-1)*2+2)])

                print(sum(markov[q,order,1:(2^order)]))
            }

        }

    }

}



markov_chain_estimation(seasons=2,orderTot=4)