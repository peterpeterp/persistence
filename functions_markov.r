
markov_chain_estimation <- function(dataset="_TMean",trendID="91_5",additional_style="",period="1950-2014",ID_select=488){
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")


}

markov_chain_estimation()