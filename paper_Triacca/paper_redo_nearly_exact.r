# teste teste


source("../functions_persistence.r")



if (2==2){

	if (8==9) {
		laenge<-250
		t<-seq(1,laenge,1)
		x<-seq(1,laenge,1)
		y<-seq(1,laenge,1)
		ut<-seq(1,laenge,1)

		for (i in 0:laenge){
			ut[i]=runif(1,-1,1)
		}


		for (i in 0:laenge){
			if(i<3){
				y[i]=ut[i]
				x[i]=ut[i]
			}
			if(i>2){
				x[i]=1.7*x[i-1]-0.7*x[i-2]+ut[i]
				y[i]=0.3*y[i-1]+0.7*y[i-2]+ut[i]
			}
			
		}
		first_dif=x

	}


	if (2==2) {
		tmp = read.table("land_anomalie.txt",sep="\t")
		data=tmp[,2]
		time=tmp[,1]
		size=length(data)


		first_dif=data*NA
		for (i in 2:size){
			first_dif[i]=data[i]-data[i-1]
		}

		pdf(file="land")
		plot(time,data,col="white")
		lines(time,data)
		#lines(time,first_dif,col="red")
		graphics.off()

		first_dif=diff(data)
		tmp=shock_ar(as.vector(first_dif),2)
		print(tmp)

	}

	if (7==7){
		shock=array(NA,10)
		bic=array(NA,10)
		for (i in 1:5){
			time0=proc.time()[1]
			tmp=shock_ar(as.vector(data),i)
			shock[i+5]=tmp$P_w
			bic[i+5]=tmp$bic
			print(proc.time()[1]-time0)
			cat("\n ar ",i, " bic = ",bic[i+5], " per = ",shock[i+5],"\n")

		}		

		for (i in 1:5){
			time0=proc.time()[1]
			tmp=shock_ma(as.vector(data),i)
			shock[i]=tmp$P_w
			bic[i]=tmp$bic
			print(proc.time()[1]-time0)
			cat("\n ma ",i, " bic = ",bic[i], " per = ",shock[i],"\n")
		}		


		print(shock)
		print(bic)

	}
}


