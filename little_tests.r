

transistion_probabilities_mean_duration <- function(p=0.77){
	for (j in c(10,20,30,40,50,100,1000)){
		expt_d=0
		for (i in 1:j){
			expt_d=expt_d+i*(p^(i-1)*(1-p))
		}
		print(expt_d)
	}
}

transistion_probabilities_mean_duration()