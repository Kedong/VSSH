proportion_heredity <- function(coef, q){
  # q is the # of main effects
  
  model.vio = numeric(dim(coef)[1])

if (dim(coef)[2] > (q*(q+1)/2)){
 	for (sim in 1:nrow(coef)){
		coef_m = coef[sim,1:q]
		coef_2fi = coef[sim,(q+1):((q+1)*q/2)]
		coef_q = coef[sim,(q*(q+1)/2+1):(q*(q+3)/2)]
		
		matrix4main = matrix(0,q,q)
		k = 1
		for (i in 1:q){
			matrix4main[i,i] = coef_q[i]
			for (j in i:q){
				if (i!=j) {
					matrix4main[j,i] = coef_2fi[k]
					k = k+1
				}
			}
		}
		# main effects that should have appeared - A#
		A1 = c(1:q)[apply(1*(matrix4main!=0),1,sum)!=0]
		A2 = c(1:q)[apply(1*(matrix4main!=0),2,sum)!=0]
		A = union(A1, A2)

		# main effects that do not appear - B#
		B = c(1:q)[coef_m==0]

		# take a intersection of A and B #
		# divide it by the total should have appeared #
		C = intersect(A, B)
		if (length(A) == 0){
			model.vio[sim] = 0
		} else {
			model.vio[sim] = length(C)/length(A)
		}
	}
} else {
 	for (sim in 1:nrow(coef)){
		coef_m = coef[sim,1:q]
		coef_2fi = coef[sim,(q+1):((q+1)*q/2)]
		
		matrix4main = matrix(0,q,q)
		k = 1
		for (i in 1:q){
			for (j in i:q){
				if (i!=j) {
					matrix4main[j,i] = coef_2fi[k]
					k = k+1
				}
			}
		}
		# main effects that should have appeared - A#
		A1 = c(1:q)[apply(1*(matrix4main!=0),1,sum)!=0]
		A2 = c(1:q)[apply(1*(matrix4main!=0),2,sum)!=0]
		A = union(A1, A2)

		# main effects that do not appear - B#
		B = c(1:q)[coef_m==0]

		# take a intersection of A and B #
		# divide it by the total should have appeared #
		C = intersect(A, B)
		if (length(A) == 0){
			model.vio[sim] = 0
		} else {
			model.vio[sim] = length(C)/length(A)
		}
	}
}
return(model.vio)
}