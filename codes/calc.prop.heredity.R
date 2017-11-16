proportion_heredity <- function(x, q){
	# q is the # of main effects

	if (dim(x)[2] > (q*(q+1)/2)){
		main = x[,1:q]
		main_int = t(apply(main, 1, combn, 2, prod))
		interaction = x[,(q+1):(q*(q+1)/2)]
		quad = x[,(q*(q+1)/2+1):(q*(q+3)/2)]
		
		total.nonzero.2fi = numeric(dim(x)[1])
		vio.2fi = numeric(dim(x)[1])
		vio.2fi.ratio = numeric(dim(x)[1])
		total.nonzero.quad = numeric(dim(x)[1])
		vio.quad = numeric(dim(x)[1])
		vio.quad.ratio = numeric(dim(x)[1])
		model.vio = numeric(dim(x)[1])

		for (i in 1:(dim(x)[1])){
			total.nonzero.2fi[i] = sum(1*(interaction[i,]!=0))
			vio.2fi[i] = sum(1*(main_int[i, which(interaction[i,]!=0)]==0))
			vio.2fi.ratio[i] = vio.2fi[i]/total.nonzero.2fi[i]
			total.nonzero.quad[i] = sum(1*(quad[i,]!=0))
			vio.quad[i] = sum(1*(main[i, which(quad[i,]!=0)]==0))
			vio.quad.ratio[i] = vio.quad[i]/total.nonzero.quad[i]
			if (vio.2fi[i] || vio.quad[i]) {
				model.vio[i]=1
			} else {
				model.vio[i]=0
			}
		}
		
		prop.results = cbind(vio.2fi, total.nonzero.2fi, vio.2fi.ratio, vio.quad, total.nonzero.quad, vio.quad.ratio, model.vio)
		return(prop.results)

	} else {

		main = x[,1:q]
		main_int = t(apply(main, 1, combn, 2, prod))
		interaction = x[,(q+1):(q*(q+1)/2)]
		
		total.nonzero.2fi = numeric(dim(x)[1])
		vio.2fi = numeric(dim(x)[1])
		vio.2fi.ratio = numeric(dim(x)[1])
		model.vio = numeric(dim(x)[1])

		for (i in 1:(dim(x)[1])){
			total.nonzero.2fi[i] = sum(1*(interaction[i,]!=0))
			vio.2fi[i] = sum(1*(main_int[i, which(interaction[i,]!=0)]==0))
			vio.2fi.ratio[i] = vio.2fi[i]/total.nonzero.2fi[i]
			if (vio.2fi[i]) {
				model.vio[i]=1
			} else {
				model.vio[i]=0
			}
		}
		prop.results = cbind(vio.2fi, total.nonzero.2fi, vio.2fi.ratio, model.vio)
		return(prop.results)

	}
}