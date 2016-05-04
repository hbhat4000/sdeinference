# This code creates a matrix contains Hermite functions
# evaluated on a grid of points. Hermite functions are stored
# along each column starting from zeroth hermite function.

library('EQL')
myhermfun <- function(numterms,xgrid)
{
	nvec = c(0:(numterms-1))

	myhermpol <- function(n) {
	out = hermite(xgrid, n, prob = FALSE)
	return(out)
	}

	# hermpolymat is a matrix of Hermite polynomials stored along the
	# columns. Along the raws the degree of the oloynomial changes
	hermpolmat = sapply(FUN = myhermpol, X=nvec)

	normalizemvec1 = (2^(nvec)*factorial(nvec)*sqrt(pi))^(-1/2)
	normalizemvec2 = exp(-xgrid^2/2)

	normalizationmat = outer(normalizemvec2,normalizemvec1)
	# hermfunmat is a matrix of Hermite functions
	hermfunmat = normalizationmat*hermpolmat

	return(hermfunmat)
}