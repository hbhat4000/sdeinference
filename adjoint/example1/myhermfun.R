# This code creates a matrix containing Hermite functions
# evaluated on a grid of points. Hermite functions are stored
# along each column starting from zeroth Hermite function.

library('EQL')
myhermfun <- function(numterms, xgrid)
{
	nvec = c(0:(numterms-1))

	myhermpol <- function(n) 
	{
		out = hermite(xgrid, n, prob = FALSE)
		return(out)
	}

	# hermpolmat is a matrix of Hermite polynomials stored along the columns.
	# Along the rows the degree of the Hermite function increases starting from zero th degree to length(c0)-1
	hermpolmat = sapply(FUN = myhermpol, X = nvec)

	# Normalized Hermite polynomial
	# phi(x) = exp(-(x^2)/2) / (sqrt(n! 2^n sqrt(pi))) * H_n(x)
	normalizemvec1 = (2^(nvec)*factorial(nvec)*sqrt(pi))^(-1/2)
	normalizemvec2 = exp(-(xgrid^2) / 2)

	normalizationmat = outer(normalizemvec2, normalizemvec1)
	# hermfunmat is a matrix of Hermite functions
	hermfunmat = normalizationmat * hermpolmat

	return(hermfunmat)
}