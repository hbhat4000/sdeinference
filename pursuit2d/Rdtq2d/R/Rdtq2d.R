Rdtq2d <- function (thetavec, xc, yc, xr, yr, h, numsteps, k, yM)
{
    mylist = .Call("dtq2dCPP", thetavec, xc, yc, xr, yr, h, numsteps, k, yM, PACKAGE = "Rdtq2d")
    return(mylist)
}


