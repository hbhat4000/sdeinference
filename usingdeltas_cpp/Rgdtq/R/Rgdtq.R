Rgdtq <- function (thetavec, c0, h, numsteps, k, yM)
{
    mylist = .Call("gdtqCPP", thetavec, c0, h, numsteps, k, yM, PACKAGE = "Rgdtq")
    return(mylist)
}


