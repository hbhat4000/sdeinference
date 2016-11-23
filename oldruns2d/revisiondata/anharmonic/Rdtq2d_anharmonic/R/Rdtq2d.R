Rdtq2d <- function (thetavec, c1, c2, h, numsteps, k, yM)
{
    mylist = .Call("dtq2dCPP", thetavec, c1, c2, h, numsteps, k, yM, PACKAGE = "Rdtq2d")
    return(mylist)
}
