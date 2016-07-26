Rdtq2d <- function (thetavec, gammavec, runner, chaser, h, numsteps, k, yM)
{
    mylist = .Call("dtq2dCPP", thetavec, gammavec, runner, chaser, h, numsteps, k, yM, PACKAGE = "Rdtq2d")
    return(mylist)
}


