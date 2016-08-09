PDFcheck <- function (thetavec, gammavec, runner, h, k, yM)
{
    mylist = .Call("GCPP", thetavec, gammavec, runner, h, k, yM, PACKAGE = "Rdtq2d")
    return(mylist)
}

