PDFcheck <- function (thetavec, h, k, yM)
{
    mylist = .Call("GCPP", thetavec, h, k, yM, PACKAGE = "Rdtq2d")
    return(mylist)
}