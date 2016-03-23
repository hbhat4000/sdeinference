PDFcheck <- function (thetavec, h, k, yM)
{
    mylist = .Call("GCPP", thetavec, h, k, yM, PACKAGE = "Rgdtq")
    return(mylist)
}

