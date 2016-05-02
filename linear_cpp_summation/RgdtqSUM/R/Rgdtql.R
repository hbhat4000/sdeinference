Rgdtql <- function (thetavec, h, k, M, littlet, init_data)
{
    mylist = .Call("gdtqCPP_linear", thetavec, h, k, M, littlet, init_data, PACKAGE = "RgdtqSUM")
    return(mylist)
}