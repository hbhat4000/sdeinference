Rgdtql <- function (thetavec, h, k, M, littlet, init_data)
{
    mylist = .Call("gdtqCPP", thetavec, h, k, M, littlet, init_data, PACKAGE = "Rgdtql")
    return(mylist)
}


