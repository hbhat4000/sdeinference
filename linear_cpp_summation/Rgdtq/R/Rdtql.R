Rdtql <- function (thetavec, h, k, M, littlet, init_data)
{
    mylist = .Call("dtqCPP_linear", thetavec, h, k, M, littlet, init_data, PACKAGE = "Rgdtq")
    return(mylist)
}
