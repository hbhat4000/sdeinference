Rdtq <- function (thetavec, h, k, M, littlet, init_data)
{
    mylist = .Call("dtqCPP", thetavec, h, k, M, littlet, init_data, PACKAGE = "RgdtqSUM")
    return(mylist)
}
