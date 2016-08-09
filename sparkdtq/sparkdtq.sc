import breeze.linalg._
import breeze.stats.distributions._
import scala.math
import java.io._
import org.apache.spark.mllib.optimization._

// global parameters
val h: Double = 0.01
val k: Double = math.pow(h,0.75)
val M: Int = (math.Pi/math.pow(k,1.5)).toInt
val grid = k*linspace(-M,M,2*M+1)
val irange = convert(linspace(0,2*M,2*M+1),Int) 
val gamma: Int = 25

def drift(theta: DenseVector[Double], x: Double) = {
    val f = theta(0)*(theta(1) - x)
    f
}

def driftV(theta: DenseVector[Double], x: DenseVector[Double]) = {
    val f = theta(0)*(theta(1) - x)
    f
}

def diffusion(theta: DenseVector[Double], x: Double) = {
    val g = theta(2)
    g
}

def diffusionV(theta: DenseVector[Double], x: DenseVector[Double]) = {
    val g = (DenseVector.ones[Double](x.size))*theta(2)
    g
}

def firststep(theta: DenseVector[Double], y: Double) = {
    val mu = y + h*drift(theta,y)
    val sigma = math.abs(diffusion(theta,y))*math.sqrt(h)
    val g = Gaussian(mu,sigma)
    g.pdf(grid)
}

// pdf with vectors of means and variances, evaluated at one grid pt
def ourgauss(theta: DenseVector[Double], x: Double, y: DenseVector[Double]) = {
    val mu = y + h*driftV(theta,y)
    val thisdiff = diffusionV(theta,y)
    val sigma2 = (thisdiff :* thisdiff)*h 
    val xmu = x - mu
    val xmu2 = xmu :* xmu 
    val denom = math.sqrt(2.0*math.Pi*h)*thisdiff.map(math.abs)
    val exparg = -(xmu2 :/ (2.0*sigma2))
    val pdf = exparg.map(math.exp) :/ denom
    pdf
}

// generate one row of the "K" matrix
def genonerow(theta: DenseVector[Double], i: Int, gamma: Int) = {
    val j = i - M
    val x = j*k
    val y = k*linspace(j-gamma,j+gamma,(2*gamma+1))
    k*ourgauss(theta,x,y)
}

// intelligent subsetting of a vector
def mysubset(x: DenseVector[Double], i: Int) = {
    if (i < 0) 0.0
    else if (i >= x.size) 0.0
    else x(i)
}

// given a vector of length N, create i-th vector of length 2*gamma+1
def gammawindow(x: DenseVector[Double], i: Int, gamma: Int) = {
    val grange = convert(linspace(i-gamma,i+gamma,(2*gamma+1)),Int)
    val outvec = grange.map(mysubset(x,_:Int))
    outvec
}

// given the parameter vector \theta and propagator
// given two data points, (t_j,x_j) and (t_{j+1},x_{j+1}), this function
// starts off with a delta function centered at x_j at time t_j, and
// steps the density forward in time until t_{j+1}, where it is evaluated
// at the point x_{j+1}
def denevolve(propagator: DenseVector[DenseVector[Double]], theta: DenseVector[Double], txvec: DenseVector[(Double,Double)]) = {
    var px = firststep(theta,txvec(0)._2)
    val T: Double = txvec(1)._1 - txvec(0)._1
    val numsteps: Int = (math.ceil(T/h).toInt) - 2
    for (i <- 1 to numsteps) {
      val allwins = irange.map(gammawindow(px,_:Int,gamma))
      px = propagator dot allwins
    }
    val finalrow = k*ourgauss(theta,txvec(1)._2,grid)
    finalrow dot px
}

// this function takes x and theta and computes the log likelihood
def loglik(theta: DenseVector[Double], txvec: DenseVector[(Double,Double)]) = {
    val trange = convert(linspace(0,txvec.size-2,txvec.size-1),Int)
    val tslices = trange.map(x => txvec(x to (x+1))).toArray
    val tslicesRDD = sc.parallelize(tslices)

    // compute the propagator
    val prop = irange.map(genonerow(theta,_:Int,gamma))

    // apply denevolve to each consecutive pair of time series points
    val indivdens = tslicesRDD.map(denevolve(prop,theta,_))
    // val indivdens = tslices.map(denevolve(prop,theta,_))

    // compute overall log likelihood
    sum(DenseVector(indivdens.collect()).map(math.log))
    // sum(DenseVector(indivdens).map(math.log))
}

// this function computes p(y | x, sigma_eps^2)
def filtlik(tyvec: DenseVector[(Double,Double)], txvec: DenseVector[(Double,Double)], sigeps2: Double) = {
    val y = tyvec.map(x => x._2)
    val x = txvec.map(x => x._2)

    val normcon = -0.5*math.log(2.0*math.Pi*sigeps2)
    val xmy = x - y
    val xmy2 = xmy :* xmy
    val normmain = -xmy2/(2*sigeps2)
    
    normcon*normmain.size + sum(normmain)
}

// prior for theta
def thetaprior(theta: DenseVector[Double]) = {
    val prior1 = Gaussian(0.5,1)
    val prior2 = Gaussian(2.0,10.0)
    prior1.logPdf(theta(0)) + prior2.logPdf(theta(1))
}

// prior for sigeps2
def sigeps2prior(sigeps2: Double) = {
    val prior = Exponential(1.0)
    prior.logPdf(sigeps2)
}

// auxiliary function that takes two vectors and creates a vector of tuples
def vec2tuples(v1: DenseVector[Double], v2: DenseVector[Double]) = {
    val bigvec = DenseVector.vertcat(v1,v2)
    val txs = v1.size
    val myrange = convert(linspace(0,txs-1,txs),Int)
    myrange.map( i => (bigvec(i),bigvec(i+txs)) )
}

def fulllik(tyvec: DenseVector[(Double,Double)], txvec: DenseVector[(Double,Double)], theta: DenseVector[Double], sigeps2: Double) = {
    var lik = loglik(theta, txvec)
    lik += filtlik(tyvec, txvec, sigeps2)
    lik += thetaprior(theta)
    lik += sigeps2prior(sigeps2)
    lik
}

val xvecproposal = new Gaussian(0.0,0.02)
val thetaproposal = new Gaussian(0.0,0.05)
val sigeps2proposal = new Gaussian(0.0,0.02)
val metro = new Uniform(0,1)

// generate one metropolis sample
// must pass in tyvec, the data
// and also the previous iteration's values for txvec, theta, sigeps2
// and the old likelihood 
// def metropolis(tyvec: DenseVector[(Double,Double)], txvec: DenseVector[(Double,Double)], theta: DenseVector[Double], sigeps2: Double, oldlik: Double) = {

//     // create proposal
//     val xvec = txvec.map(x => x._2)
//     val xvecstar = xvec + DenseVector(xvecproposal.sample(xvec.size).toArray)
//     val txvecstar = vec2tuples(txvec.map(x=>x._1),xvecstar)
//     var thetastar = DenseVector[Double](theta(0),theta(1),math.log(theta(2)))
//     thetastar = thetastar + DenseVector(thetaproposal.sample(theta.size).toArray)
//     thetastar(2) = 0.25
//     val sigeps2star = math.exp(math.log(sigeps2) + sigeps2proposal.sample(1)(0))

//     // evaluate likelihood
//     val likstar = fulllik(tyvec, txvecstar, thetastar, sigeps2star)

//     // accept/reject step
//     val u = metro.sample(1)(0)
//     val ratio = math.exp(likstar - oldlik)
//     if (ratio > u)
//         (1,txvecstar,thetastar,sigeps2star,likstar)
//     else
//         (0,txvec,theta,sigeps2,oldlik)
// }

// aux function to output vector to file
def outvec(v: DenseVector[Double], fname: String) = {
    val pw = new PrintWriter(new FileOutputStream(new File(fname),true))
    pw.write(v.foldLeft("")((a,b) => a+b.toString+',').toString.stripSuffix(","))
    pw.write("\n")
    pw.close
    0
}

// mcmc loop
// def mcmc(timeseries: DenseVector[(Double,Double)], numsamples: Int) = {
//     // initial value of theta
//     var theta: DenseVector[Double] = DenseVector(1.0,0.1,0.25)

//     // initial value of txvec
//     var txvec: DenseVector[(Double,Double)] = timeseries.copy

//     // initial value of sigeps2
//     var sigeps2: Double = 1.0

//     // compute initial likelihood
//     var lik = fulllik(timeseries, txvec, theta, sigeps2)

//     // initial value of accept/ratio flag
//     var accept = 1
   
//     for (i <- 1 to numsamples) {
//         // concatenate everything and save to disk
//         var everything = DenseVector.vertcat(theta,txvec.map(x => x._2))
//         everything = DenseVector.vertcat(everything,DenseVector(sigeps2))
//         everything = DenseVector.vertcat(everything,DenseVector(lik))
//         everything = DenseVector.vertcat(everything,DenseVector(accept))
//         val tmp = outvec(everything,"mcmc.out")

//         // take a metropolis step
//         val metrostep = metropolis(timeseries, txvec, theta, sigeps2, lik)
//         accept = metrostep._1
//         if (accept == 1) {
//             txvec = metrostep._2
//             theta = metrostep._3
//             sigeps2 = metrostep._4
//             lik = metrostep._5
//         }
//     }
//     0
// }

// fake data created in R
import scala.io.Source
val tvecarr = Source.fromFile("tvec.csv").getLines.map(_.split(",")).toArray
val t = DenseVector[Double](tvecarr(0).map(_.toDouble))
val yvecarr = Source.fromFile("xvec.csv").getLines.map(_.split(",")).toArray
val y = DenseVector[Double](yvecarr(0).map(_.toDouble))

val timeseries = vec2tuples(t,y)

// run mcmc
val mcmcout = mcmc(timeseries,1000)

System.exit(0)



