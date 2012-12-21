createSplsTuneGrid<-function (Kmin=1, Kmax=10, Ks = NULL,etamin =0.1, etamax = 0.99, etas = NULL, kappas = NULL, 
          kappamin = 0.1, verbose = TRUE) 
{
  if(is.null(Ks)){
    Ks <- seq(1:Kmax)
  }
  if (is.null(kappas))
    kappas <- c(0.0001,0.001,0.1,0.5)
  if (is.null(etas)) {
    etas <- c(0.01,seq(etamin,0.99, by = 0.1),0.99)
  }
  message(paste("Testing", length(Ks), "values of K(hidden components) from", 
                round(min(Ks), 5), "to", round(max(Ks), 5)))
  message(paste("Testing", length(kappas), "values of kappa from", 
                round(min(kappas), 5), "to", round(max(kappas), 5)))
  message(paste("Testing", length(etas), "values of eta from", 
                round(min(etas), 3), "to", round(max(etas), 3)))
  Grid <- expand.grid(.K = Ks,.kappa = kappas, .eta = etas)
  Grid
}
