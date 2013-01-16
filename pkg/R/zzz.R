
.onAttach <- function(libname, pkgname){
pkgversion <- packageDescription("forensim", fields = "Version")
m0<-paste("## forensim ",pkgversion," is loaded ###")
  packageStartupMessage(m0)
}

