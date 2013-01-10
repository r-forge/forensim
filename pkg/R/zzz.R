
.onAttach <- function(libname, pkgname){
pkg.version <- packageDescription("forensim", fields = "Version")
m0<-"## forensim 3.2.1 is loaded ###"
  packageStartupMessage(m0)
}