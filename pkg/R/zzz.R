
.First.lib <-function (lib, pkg) {
library.dynam("forensim", pkg, lib)
m0<-"## forensim 3.2.1 is loaded ###"
packageStartupMessage(m0)
}
