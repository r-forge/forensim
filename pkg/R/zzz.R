
.First.lib <-function (lib, pkg) {
library.dynam("forensim", pkg, lib)
packageStartupMessage('## forensim 3.0 is loaded ###' )
}