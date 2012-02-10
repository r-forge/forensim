.First.lib <- function (lib, pkg){
    library.dynam("forensim", pkg, lib)
    pkg.version <- packageDescription("forensim", fields = "Version")

    startup.txt <- paste("   ==========================\n    adegenet", pkg.version, "is loaded\n   ==========================\n\n - to start, type '?adegenet'\n - to browse adegenet website, type 'adegenetWeb()'\n - to post questions/comments: adegenet-forum@lists.r-forge.r-project.org\n\n")

    packageStartupMessage(startup.txt)
}