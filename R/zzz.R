.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)

  packageStartupMessage(
    paste0(
      "Package ", pkgname, ": Tools for Geometric Processing and ",
      "Dimensional Analysis,\n",
      "version ", version, "\n",
      "Developed by the MODES-CEMI group.\n",
      "Type 'help(package = \"", pkgname, "\")' for an overview of the package or\n",
      "visit https://modes-cemi.github.io/dimControl/."
    )
  )
}
