#' CAD Object of a Steel Panel
#'
#' This object corresponds to a portion of the theoretical CAD model of a steel
#' panel, stored as a `mesh3d` object.
#'
#' @format A `mesh3d` object (list) with three main components:
#' \describe{
#'   \item{vb}{A 4 x 267 matrix containing vertex coordinates in homogeneous form.}
#'   \item{it}{A 3 x 343 matrix containing triangle indices.}
#'   \item{normals}{A 4 x 267 matrix containing vertex normals.}
#' }
#'
#' @source Theoretical CAD model of a steel panel.
#'
"cad"
