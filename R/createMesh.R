#' Reconstruct a Triangular Mesh from a 3D Point Cloud
#'
#' Reconstructs a triangular surface from a 3D point cloud and returns it as a
#' `mesh3d` object. The function performs linear binning using `npsp::binning()`,
#' optionally adds a zero-valued border around the point cloud in the X and Z
#' directions, removes low-density nodes, truncates high-density nodes and
#' extracts the isosurface using `misc3d::contour3d()`.
#'
#' Disconnected surface components are sorted by size and progressively merged
#' until the reconstructed surface reaches a specified proportion of a reference
#' bounding box. If `mesh_cad` is provided, the CAD model dimensions are used as
#' reference; otherwise, the dimensions of the input point cloud are used.
#'
#' @param panel A numeric matrix or data frame with three columns containing the X,
#' Y and Z coordinates of the point cloud.
#' @param resol_nbin A numeric vector of length three defining the approximate
#' binning resolution in the X, Y and Z directions. Default is `c(4, 4, 2)`.
#' @param mesh_cad Optional. A `mesh3d` object used as a reference for determining
#' the minimum spatial coverage of the reconstructed surface. If `NULL`, the
#' dimensions of the input point cloud are used instead. Default is `NULL`.
#' @param zero_border Logical. If `TRUE`, one zero-valued layer is added on the
#' lower and upper sides of the X and Z directions before extracting the
#' isosurface. If `FALSE`, no external zero-valued layers are added. Default is
#' `FALSE`.
#' @param k_factor Numeric. Factor used to compute the tolerance for removing
#' low-density nodes. If `NULL`, it is automatically assigned according to
#' `resol_nbin`. Default is `NULL`.
#' @param trunc_factor Numeric. Factor used to truncate high-density node weights.
#' Default is `3.5`.
#' @param level_factor Numeric. Factor used to define the isosurface extraction
#' level from the truncated weights. Default is `0.99`.
#' @param coverage_tol Numeric. Minimum proportion of the reference bounding box
#' dimensions that must be covered by the reconstructed surface. The reference is
#' taken from `mesh_cad` when provided, or from the input point cloud otherwise.
#' Default is `0.99`.
#'
#' @returns A list containing:
#' \itemize{
#'   \item `mesh`: the reconstructed surface as an object of class `mesh3d`.
#'   \item `bin`: the processed binning object.
#'   \item `contours`: the disconnected triangular surfaces returned by
#'   `misc3d::contour3d()`, ordered from largest to smallest.
#'   \item `params`: a list containing the parameters used during the
#'   reconstruction process.
#' }
#'
#' @details
#' Linear binning is computed using `npsp::binning()`, and the corresponding grid
#' coordinates are obtained with `npsp::coordvalues()`. The isosurface is
#' extracted using `misc3d::contour3d()`, which implements the Marching Cubes
#' algorithm.
#'
#' When `zero_border = TRUE`, one external zero-valued layer is added on both
#' sides of the X axis and both sides of the Z axis. This creates a transition
#' between occupied and empty regions and can facilitate surface closure during
#' isosurface extraction.
#'
#' Low-density nodes are removed using a tolerance based on the mean number of
#' points per occupied bin. High-density node weights are subsequently truncated,
#' and the Marching Cubes extraction level is defined from the truncated value.
#'
#' When several disconnected surfaces are generated, they are ordered from
#' largest to smallest and merged until the reconstructed surface reaches
#' `coverage_tol` times the dimensions of the selected reference bounding box.
#'
#' @seealso `npsp::binning()`, `npsp::coordvalues()`,
#' `misc3d::contour3d()`, `rgl::tmesh3d()`
#'
#' @examples
#' \dontrun{
#' # Reconstruct a mesh without adding zero-valued border layers
#' result <- createMesh(
#'   panel = panel,
#'   resol_nbin = c(4, 4, 2),
#'   zero_border = FALSE
#' )
#' mesh <- result$mesh
#'
#' # Reconstruct a mesh adding one zero-valued layer on the X and Z borders
#' result <- createMesh(
#'   panel = panel,
#'   resol_nbin = c(4, 4, 2),
#'   zero_border = TRUE
#' )
#' mesh <- result$mesh
#'
#' # Use a CAD mesh as reference for the coverage criterion
#' result <- createMesh(
#'   panel = panel,
#'   resol_nbin = c(4, 4, 2),
#'   mesh_cad = mesh_cad,
#'   zero_border = TRUE
#' )
#' mesh <- result$mesh
#'
#' rgl::open3d()
#' rgl::shade3d(mesh, col = "lightgray")
#' rgl::decorate3d()
#' }
#'
#' @importFrom npsp binning coordvalues
#' @importFrom rgl tmesh3d
#' @importFrom utils capture.output getFromNamespace
#'
#' @export
createMesh <- function(panel,
                               resol_nbin = c(4, 4, 2),
                               mesh_cad = NULL,
                               zero_border = FALSE,
                               k_factor = NULL,
                               trunc_factor = 3.5,
                               level_factor = 0.99,
                               coverage_tol = 0.99) {

  panel <- as.matrix(panel)

  if (!is.numeric(panel))
    stop("Argument 'panel' must be numeric")

  if (ncol(panel) != 3)
    stop("Argument 'panel' must have three columns: x, y and z")

  if (anyNA(panel))
    stop("Argument 'panel' contains NA values")

  if (length(resol_nbin) != 3 || any(resol_nbin <= 0))
    stop("Argument 'resol_nbin' must be a positive numeric vector of length 3")

  if (!is.logical(zero_border) || length(zero_border) != 1)
    stop("Argument 'zero_border' must be TRUE or FALSE")

  if (coverage_tol <= 0 || coverage_tol > 1)
    stop("Argument 'coverage_tol' must be greater than 0 and less than or equal to 1")

  # Point cloud dimensions
  dimlen <- apply(panel, 2, function(x) diff(range(x)))

  # Number of bins in each spatial direction
  nbin <- trunc(dimlen / resol_nbin)
  nbin <- pmax(nbin, 2)

  # Linear binning
  bin <- npsp::binning(
    x = panel,
    nbin = nbin,
    type = "linear"
  )
  bin$data <- NULL

  # Add one external zero-valued layer on the four sides:
  # left/right (X) and bottom/top (Z)
  if (zero_border) {

    dims <- bin$grid$n
    new_dims <- dims + c(2, 0, 2)

    binw_ext <- array(
      0,
      dim = new_dims
    )

    binw_ext[
      2:(dims[1] + 1),
      ,
      2:(dims[3] + 1)
    ] <- bin$binw

    bin$binw <- binw_ext
    bin$grid$n <- new_dims

    # Extend X limits
    bin$grid$min[1] <- bin$grid$min[1] - bin$grid$lag[1]
    bin$grid$max[1] <- bin$grid$max[1] + bin$grid$lag[1]

    # Extend Z limits
    bin$grid$min[3] <- bin$grid$min[3] - bin$grid$lag[3]
    bin$grid$max[3] <- bin$grid$max[3] + bin$grid$lag[3]
  }

  bin$grid$dimnames <- c("X", "Y", "Z")

  # Automatic tolerance factor
  if (is.null(k_factor)) {

    if (all(resol_nbin == c(4, 4, 4))) {

      k_factor <- 0.1

    } else if (
      all(resol_nbin == c(2, 2, 2)) ||
      all(resol_nbin == c(4, 4, 2))
    ) {

      k_factor <- 0.05

    } else {

      k_factor <- 0.05

      warning(
        "No specific rule for 'resol_nbin'. Using k_factor = 0.05"
      )
    }
  }

  # Non-empty binning nodes
  w <- bin$binw[bin$binw > 0]

  if (length(w) == 0)
    stop("The binning object has no positive weights")

  # Remove low weights
  tol <- (nrow(panel) / length(w)) * k_factor
  bin$binw[bin$binw <= tol] <- 0

  # Truncate high weights
  trunc_value <- trunc_factor * tol
  bin$binw[bin$binw > trunc_value] <- trunc_value

  # Marching Cubes density level
  level_bin <- level_factor * trunc_value

  # Binning coordinates
  coorval <- npsp::coordvalues(bin)

  # Marching Cubes
  utils::capture.output({

    contours <- with(
      coorval,
      misc3d::contour3d(
        bin$binw,
        level = level_bin,
        X, Y, Z,
        draw = FALSE,
        separate = TRUE
      )
    )

  })

  if (length(contours) == 0)
    stop("Marching Cubes did not generate any mesh")

  # Bounding box of a Triangles3D object
  triangleBBox <- function(x) {

    vertices <- rbind(
      x$v1,
      x$v2,
      x$v3
    )

    apply(vertices, 2, range)
  }

  # Order disconnected components from largest to smallest
  contours <- contours[
    order(
      sapply(contours, function(x) nrow(x$v1)),
      decreasing = TRUE
    )
  ]

  # Dimensions used as reference for mesh coverage
  if (is.null(mesh_cad)) {

    ref_dim <- dimlen

  } else {

    if (!inherits(mesh_cad, "mesh3d"))
      stop("Argument 'mesh_cad' must be an object of class 'mesh3d'")

    cad_vertices <- t(
      mesh_cad$vb[1:3, , drop = FALSE]
    )

    cad_bbox <- apply(
      cad_vertices,
      2,
      range
    )

    ref_dim <- cad_bbox[2, ] - cad_bbox[1, ]
  }

  # Merge disconnected components until the required coverage is reached
  surface <- NULL

  for (current_surface in contours) {

    if (is.null(surface)) {

      surface <- current_surface

    } else {

      surface$v1 <- rbind(
        surface$v1,
        current_surface$v1
      )

      surface$v2 <- rbind(
        surface$v2,
        current_surface$v2
      )

      surface$v3 <- rbind(
        surface$v3,
        current_surface$v3
      )
    }

    current_bbox <- triangleBBox(surface)
    current_dim <- current_bbox[2, ] - current_bbox[1, ]

    if (all(current_dim >= ref_dim * coverage_tol))
      break
  }

  # Convert Triangles3D to mesh3d
  t2ve <- utils::getFromNamespace(
    "t2ve",
    "misc3d"
  )

  mesh0 <- t2ve(surface)

  mesh <- rgl::tmesh3d(
    vertices = mesh0[["vb"]],
    indices = mesh0[["ib"]]
  )

  return(
    list(
      mesh = mesh,
      bin = bin,
      contours = contours,
      params = list(
        resol_nbin = resol_nbin,
        nbin = nbin,
        zero_border = zero_border,
        k_factor = k_factor,
        tol = tol,
        trunc_factor = trunc_factor,
        trunc_value = trunc_value,
        level_factor = level_factor,
        level_bin = level_bin,
        coverage_tol = coverage_tol
      )
    )
  )
}
