#' Prepare a Theoretical Panel Mesh for Point Cloud Simulation
#'
#' Prepares a theoretical steel panel mesh from a CAD `mesh3d` object for use
#' in point cloud simulation. The function removes the lower face of the panel,
#' identifies the base surface, determines a reference corner from its boundary,
#' translates the mesh to the origin and removes any remaining triangles below
#' the Z = 0 plane.
#'
#' The resulting mesh can be passed directly to `simulatePanelFloorCloud()` to
#' generate a simulated 3D point cloud representing a scanned steel panel.
#'
#' @param mesh A `mesh3d` object containing the CAD geometry of the steel panel.
#' The mesh must contain a triangular `it` matrix.
#' @param lower_angle Numeric. Angular threshold in degrees used to identify and
#' remove the lower face of the panel. Default is `120`.
#' @param horizontal_angle Numeric. Maximum angle in degrees used to identify
#' approximately horizontal faces. Default is `40`.
#' @param boundary_strip Numeric. Width of the strip, measured along the Y
#' direction, used to identify the reference corner of the panel base.
#' Default is `10`.
#'
#' @returns A `mesh3d` object containing the processed theoretical panel mesh,
#' translated to the origin and ready for point cloud simulation.
#'
#' @seealso `simulatePanelFloorCloud()`, `getBoundarySegments()`,
#' `angleFromAxis()`
#'
#' @examples
#' \dontrun{
#' data("cad", package = "dimControl")
#'
#' mesh_teor <- preparePanelMesh(cad)
#'
#' rgl::open3d()
#' rgl::shade3d(mesh_teor, col = "gray")
#' rgl::decorate3d()
#' }
#'
#' @export
preparePanelMesh <- function(mesh,
                             lower_angle = 120,
                             horizontal_angle = 40,
                             boundary_strip = 10) {

  if (!inherits(mesh, "mesh3d"))
    stop("Argument 'mesh' must be an object of class 'mesh3d'")

  if (is.null(mesh$it))
    stop("Argument 'mesh' must contain a triangular 'it' matrix")

  # ---------------------------------------------------------------------------
  # Remove lower face of the panel
  # ---------------------------------------------------------------------------

  normals <- Rvcg::vcgFaceNormals(mesh)

  angle_z <- angleFromAxis(
    normals,
    dim = 3,
    deg = TRUE
  )

  lower_face <- angle_z > lower_angle

  mesh_theoretical <- mesh
  mesh_theoretical$it <- mesh$it[, !lower_face, drop = FALSE]

  # Remove attributes that are no longer consistent with the modified mesh
  mesh_theoretical$normals <- NULL
  mesh_theoretical$tags <- NULL

  # Clean unused vertices
  mesh_theoretical <- cleanMesh3d(mesh_theoretical)

  # ---------------------------------------------------------------------------
  # Recalculate normals and barycenters
  # ---------------------------------------------------------------------------

  normals <- Rvcg::vcgFaceNormals(mesh_theoretical)

  angle_z <- angleFromAxis(
    normals,
    dim = 3,
    deg = TRUE
  )

  bary <- Rvcg::vcgBary(mesh_theoretical)

  # ---------------------------------------------------------------------------
  # Identify horizontal components and panel base
  # ---------------------------------------------------------------------------

  horizontal <- angle_z < horizontal_angle

  z_limit <- mean(bary[, 3])

  base <- horizontal & bary[, 3] < z_limit

  base_mesh <- mesh_theoretical
  base_mesh$it <- mesh_theoretical$it[, base, drop = FALSE]

  # Remove normals because they correspond to the complete mesh
  base_mesh$normals <- NULL

  # ---------------------------------------------------------------------------
  # Extract base boundary
  # ---------------------------------------------------------------------------

  boundary_mesh <- getBoundarySegments(
    base_mesh,
    malla = TRUE,
    simplify = TRUE
  )

  boundary_coordinates <- as.data.frame(
    t(boundary_mesh$vb[1:3, , drop = FALSE])
  )

  names(boundary_coordinates) <- c("X", "Y", "Z")

  # ---------------------------------------------------------------------------
  # Identify reference corner
  # ---------------------------------------------------------------------------

  lower_strip <- boundary_coordinates[
    boundary_coordinates$Y <=
      min(boundary_coordinates$Y) + boundary_strip,
    ,
    drop = FALSE
  ]

  corner <- lower_strip[
    which.min(lower_strip$X),
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------------------------
  # Translate theoretical mesh to the origin
  # ---------------------------------------------------------------------------

  mesh_theoretical <- rgl::translate3d(
    mesh_theoretical,
    -corner$X,
    -corner$Y,
    -corner$Z
  )

  # ---------------------------------------------------------------------------
  # Remove triangles containing vertices below Z = 0
  # ---------------------------------------------------------------------------

  idx <- which(mesh_theoretical$vb[3, ] >= 0)

  keep <- apply(
    mesh_theoretical$it,
    2,
    function(face) all(face %in% idx)
  )

  mesh_theoretical$it <- mesh_theoretical$it[
    ,
    keep,
    drop = FALSE
  ]

  # Remove attributes that are no longer valid
  mesh_theoretical$normals <- NULL
  mesh_theoretical$tags <- NULL

  # Clean unused vertices
  mesh_theoretical <- cleanMesh3d(mesh_theoretical)

  return(mesh_theoretical)
}


#' Simulate a Steel Panel and Floor Point Cloud
#'
#' Generates a simulated 3D point cloud representing a steel panel and the
#' surrounding floor, reproducing a typical acquisition scenario in which a
#' 3D scanner captures points belonging both to the inspected panel and to the
#' supporting environment.
#'
#' Panel points are sampled from a theoretical `mesh3d` surface using
#' `sampleMesh()` and perturbed with truncated normal noise. Floor points are
#' generated in four external regions surrounding the panel in the XY plane
#' and are placed below the panel with normally distributed vertical noise.
#' Both sets of points are combined into a single point cloud.
#'
#' The resulting data frame represents a simulated 3D scan containing points
#' from both the steel panel and the surrounding floor.
#'
#' @param mesh A theoretical `mesh3d` object representing the steel panel,
#' typically obtained with `preparePanelMesh()`.
#' @param n_panel Number of points sampled from the panel surface.
#' Default is `1e7`.
#' @param n_floor Number of points generated for the surrounding floor.
#' Default is `1e6`.
#' @param prec Positive numeric value defining the lower and upper bounds of
#' the truncated normal noise applied to panel coordinates. Noise is restricted
#' to the interval `[-prec, prec]`. Default is `0.5`.
#' @param sd_panel Numeric. Standard deviation of the truncated normal noise
#' applied independently to the X, Y and Z coordinates of the panel points.
#' Default is `0.2`.
#' @param sd_floor_z Numeric. Standard deviation of the vertical noise applied
#' to the floor points. Default is `0.2`.
#' @param margin_xy Numeric. Additional XY distance used to extend the simulated
#' floor beyond the panel bounding box. Default is `20`.
#' @param gap_z Numeric. Vertical distance between the lowest point of the panel
#' and the mean height of the simulated floor. Default is `2`.
#' @param seed Integer used to initialize the random number generator and make
#' the simulation reproducible. Default is `1`.
#'
#' @returns A data frame with three numeric columns, `X`, `Y` and `Z`,
#' containing the simulated panel and floor points. The first `n_panel` rows
#' correspond to panel points and the remaining `n_floor` rows correspond to
#' floor points.
#'
#' @seealso `preparePanelMesh()`, `sampleMesh()`
#'
#' @examples
#' \dontrun{
#' data("cad", package = "dimControl")
#'
#' mesh_teor <- preparePanelMesh(cad)
#'
#' data <- simulatePanelFloorCloud(
#'   mesh = mesh_teor,
#'   n_panel = 1e7,
#'   n_floor = 1e6,
#'   seed = 1
#' )
#'
#' dim(data)
#' head(data)
#'
#' rgl::open3d()
#' rgl::points3d(data, col = "gray")
#' rgl::decorate3d()
#' }
#'
#' @export
simulatePanelFloorCloud <- function(mesh,
                                    n_panel = 1e7,
                                    n_floor = 1e6,
                                    prec = 0.5,
                                    sd_panel = 0.2,
                                    sd_floor_z = 0.2,
                                    margin_xy = 20,
                                    gap_z = 2,
                                    seed = 1) {

  if (!inherits(mesh, "mesh3d"))
    stop("Argument 'mesh' must be an object of class 'mesh3d'")

  if (n_panel <= 0)
    stop("Argument 'n_panel' must be greater than 0")

  if (n_floor < 0)
    stop("Argument 'n_floor' must be greater than or equal to 0")

  if (prec <= 0)
    stop("Argument 'prec' must be greater than 0")

  if (sd_panel < 0 || sd_floor_z < 0)
    stop("Standard deviations must be greater than or equal to 0")

  # ---------------------------------------------------------------------------
  # Internal function: truncated normal random values
  # ---------------------------------------------------------------------------

  rtnormInternal <- function(n,
                             mean = 0,
                             sd = 1,
                             lower = -Inf,
                             upper = Inf) {

    p_lower <- stats::pnorm(
      lower,
      mean = mean,
      sd = sd
    )

    p_upper <- stats::pnorm(
      upper,
      mean = mean,
      sd = sd
    )

    stats::qnorm(
      stats::runif(
        n,
        min = p_lower,
        max = p_upper
      ),
      mean = mean,
      sd = sd
    )
  }

  set.seed(seed)

  n_panel <- as.integer(n_panel)
  n_floor <- as.integer(n_floor)

  # ---------------------------------------------------------------------------
  # Simulate panel points
  # ---------------------------------------------------------------------------

  panel <- sampleMesh(
    mesh,
    n_panel
  )

  panel <- t(panel)

  # Add truncated normal noise independently to X, Y and Z
  for (j in seq_len(3)) {

    panel[, j] <- panel[, j] +
      rtnormInternal(
        n_panel,
        mean = 0,
        sd = sd_panel,
        lower = -prec,
        upper = prec
      )
  }

  panel <- as.data.frame(panel)

  names(panel) <- c("X", "Y", "Z")

  # ---------------------------------------------------------------------------
  # Panel bounding box
  # ---------------------------------------------------------------------------

  xr <- range(panel$X)
  yr <- range(panel$Y)
  zr <- range(panel$Z)

  x_min_ext <- xr[1] - margin_xy
  x_max_ext <- xr[2] + margin_xy

  y_min_ext <- yr[1] - margin_xy
  y_max_ext <- yr[2] + margin_xy

  # ---------------------------------------------------------------------------
  # Number of floor points in each external region
  # ---------------------------------------------------------------------------

  n_left <- round(n_floor * 0.25)
  n_right <- round(n_floor * 0.25)
  n_bottom <- round(n_floor * 0.25)

  n_top <- n_floor -
    n_left -
    n_right -
    n_bottom

  # ---------------------------------------------------------------------------
  # Simulate floor points
  # ---------------------------------------------------------------------------

  floor_left <- data.frame(
    X = stats::runif(
      n_left,
      x_min_ext,
      xr[1]
    ),
    Y = stats::runif(
      n_left,
      y_min_ext,
      y_max_ext
    )
  )

  floor_right <- data.frame(
    X = stats::runif(
      n_right,
      xr[2],
      x_max_ext
    ),
    Y = stats::runif(
      n_right,
      y_min_ext,
      y_max_ext
    )
  )

  floor_bottom <- data.frame(
    X = stats::runif(
      n_bottom,
      xr[1],
      xr[2]
    ),
    Y = stats::runif(
      n_bottom,
      y_min_ext,
      yr[1]
    )
  )

  floor_top <- data.frame(
    X = stats::runif(
      n_top,
      xr[1],
      xr[2]
    ),
    Y = stats::runif(
      n_top,
      yr[2],
      y_max_ext
    )
  )

  floor <- rbind(
    floor_left,
    floor_right,
    floor_bottom,
    floor_top
  )

  # ---------------------------------------------------------------------------
  # Floor height
  # ---------------------------------------------------------------------------

  z_floor <- zr[1] - gap_z

  floor$Z <- stats::rnorm(
    nrow(floor),
    mean = z_floor,
    sd = sd_floor_z
  )

  # ---------------------------------------------------------------------------
  # Combine panel and floor points
  # ---------------------------------------------------------------------------

  data <- rbind(
    panel,
    floor
  )

  return(data)
}
