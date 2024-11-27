crowding_index_trajectory <- function(movement_models, times, scale=0.) {
  # Calculate crowding index for each time point
  crowding_index <- sapply(times, function(t) {
    # Get the mean positions and covariance matrices for each time point
    position_mean <- lapply(movement_models, function(model) {
      return(as.numeric(as.telemetry(model, t)))
    })
    cov_matrix <- lapply(movement_models, function(model) {
      return(as.matrix(as.telemetry(model, t, "covariance")))
    })
    # Calculate the crowding index
    return(calculate_crowding_index(position_mean, cov_matrix, scale))
  })
  return(crowding_index)
}


#' Calculate Crowding Index
#'
#' This function calculates the crowding index based on the provided mean positions and covariance matrices.
#'
#' @param position_mean A list of mean positions.
#' @param cov_matrix A list of covariance matrices corresponding to the mean positions.
#' @param scale A scaling factor for the covariance matrices. Default is 1.
#'
#' @return The average crowding index.
#'
#' @details
#' The function first checks if the dimensions of the mean positions and covariance matrices match.
#' It then ensures that there are at least two mean positions. The function calculates the differences
#' between the mean positions and the summed covariance matrices, and then computes the pairwise crowding
#' indices using the `gauss_local_crowd` function. Any NA values (resulting from ill defined covariance)
#' are replaced with zero. The final result is the average of the pairwise crowding indices.
#'
#' @examples
#' position_mean <- list(c(1, 2), c(3, 4))
#' cov_matrix <- list(matrix(c(1, 0, 0, 1), nrow = 2), matrix(c(2, 0, 0, 2), nrow = 2))
#' calculate_crowding_index(position_mean, cov_matrix)
#'
#' @export
calculate_crowding_index <- function(position_mean, cov_matrix, scale = 1) {
  # Check if dimensions match
  if (length(position_mean) != length(cov_matrix)) {
    print("Number of means:")
    stop("The number of means and covariance matrices must match")
  }
  n <- length(position_mean)

  # Check if there are at least 2 means
  if (n < 2) {
    stop("At least two means are required")
  }

  # calculate distances between means and summed covariance matrices
  summed_cov <- list()
  delta_mean <- list()
  k <- 1
  scale_matrix <- scale * scale * diag(2)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i <= j) {
        next
      }
      summed_cov[[k]] <- cov_matrix[[i]] + cov_matrix[[j]] + scale_matrix
      delta_mean[[k]] <- position_mean[[i]] - position_mean[[j]]
      k <- k + 1
    }
  }

  n_pairs <- k - 1

  # Calculate pairwise crowding
  crowd_idxs <- sapply(1:n_pairs, function(x) {
    return(gauss_local_crowd(delta_mean[[x]], summed_cov[[x]]))
  })
  # change NA values (det=0) to 0
  crowd_idxs[is.na(crowd_idxs)] <- 0.

  # Return the average crowding index
  return(mean(crowd_idxs))
}

#' Calculate Gaussian Local Crowding Index
#'
#' This function calculates the crowding index between two organisms
#' given the mean displacement and covariance matrix using a Gaussian kernel.
#'
#' @param position_mean A numeric vector representing the mean displacement.
#' @param cov_matrix A numeric matrix representing the covariance in the position.
#'
#' @warning The function will issue a warning if the determinant of the covariance matrix is not positive.
#' In this case, the function will return NA.
#' 
#' @return A numeric value representing the crowding index. If the determinant of the covariance
#' matrix is not positive, a warning is issued and the function returns NA.
#'
#' @examples
#' position_mean <- c(1, 2)
#' cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' gauss_local_crowd(position_mean, cov_matrix)
#'
#' @export
gauss_local_crowd <- function(position_mean, cov_matrix) {
  mat_det <- det(cov_matrix)
  if (mat_det <= 0) {
    warning("Determinant of the covariance matrix is not positive, please check the matrix")
    return(NA)
  }
  # crowd rate for a gaussian kernel
  crowding_idx <- 1 / (2 * pi * sqrt(mat_det)) * exp(-0.5 * t(position_mean) %*% solve(cov_matrix) %*% position_mean)
  return(crowding_idx)
}


