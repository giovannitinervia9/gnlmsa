


#-------------------------------------------------------------------------------



#' Options for Simulated Annealing algorithm
#'
#' @param iterations integer indicating the number of iterations to perform.
#' @param compute_v how many iterations wait in order to recompute variance-covariance matrix.
#' @param initial_temperature initial temperature of the algorithm.
#' @param final_temperature final temperature of the algorithm.
#' @param restart_if_stuck how many iterations wait before restarting from last best parameters value if the algorithm stucks.
#' @param save_history whether to save full history of the algorithm.
#'
#' @return a list containing the specified options.
#' @export
sa_control <- function(iterations = 1000,
                       compute_v = floor(iterations/10),
                       initial_temperature = 100,
                       final_temperature = 1,
                       restart_if_stuck = floor(iterations*0.33),
                       save_history = FALSE
) {

  if (compute_v > iterations) {
    warning("compute_v must be less or equal to iterations, setting compute_v equal to iterations")
    compute_v <- iterations
  }

  if (initial_temperature < final_temperature) {
    stop("initial temperature must be greater than final temperature")
  }

  list(iterations = iterations,
       compute_v = compute_v,
       initial_temperature = initial_temperature,
       final_temperature = final_temperature,
       restart_if_stuck = restart_if_stuck,
       save_history = save_history)
}
