#' Transform age at time of screening data to intervals.
#' 
#' @noRd
#' @param data.obj A list object containing
#'   `$endpoint_type`, the censoring type (1, 2, 3, 4); 
#'    1: preclinical; 2: censored w/ screens; 3: clinical; 4: censored w/out screens
#'   `irateP`, the index of the rate_P parameters;
#'   `$n`, the number of cases;
#'   `$endpoint_time`, the censoring times;
#'   `$age_entry`, the ages at time of entry;
#'   `$ages_screen` a list, each element a vector of ages at time of screening; and 
#'   `$n_screen_positive`, an indicator of a positive screen.
#' @param t0 A scalar numeric object. The risk onset age.
#' 
#' @return A list object containing `$values`, a vector of 
#'   endpoints; `$starts` a vector of length n giving the first index of 
#'   `$values` pertaining to each case; `$ends`, a vector of length n giving 
#'   the last index of `$values` pertaining to each case; and `$lengths`, a 
#'   vector of length n giving the number of values of `$values` that pertain 
#'   to each individual. Note that the `$start` and `$stop` values are prepared 
#'   for use in C++, which indexes from 0.
#' 
#' @importFrom utils head
#' @keywords internal
.computeEndpoints <- function(data.obj, t0) {
  
  required_names <- c("endpoint_type", "n", "endpoint_time", "ages_screen", 
                      "n_screen_positive", "irateP", "ids")
  
  stopifnot(
    "`data.obj` does not contain required data" = !missing(data.obj) &&
      is.list(data.obj) && all(names(data.obj) %in% required_names),
    "`t0` must be a scalar numeric" = !missing(t0) && is.numeric(t0) &&
      is.vector(t0) && length(t0) == 1L
  )
  
  if (data.obj$endpoint_type == 1L) {
    screen <- lapply(seq_len(data.obj$n),
                     FUN = function(i, data, t0) {
                       idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                       if (any(data$ages_screen$values[idx] < t0)) {
                         stop("assessment ages younger than t0 encountered", 
                              call. = FALSE)
                       }
                       if (any(data$ages_screen$values[idx] > data$endpoint_time[i])) {
                         stop("assessments after endpoint encountered", call. = FALSE)
                       }
                       c(t0, data$ages_screen$values[idx]) |> unique() |> sort()
                     },
                     data = data.obj, t0 = t0)
    
    return(list("values"  = unlist(screen),
                "starts"  = utils::head(c(1L, cumsum(lengths(screen)) + 1L), -1L) - 1L,
                "ends"    = cumsum(lengths(screen)) - 1L,
                "lengths" = lengths(screen)))
  }
  
  if (data.obj$endpoint_type == 2L) {
    censored <- lapply(seq_len(data.obj$n),
                       FUN = function(i, data, t0) {
                         idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                         if (any(data$ages_screen$values[idx] < t0)) {
                           stop("assessment ages younger than t0 encountered", call. = FALSE)
                         }
                         if (any(data$ages_screen$values[idx] > data$endpoint_time[i])) {
                           stop("assessments after endpoint encountered", 
                                call. = FALSE)
                         }
                         c(t0, 
                           data$ages_screen$values[idx],
                           data$endpoint_time[i], Inf) |> unique()
                       },
                       data = data.obj, t0 = t0)
    
    return(list("values"  = unlist(censored),
                "starts"  = utils::head(c(1L, cumsum(lengths(censored)) + 1L), -1L) - 1L,
                "ends"    = cumsum(lengths(censored)) - 1L,
                "lengths" = lengths(censored)))
  }
           
  if (data.obj$endpoint_type == 3L) {
    if (length(data.obj$ages_screen) > 0L) {
      clinical <- lapply(seq_len(data.obj$n),
                         FUN = function(i, data, t0) {
                           idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                           if (any(data$ages_screen$values[idx] < t0)) {
                             stop("assessment ages younger than t0 encountered", 
                                  call. = FALSE)
                           }
                           if (any(data$ages_screen$values[idx] > data$endpoint_time[i])) {
                             stop("assessments after endpoint encountered", call. = FALSE)
                           }
                           c(t0, 
                             data$ages_screen$values[idx],
                             data$endpoint_time[i]) |> unique()
                         },
                         data = data.obj, t0 = t0)
    } else {
      clinical <- lapply(seq_len(data.obj$n),
                         FUN = function(i, data, t0) {
                           c(t0, 
                             data$endpoint_time[i]) |> unique()
                         },
                         data = data.obj, t0 = t0)
    }
    
    return(list("values"  = unlist(clinical),
                "starts"  = utils::head(c(1L, cumsum(lengths(clinical)) + 1L), -1L) - 1L,
                "ends"    = cumsum(lengths(clinical)) - 1L,
                "lengths" = lengths(clinical)))
  }
  
  if (data.obj$endpoint_type == 4L) {
    # censored without screens
    censored <- lapply(seq_len(data.obj$n),
                       FUN = function(i, data, t0) {
                         c(t0, data$endpoint_time[i], Inf) |> unique()
                       },
                       data = data.obj, t0 = t0)
  
    return(list("values"  = unlist(censored),
                "starts"  = utils::head(c(1L, cumsum(lengths(censored)) + 1L), -1L) - 1L,
                "ends"    = cumsum(lengths(censored)) - 1L,
                "lengths" = lengths(censored)))
  }
  

  NULL
}

# internal function to vectorize a jagged list
#
# `age.screen` is a list. Each element of the list provides the ages at
#   the time of screening for a single participant.
# The returned list contains: `$values`, a numeric vector, the input ages 
#   vectorized; `$starts`, an integer vector of length n, each element, i,
#   is the first index of `$values` pertaining to participant i; `$ends`,
#   an integer vector of length n, each element, i, is the last index of
#   `$values` pertaining to participant i; and `$lengths`, an integer
#   vector of length n, each element, i, is the number of elements of
#   `$value` that pertain to participant i. 
# NOTE: The indices provided in `$start` and `$stop` are generated for use 
#   in C++, which indexes vectors starting at 0. To use in R, must add
#   1L.
.defineAges <- function(age.screen, screen.types) {
  list("values"  = unlist(age.screen) |> unname(),
       "starts"  = unname(head(c(1L, cumsum(lengths(age.screen)) + 1L), -1L) - 1L),
       "ends"    = unname(cumsum(lengths(age.screen)) - 1L),
       "lengths" = unname(lengths(age.screen)),
       "types" = unname(unlist(screen.types)) - 1L)
}

.basicDataObject <- function(grp, endpoint.type,
                             data.clinical,
                             age.screen,
                             screen.types,
                             irateP, t0, idx) {
  
  tmp_list <- list("endpoint_type" = endpoint.type,
                   "irateP" = irateP,
                   "n" = sum(grp),
                   "ids" = data.clinical$id[grp],
                   "endpoint_time" = data.clinical$age_endpoint[grp],
                   "ages_screen" = if (!any(is.na(idx[grp]))) {
                                     .defineAges(age.screen[idx[grp]],
                                                 screen.types[idx[grp]])
                                   } else {
                                     list()
                                   },
                   "n_screen_positive" = rep(as.integer(endpoint.type == 1L), 
                                             sum(grp)))
  tmp_list$endpoints <- .computeEndpoints(tmp_list, t0)
  age_entry <- table(data.clinical$age_entry[grp]) |> as.data.frame()
  tmp_list$age_entry <- cbind(as.numeric(levels(age_entry[, 1L])[age_entry[, 1L]]),
                              age_entry[, 2L])
  
  tmp_list
  
}

#' @noRd
#' @param data.clinical A data.frame object. The clinical dataset.
#' @param age.screen A list object. Each element contains the ages at screenings
#'   for a single participant
#' @param screen.types A list object. Each element contains the types of 
#'   screens for a single participant
#' @param group.rateP An integer vector. The values indicate groups with common
#'   sojourn distributions.
#' @param idx An integer vector. The indices of age.screen and screen.types
#'   corresponding to the participants in data.clinical. NA indicates no
#'   screening history
#' @param t0 A scalar numeric. The risk onset age.
#' @keywords internal
.makeDataObjectsFull <- function(data.clinical, data.assess, groups.rateP, t0, all.types = TRUE) {
  
  screens  <- data.clinical$endpoint_type == "preclinical"
  censored <- data.clinical$endpoint_type == "censored"
  clinical <- data.clinical$endpoint_type == "clinical"
  
  if ({!any(screens) || !any(censored) || !any(clinical)} && all.types) {
    stop("all compartments must be represented in the data; ",
         "verify `data.clinical$endpoint_type", call. = FALSE)
  }
  
  # data.assess had been ordered by id and ages
  # extract screen ages and screen types for each participant as lists
  age_screen <- split(data.assess$age_assess, f = data.assess$id)
  screen_type <- split(data.assess$screen_type, f = data.assess$id)
  
  # there will be NAs if some participants do not have screens
  idx <- match(data.clinical$id, names(age_screen))
  
  data.obj <- list()
  
  study_arm_groups <- data.clinical$arm |> unique() |> sort()
  
  .message <- function(n, N, i, j, endpoint) {
    message(format(n, width = max(ceiling(log10(N)), 2)), 
            " participants in study arm ", study_arm_groups[i],
            ", rate_P group", names(groups.rateP)[j], " ", endpoint)
    
  }
  
  .stopMixed <- function(i, j) {
    stop("study arm ", study_arm_groups[i], 
         ", rate_P group", names(groups.rateP)[j], 
         " has missing screens for some participants;\n\t", 
         "mixed screening history within a group is not supported\n\t", 
         call. = FALSE)
  }
  
  cnt <- 1L
  
  for (i in 1L:length(study_arm_groups)) {
    for (j in 1L:length(groups.rateP)) {
      
      in_group <- data.clinical$arm == study_arm_groups[i] & 
                  data.clinical$grp.rateP == j
      
      if (!any(in_group)) {
        .message(0, nrow(data.clinical), i, j, "")
        next
      }
      
      if (any(in_group & screens)) {
        grp <- in_group & screens
        if (any(is.na(idx[grp]))) {
          stop("a participant marked as screen detected does not have screening data",
               call. = FALSE)
        }
        
        data.obj[[cnt]] <- .basicDataObject(grp = grp, 
                                            endpoint.type = 1L,
                                            data.clinical = data.clinical,
                                            age.screen = age_screen,
                                            screen.types = screen_type,
                                            irateP = j - 1L, 
                                            t0 = t0,
                                            idx = idx)
        
        .message(data.obj[[cnt]]$n, nrow(data.clinical), i, j, "are preclinical")
        cnt <- cnt + 1L
      } else {
        .message(0, nrow(data.clinical), i, j, "are preclinical")
      }
      
      if (any(in_group & censored)) {
        grp <- in_group & censored
        # all participants have a screen or all participants do not have a screen
        if (!{all(is.na(idx[grp])) || all(!is.na(idx[grp]))}) .stopMixed(i, j)

        if (all(!is.na(idx[grp]))) {
          # all have screening history
          data.obj[[cnt]] <- .basicDataObject(grp = grp, 
                                              endpoint.type = 2L,
                                              data.clinical = data.clinical,
                                              age.screen = age_screen,
                                              screen.types = screen_type,
                                              irateP = j - 1L, 
                                              t0 = t0,
                                              idx = idx)
          .message(data.obj[[cnt]]$n, nrow(data.clinical), i, j, 
                   "are censored with screening data")
        } else {
          # accommodate control group with no screens
          data.obj[[cnt]] <- .basicDataObject(grp = grp, 
                                              endpoint.type = 4L,
                                              data.clinical = data.clinical,
                                              age.screen = list(),
                                              screen.types = list(),
                                              irateP = j - 1L, 
                                              t0 = t0,
                                              idx = idx)
          .message(data.obj[[cnt]]$n, nrow(data.clinical), i, j, 
                   "are censored with no screening data")
        }
        
        cnt <- cnt + 1L
      } else {
        .message(0, nrow(data.clinical), i, j, "are censored")
      }
      
      if (any(in_group & clinical)) {
        grp <- in_group & clinical
        # all participants have a screen or all participants do not have a screen
        if (!{all(is.na(idx[grp])) || all(!is.na(idx[grp]))}) .stopMixed(i, j)
        
        data.obj[[cnt]] <- .basicDataObject(grp = grp, 
                                            endpoint.type = 3L,
                                            data.clinical = data.clinical,
                                            age.screen = age_screen,
                                            screen.types = screen_type,
                                            irateP = j - 1L, 
                                            t0 = t0,
                                            idx = idx)
        
        .message(data.obj[[cnt]]$n, nrow(data.clinical), i, j, "are clinical")
        cnt <- cnt + 1L
      } else {
        .message(0, nrow(data.clinical), i, j, "are clinical")
      }
    }
  }
  
  data.obj
}
