PepIsotopeDistribution <- function (sequence, IAA = TRUE, charge = 1, series = NULL) {
  # taken from package OrgMassSpec function IsotopeDistrubtionN and modified by Rico Derks
  if(sequence == "") {
    return(data.frame(mz = 0,
                      percent = 0))
  }
  
  if (charge < 1 | charge > 3) 
    stop("charge must be 1, 2, or 3")
  seq_vector <- strsplit(sequence, split = "")[[1]]
  x <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)
  for (i in 1:(length(seq_vector))) {
    if (seq_vector[i] == "A") 
      x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "R") 
      x <- x + c(C = 6, H = 12, N = 4, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "N") 
      x <- x + c(C = 4, H = 6, N = 2, O = 2, S = 0, P = 0)
    if (seq_vector[i] == "D") 
      x <- x + c(C = 4, H = 5, N = 1, O = 3, S = 0, P = 0)
    if (seq_vector[i] == "E") 
      x <- x + c(C = 5, H = 7, N = 1, O = 3, S = 0, P = 0)
    if (seq_vector[i] == "Q") 
      x <- x + c(C = 5, H = 8, N = 2, O = 2, S = 0, P = 0)
    if (seq_vector[i] == "G") 
      x <- x + c(C = 2, H = 3, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "H") 
      x <- x + c(C = 6, H = 7, N = 3, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "I") 
      x <- x + c(C = 6, H = 11, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "L") 
      x <- x + c(C = 6, H = 11, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "K") 
      x <- x + c(C = 6, H = 12, N = 2, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "M") 
      x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 1, P = 0)
    if (seq_vector[i] == "F") 
      x <- x + c(C = 9, H = 9, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "P") 
      x <- x + c(C = 5, H = 7, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "S") 
      x <- x + c(C = 3, H = 5, N = 1, O = 2, S = 0, P = 0)
    if (seq_vector[i] == "T") 
      x <- x + c(C = 4, H = 7, N = 1, O = 2, S = 0, P = 0)
    if (seq_vector[i] == "W") 
      x <- x + c(C = 11, H = 10, N = 2, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "Y") 
      x <- x + c(C = 9, H = 9, N = 1, O = 2, S = 0, P = 0)
    if (seq_vector[i] == "V") 
      x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 0, P = 0)
    if (seq_vector[i] == "C" & IAA == TRUE) 
      x <- x + c(C = 5, H = 8, N = 1, O = 2, S = 1, P = 0)
    if (seq_vector[i] == "C" & IAA == FALSE) 
      x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 1, P = 0)
  }
  elements <- x + c(C = 0, H = 2, N = 0, O = 1, S = 0, P = 0)
  simulation <- function(elements) {
    mz <- vector(mode = "double")
    mc <- sum(sample(c(12, 13.0033548378), size = elements["C"], 
                     replace = TRUE, prob = c(0.9893, 0.0107)))
    mh <- sum(sample(c(1.0078250321, 2.014101778), size = elements["H"], 
                     replace = TRUE, prob = c(0.999885, 0.000115)))
    mn <- sum(sample(c(14.0030740052, 15.0001088984), size = elements["N"], 
                     replace = TRUE, prob = c(0.99636, 0.00364)))
    if (IAA == TRUE) {
      mn_iaa <- sum(sample(c(14.0030740052, 15.0001088984), 
                           size = length(grep("C", seq_vector)), replace = TRUE, 
                           prob = c(0.99632, 0.00368)))
    }
    else mn_iaa <- 0
    mo <- sum(sample(c(15.9949146221, 16.9991315, 17.9991604), 
                     size = elements["O"], replace = TRUE, prob = c(0.99757, 
                                                                    0.00038, 0.00205)))
    ms <- sum(sample(c(31.97207069, 32.9714585, 33.96786683, 
                       35.96708088), size = elements["S"], replace = TRUE, 
                     prob = c(0.9493, 0.0076, 0.0429, 2e-04)))
    mch <- sum(sample(c(1.0072764522, 2.0135531981), size = charge, 
                      replace = TRUE, prob = c(0.999885, 0.000115)))
    mz <- sum(mc, mh, mn, mo, ms, mch, mn_iaa)/charge
    return(mz)
  }
  sim <- replicate(10000, expr = simulation(elements))
  b <- seq(from = min(sim) - (1/(2 * charge)), to = max(sim) + 
             1, by = 1/charge)
  bins <- cut(sim, breaks = b)
  mz <- tapply(sim, bins, mean)
  if(series == "b") {
    mz <- mz - 18.010565 / charge
  }
  intensity <- as.vector(table(bins))
  spec <- data.frame(mz, intensity)
  spec <- subset(spec, intensity != 0)
  spec <- transform(spec, percent = round(intensity/max(intensity) * 
                                            100, digits = 2))
  row.names(spec) <- 1:(nrow(spec))
  
  spec <- spec[1:4, ]
  return(spec)
}