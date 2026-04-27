###############################################################
# Utility Function: Filter files by numeric suffix
###############################################################

filter_paths <- function(paths, keep, digits = 6) {
  ids <- sub(".*(?:\\.|_)(\\d+)\\.RData$", "\\1", paths, perl = TRUE)
  ids[grepl("^[0-9]+$", ids) == FALSE] <- NA
  ids_num <- as.numeric(ids)
  mod <- 10^digits
  last_digits <- ids_num %% mod
  idx <- !is.na(last_digits) & (last_digits %in% keep)
  paths[idx]
}

mysum <- function(x){y <- sum(x, na.rm = TRUE); return(y)}

fix_LC_arrays <- function(output_netCDFs, verbose = TRUE) {
  stopifnot(is.list(output_netCDFs))
  
  for (exp_name in names(output_netCDFs)) {
    if (verbose) message("Fixing: ", exp_name)
    arr <- output_netCDFs[[exp_name]]
    
    if (!is.array(arr) || length(dim(arr)) < 4) {
      stop("Element '", exp_name, "' is not a 4D array.")
    }
    
    dn <- dimnames(arr)
    dims <- dim(arr)
    nx <- dims[1]; ny <- dims[2]; nclass <- dims[3]; nyears <- dims[4]
    
    # Work on a copy
    A <- arr
    
    # 1) Clip everything to [0,1] (keep NAs)
    # We must be careful to not turn NA into 0:
    not_na_idx <- !is.na(A)
    A[not_na_idx] <- pmin(1, pmax(0, A[not_na_idx]))
    
    # 2) For each year: compute pixel sums (respect NA-only)
    for (t in seq_len(nyears)) {
      slice <- A[,,,t, drop = FALSE] # keep dims
      
      # detect pixels where ALL class values are NA
      all_na_mat <- apply(slice, c(1,2), function(x) all(is.na(x)))
      
      # compute pixel sums, returning NA where all values are NA
      pix_sum <- apply(slice, c(1,2), function(x) {
        if (all(is.na(x))) return(NA_real_)
        sum(x, na.rm = TRUE)
      })
      
      # if any pixel sum > 1, scale that pixel so the sum becomes 1
      exceed_mask <- (!is.na(pix_sum)) & (pix_sum > 1)
      if (any(exceed_mask)) {
        # create scale factor matrix f: 1/pix_sum where exceed, else 1
        f <- matrix(1, nrow = nx, ncol = ny)
        f[exceed_mask] <- 1 / pix_sum[exceed_mask]
        
        # vectorized multiplication across classes: multiply each class layer by f
        for (ccc in seq_len(nclass)) {
          layer <- A[,,ccc,t]
          # only multiply non-NA entries (NA remain NA)
          layer[!is.na(layer)] <- layer[!is.na(layer)] * f[!is.na(layer)]
          A[,,ccc,t] <- layer
        }
        
        # recompute pix_sum for diagnostics if needed (not strictly necessary)
        pix_sum <- apply(A[,,,t, drop = FALSE], c(1,2), function(x) {
          if (all(is.na(x))) return(NA_real_)
          sum(x, na.rm = TRUE)
        })
      }
      
      # safety: enforce that numeric pixel sums are <= 1 + tiny tol
      tol <- 1e-10
      if (any(pix_sum > 1 + tol, na.rm = TRUE)) {
        warning(sprintf("After scaling, some pixel sums still > 1 in %s year index %s (max = %g).",
                        exp_name, t, max(pix_sum, na.rm = TRUE)))
        # final hard clip per-pixel: scale any remaining by their sum
        exceed_mask2 <- (!is.na(pix_sum)) & (pix_sum > 1 + tol)
        if (any(exceed_mask2)) {
          f2 <- matrix(1, nrow = nx, ncol = ny)
          f2[exceed_mask2] <- 1 / pix_sum[exceed_mask2]
          for (ccc in seq_len(nclass)) {
            layer <- A[,,ccc,t]
            layer[!is.na(layer)] <- layer[!is.na(layer)] * f2[!is.na(layer)]
            A[,,ccc,t] <- layer
          }
        }
      }
      
      # ensure NA-only pixels remain NA (they should already)
      for (ccc in seq_len(nclass)) {
        layer <- A[,,ccc,t]
        layer[all_na_mat] <- NA_real_
        A[,,ccc,t] <- layer
      }
    } # end year loop
    
    # restore dimnames
    dimnames(A) <- dn
    output_netCDFs[[exp_name]] <- A
  } # end experiments loop
  
  invisible(output_netCDFs)
}



test_max_share_sums <- function(output_netCDFs, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  
  results <- data.frame(
    experiment = character(),
    year = character(),
    quantile = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (exp_name in names(output_netCDFs)) {
    
    arr <- output_netCDFs[[exp_name]]
    years <- dimnames(arr)[[4]]
    nyears <- dim(arr)[4]
    
    for (t in seq_len(nyears)) {
      
      # sum over land-cover classes (dim 3), keep NA pixels handled
      pix_sum <- apply(arr[,,,t, drop = FALSE], c(1,2), function(x) {
        if (all(is.na(x))) return(NA_real_)
        sum(x, na.rm = TRUE)
      })
      
      # compute quantiles, ignoring NA
      q_vals <- quantile(pix_sum, probs = probs, na.rm = TRUE, names = FALSE)
      
      # store in results
      for (i in seq_along(probs)) {
        results <- rbind(
          results,
          data.frame(
            experiment = exp_name,
            year = years[t],
            quantile = paste0(probs[i]*100, "%"),
            value = q_vals[i],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  return(results)
}
