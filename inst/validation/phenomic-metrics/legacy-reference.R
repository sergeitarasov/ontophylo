# Read-only reference formulas from Hymenoptera_paper/R/utils-ST.R.
# Kept only for numerical compatibility validation, not loaded by ontophylo.
# Source MD5: a06db555bae64f8024ef711c98208768
get_enrichment_mask <- function(rate_matrix, lower_q = 0.025, upper_q = 0.975, direction = c("both", 
    "over", "under"), return_type = c("binary", "masked")) {
    direction <- match.arg(direction)
    return_type <- match.arg(return_type)
    enrichment_mask <- matrix(0, nrow = nrow(rate_matrix), ncol = ncol(rate_matrix))
    colnames(enrichment_mask) <- colnames(rate_matrix)
    rownames(enrichment_mask) <- rownames(rate_matrix)
    for (col in colnames(rate_matrix)) {
        col_vals <- rate_matrix[, col]
        q <- quantile(col_vals, probs = c(lower_q, upper_q), na.rm = TRUE)
        if (direction == "both") {
            mask <- col_vals < q[1] | col_vals > q[2]
        }
        else if (direction == "over") {
            mask <- col_vals >= q[2]
        }
        else {
            mask <- col_vals <= q[1]
        }
        if (return_type == "binary") {
            enrichment_mask[, col] <- as.numeric(mask)
        }
        else {
            enrichment_mask[, col] <- ifelse(mask, col_vals, 0)
        }
    }
    return(enrichment_mask)
} 

diversity_metric <- function(states, type = c("Neff", "Neff_scaled", "weighted"), base = exp(1)) {
    type <- match.arg(type)
    shannon_entropy <- function(x) {
        freqs <- table(x)/length(x)
        -sum(freqs * log(freqs, base = base))
    }
    if (type == "Neff") {
        enriched_states <- states[grepl("1", states)]
        if (length(enriched_states) == 0) 
            return(0)
        H_plus <- shannon_entropy(enriched_states)
        return(base^H_plus)
    }
    if (type == "Neff_scaled") {
        enriched_states <- states[grepl("1", states)]
        prop_enriched <- mean(grepl("1", states))
        if (length(enriched_states) == 0) 
            return(0)
        H_plus <- shannon_entropy(enriched_states)
        Neff_plus <- base^H_plus
        return(Neff_plus * prop_enriched)
    }
    if (type == "weighted") {
        H <- shannon_entropy(states)
        f_enriched <- mean(grepl("1", states))
        return(H + f_enriched)
    }
} 

relative_null_metric <- function(states, type = c("not_scaled", "scaled")) {
    type <- match.arg(type)
    if (length(states) == 0) 
        return(NA)
    len <- nchar(states[1])
    enriched <- states[grepl("1", states)]
    n <- length(enriched)
    if (n == 0) 
        return(0)
    mat <- strsplit(enriched, split = "")
    mat <- do.call(rbind, mat)
    mat <- apply(mat, 2, as.integer)
    if (is.null(dim(mat))) {
        mat <- matrix(mat, nrow = 1)
    }
    d_null <- rowSums(mat)
    rel_dist <- mean(d_null/len)
    if (type == "not_scaled") {
        return(rel_dist)
    }
    if (type == "scaled") {
        prop_enriched <- mean(grepl("1", states))
        return(rel_dist * prop_enriched)
    }
} 

