# Exact symmetric binary Mk simulation with one rate per whole branch.
# Q on branch b is matrix(c(-q[b], q[b], q[b], -q[b]), 2, 2).
# q follows the input tree's edge order. Root labels are "0" and "1".
# Unlike simulate_mk_epochs(), this model does not switch within a branch.
#
# For symmetric binary Q the exit rate is state-independent. Thus the number
# of flips is Poisson(q * duration), with sorted uniform event times given
# that number. Descendants inherit their parent's ending state exactly.
# Returns tip/node states, actual branch counts, and optional true simmaps.

simulate_binary_branches <- function(tree, q, root_state, Nchar,
                                     keep_history = TRUE) {
  if (!inherits(tree, "phylo")) stop("tree must be an ape phylo.")
  if (length(Nchar) != 1 || !is.finite(Nchar) || Nchar < 1 || Nchar != floor(Nchar)) {
    stop("Nchar must be a positive integer.")
  }
  n_edge <- nrow(tree$edge)
  n_tip <- length(tree$tip.label)
  if (!is.numeric(q) || length(q) != n_edge || any(!is.finite(q)) || any(q < 0)) {
    stop("Supply one finite non-negative q per input branch.")
  }
  if (length(tree$edge.length) != n_edge || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)) stop("Invalid branch lengths.")
  if (length(keep_history) != 1 || is.na(keep_history) || !is.logical(keep_history)) {
    stop("keep_history must be TRUE or FALSE.")
  }
  root_state <- as.character(root_state)
  if (length(root_state) == 1) root_state <- rep(root_state, Nchar)
  if (length(root_state) != Nchar || anyNA(root_state) || any(!root_state %in% c("0", "1"))) {
    stop("Supply one root state, or one 0/1 root state per character.")
  }
  original_child <- tree$edge[, 2]
  tree <- ape::reorder.phylo(tree, "cladewise")
  q <- q[match(tree$edge[, 2], original_child)]
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])
  if (length(root) != 1 || root != n_tip + 1) stop("Invalid directed root numbering.")
  character_names <- paste0("CH", 1:Nchar)
  expected_counts <- q * tree$edge.length
  counts <- matrix(rpois(n_edge * Nchar, lambda = expected_counts), n_edge, Nchar,
                    dimnames = list(paste0("BR", 1:n_edge), character_names))
  node_states <- matrix(NA_character_, n_tip + tree$Nnode, Nchar,
                        dimnames = list(as.character(1:(n_tip + tree$Nnode)), character_names))
  node_states[root, ] <- root_state
  maps <- NULL
  if (keep_history) maps <- vector("list", Nchar)

  for (character_id in 1:Nchar) {
    if (keep_history) {
      branch_maps <- vector("list", n_edge)
      mapped_edge <- matrix(0, n_edge, 2, dimnames = list(NULL, c("0", "1")))
    }
    for (branch_id in 1:n_edge) {
      parent <- tree$edge[branch_id, 1]
      child <- tree$edge[branch_id, 2]
      start <- as.integer(node_states[parent, character_id])
      if (is.na(start)) stop("Parent state is unavailable during traversal.")
      n_changes <- counts[branch_id, character_id]
      node_states[child, character_id] <- as.character((start + n_changes) %% 2)
      if (keep_history) {
        events <- sort(runif(n_changes, min = 0, max = tree$edge.length[branch_id]))
        dwell <- diff(c(0, events, tree$edge.length[branch_id]))
        states <- as.character((start + 0:n_changes) %% 2)
        branch_maps[[branch_id]] <- setNames(dwell, states)
        mapped_edge[branch_id, 1] <- sum(dwell[states == "0"])
        mapped_edge[branch_id, 2] <- sum(dwell[states == "1"])
      }
    }
    if (keep_history) {
      mapped <- tree
      mapped$maps <- branch_maps
      mapped$mapped.edge <- mapped_edge
      mapped$states <- setNames(node_states[1:n_tip, character_id], tree$tip.label)
      mapped$node.states <- cbind(node_states[tree$edge[, 1], character_id],
                                  node_states[tree$edge[, 2], character_id])
      class(mapped) <- c("simmap", "phylo")
      maps[[character_id]] <- mapped
    }
  }
  tip_states <- node_states[1:n_tip, , drop = FALSE]
  rownames(tip_states) <- tree$tip.label
  if (keep_history) {
    names(maps) <- character_names
    class(maps) <- c("multiSimmap", "multiPhylo")
  }
  list(tree = tree, q = q, tip_states = tip_states, node_states = node_states,
       changes_by_branch = counts, maps = maps)
}
