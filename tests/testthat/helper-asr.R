asr_test_tree <- function() {
  ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
}

asr_test_data <- function(tree = asr_test_tree()) {
  data.frame(first = c("0", "0", "0", "0", "1", "1", "1", "1"),
             second = c("0", "0", "0", "1", "1", "1", "1", "1"),
             row.names = tree$tip.label)
}

# Independent enumeration of all internal-node states, not corHMM/phytools
# pruning. Tip observation likelihoods can include ambiguity or missing data.
asr_exact <- function(tree, tip, Q, root) {
  k <- nrow(Q)
  ntip <- length(tree$tip.label)
  node_states <- expand.grid(rep(list(seq_len(k)), tree$Nnode))
  transition <- vector("list", nrow(tree$edge))
  for (b in seq_along(transition)) {
    transition[[b]] <- expm::expm(Q * tree$edge.length[b])
  }
  weight <- numeric(nrow(node_states))
  for (i in seq_len(nrow(node_states))) {
    weight[i] <- 1
    for (b in seq_len(nrow(tree$edge))) {
      parent <- node_states[i, tree$edge[b, 1] - ntip]
      child <- tree$edge[b, 2]
      if (child <= ntip) {
        p <- sum(transition[[b]][parent, ] * tip[child, ])
      } else {
        p <- transition[[b]][parent, node_states[i, child - ntip]]
      }
      weight[i] <- weight[i] * p
    }
  }
  root_likelihood <- numeric(k)
  for (j in seq_len(k)) root_likelihood[j] <- sum(weight[node_states[, 1] == j])
  if (identical(root, "maddfitz")) root <- root_likelihood / sum(root_likelihood)
  if (identical(root, "equal")) root <- rep(1/k, k)
  if (identical(root, "yang")) {
    equations <- t(Q)
    equations[k, ] <- 1
    rhs <- c(rep(0, k - 1), 1)
    root <- solve(equations, rhs)
  }
  weight <- weight * root[node_states[, 1]]
  likelihood <- sum(weight)
  marginal <- matrix(0, tree$Nnode, k)
  for (node in seq_len(tree$Nnode)) {
    for (state in seq_len(k)) {
      marginal[node, state] <- sum(weight[node_states[, node] == state]) / likelihood
    }
  }
  list(logLik = log(likelihood), marginal = marginal)
}

# Known-Q fixture isolates conditional mapping from parameter estimation.
asr_known_fit <- function(tree, Q, tip, root = "equal") {
  index <- matrix(0L, nrow(Q), ncol(Q), dimnames = dimnames(Q))
  index[row(index) != col(index)] <- seq_len(nrow(Q) * (nrow(Q) - 1))
  character <- list(status = "ok", reason = "", Q_matrix = Q,
                    index_matrix = index, tips_prior = tip, root_p = root)
  structure(list(tree = tree, characters = list(known = character)),
             class = "ontophylo_fit")
}
