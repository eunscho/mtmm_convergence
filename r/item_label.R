item_label <- function(trait, method, item) {
  out <- vector("character", trait * method * item)
  for (i in 1:length(out)) {
    tnum <- mydiv(ceiling(i / item), trait)
    mnum <- ceiling(i / (item * trait))
    inum <- mydiv(i, item)
    out[i] <- paste0("t", tnum, "m", mnum, "_", inum)
  }
  return(out)
}