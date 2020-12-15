
# Mutiple alignment toolkit

library(doParallel)
library(foreach)

f.algin.score <- function(
  seqs, go="default", ge="default", 
  subm="default", gapval=0
){ 
  # Compute score
  sum(
    # Compute conservation matrix
    msa::msaConservationScore(
      # compute alignment using ClustalW
      msa::msaClustalW(
        seqs, substitutionMatrix = subm, 
        gapOpening = go, gapExtension = ge,
        verbose = FALSE
      ), 
      substitutionMatrix = subm, gapVsGap = gapval
    )
  )
}

linear.search <- function(
  seqs, subm="default", 
  gapval=0, # How to penalise for gaps when computing the score.
  go = seq(5, 15, 1), # interval for gap Opening penalties
  ge = seq(0.1, 1.5, 0.05) # ibidem for gap extension
){
  retval <- foreach(g_op = go, .combine = rbind) %do% {
    foreach(g_ext = ge, .combine = rbind) %do% {
      it.algn.score <- f.algin.score(seqs, g_op, g_ext, subm, gapval)
      data.frame(gop = g_op, gext = g_ext, score = it.algn.score)
    }
  }
  return(retval)
}

parallel.search <- function(
  seqs, subm, gapval, 
  go = seq(5, 15, 1), ge = seq(0.1, 1.5, 0.05), 
  nthreads=(detectCores()/2)
){
  cl <- parallel::makeCluster(nthreads, type = "FORK")
  doParallel::registerDoParallel(cl)
  
  #nreps <- length(ge) * length(go)
  #i.seqs <- iter(replicate(nreps, seqs, simplify = FALSE))
  grid.eval <- foreach(g_op = go, .combine = rbind) %dopar% {
    foreach(g_ext = ge, .combine = rbind) %dopar% {
      it.algn.score <- f.algin.score(
        seqs, g_op, g_ext, subm, gapval
      )
      #data.frame(gop = g_op, gext = g_ext)
      data.frame(gop = g_op, gext = g_ext, score = it.algn.score)
    }
  }
  stopCluster(cl) 
  return(grid.eval)
}


lin.max.score <- function(x){ x[x$score == max(x$score), ] }
lin.min.score <- function(x){ x[x$score == min(x$score), ] }

plot.align.search <- function(the.table, ttl="GridSearch"){
  gap.Open <- factor(the.table$gop)
  gap.Ext <- factor(the.table$gext)
  algn.score <- factor(the.table$score)
  ggplot2::qplot(
    x = gap.Open, y = gap.Ext, 
    fill = algn.score, geom = 'tile', main = ttl
  )
}
