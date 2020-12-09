
library(msa) # Paquet pour alignements multiples
library(aphid) # Paquet contenant les matrices de substitution
library(foreach) # Paquet pour itérations
data("substitution") # Substitution matrices
source("mAlignToolkit.R")

w.search.bloss60 <- linear.search(
  taxo.fasta, subm=substitution$BLOSUM60, 
  gapval=0, # How to penalise for gaps when computing the score.
  go = seq(5, 15, 1), # interval for gap Opening penalties
  ge = seq(0.1, 1.5, 0.1) # ibidem for gap extension
)
