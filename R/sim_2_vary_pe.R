# determine pe(230 = 1) and alpha based on pe(230 = 0)

pe_230_0 <- c(0, .1, .2, .3)

get_pe_230_1 <- function(pe_230_0, gamma, pe_overall = 0.7) {
  numerator <- log(1 - pe_overall) - (1 - gamma) * log(1 - pe_230_0)
  denominator <- gamma
  return(1 - exp(numerator / denominator))
}

pe_1s <- get_pe_230_1(pe_230_0, gamma = gammas$prop_placebo[gammas$pos == 230], pe_overall = 0.7)
round(pe_1s, 2)
round(log(1 - round(pe_1s, 2)) - log(1 - round(pe_230_0, 2)), 2)
