nsim <- 2500

auc <- vector("numeric", nsim)
accuracy <- vector("numeric", nsim)
set.seed(20230329)
for (i in 1:nsim) {
  b <- rbinom(1e5, size = 1, prob = 0.05)
  majority_probs <- rep(1 - mean(b), length(b))
  majority_classifier <- rep(0, length(b))
  accuracy[i] <- mean(majority_classifier == b)
  auc[i] <- cvAUC::cvAUC(predictions = majority_probs, labels = b)$fold.AUC
}
summary(accuracy)
summary(auc)
