expected_lst_lefttrunc_d <- function(df, mu, sigma, tau) {
# Purpose: E[Y | Y > tau] for a Student-t with df, mu, sigma (vectorised; stable).
# Inputs: df (>1), mu, sigma, tau (broadcastable).
# Outputs: u = conditional mean(s) after truncation (same shape as mu).

  # Guard rails
  sigma <- pmax(sigma, 1e-12)  
  df    <- pmax(df, 1 + 1e-6)  # mean exists only for df > 1

  # Broadcast tau to match mu/sigma shape
  if (length(tau) == 1L) {
    tau <- tau + 0 * mu
  }

  taudot <- (tau - mu) / sigma # standardise

  logC <- lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(df * pi)

  # Numerator (in log space):
  log_num <- log(df / (df - 1)) + logC + ((1 - df) / 2) * log1p((taudot^2) / df)

  # Denominator: P(T_df > taudot)
  tail_prob <- pt(taudot, df = df, lower.tail = FALSE)
  tail_prob <- pmax(tail_prob, 1e-300)  # avoid divide-by-zero

  udot <- exp(log_num) / tail_prob # E[T | T > taudot]
  u    <- mu + sigma * udot             

  return(u)
}
