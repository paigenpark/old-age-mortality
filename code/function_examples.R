age <- seq(80, 110, by = 1)
alpha <- 0.0003
beta <- 0.09

# gompertz
gompertz_fun <- function(age, alpha, beta) {
  alpha*exp(beta*age)
}

gom <- gompertz_fun(age, alpha, beta)

  pdf("gompertz.pdf", width = 3, height = 6)

plot(age, gom, type = "l", col = "blue", lwd = 2, 
     xlab = "Age", ylab = "",
     ylim = c(0, 1.5))

  dev.off()
  system("open gompertz.pdf")

kannisto_fun <- function(age, alpha, beta) {
  (alpha*exp(beta*age)) / (1 + alpha*exp(beta*age))
}

# kannisto
kan <- kannisto_fun(age, alpha, beta)

  pdf("kannisto.pdf", width = 3, height = 6)

plot(age, kan, type = "l", col = "blue", lwd = 2, 
     xlab = "Age", ylab = "",
     ylim = c(0, max(kan) * 1.1))

  dev.off()
  system("open kannisto.pdf")
  
# beard
  
beard_fun <- function(age, alpha, beta, delta) {
    (alpha*exp(beta*age)) / (1 + delta*exp(beta*age))
  }

delta <- 0.00055

  beard <- beard_fun(age, alpha, beta, delta)
  
  pdf("beard.pdf", width = 3, height = 6)
  
  plot(age, beard, type = "l", col = "blue", lwd = 2, 
       xlab = "Age", ylab = "",
       ylim = c(0, max(beard) * 1.1))
  
  dev.off()
  system("open beard.pdf")