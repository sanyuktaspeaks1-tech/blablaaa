# Chapter 14
---

**CLT for Poisson sample mean**

Simulates M = 500 replications of the sample mean from Poisson(λ=3) for n ∈ {2, 3, 5, 10, 15, 25} and overlays the N(λ, λ/n) PDF on each histogram, demonstrating the CLT.

$$\bar{X}_n \sim \mathcal{N}\\left(\lambda,\, \frac{\lambda}{n}\right) \quad \text{for large } n$$

```r
lambda = 3 # true value of lambda 
par(mfrow = c(2,3))
n_vals = c(2,3,5,10,15,25) # varying sample size
M = 500 # number of replications
for(n in n_vals){ # loop on sample size
  sample_means = numeric(length = M)
  for(i in 1:M){ # loop on number of replications
    x = rpois(n = n, lambda = lambda)
    sample_means[i] = mean(x)
  }
  hist(sample_means, probability = TRUE, main = paste("n = ", n),
       xlab = expression(bar(X[n])), breaks = 20,
       cex.lab = 1.4, cex.main = 1.4)
  curve(dnorm(x, mean = lambda, sd = sqrt(lambda/n)),
        add = TRUE, col = "red", lwd = 2)
} 
```

---

**Delta method – sampling distribution of 1/X̄_n**

Simulates the sampling distribution of Y_n = 1/X̄_n for n ∈ {2, 3, 5, 10, 15, 25} and overlays the Delta method normal approximation, showing convergence as n grows.

$$Y_n = \frac{1}{\bar{X}_n} \sim \mathcal{N}\!\left(\frac{1}{\lambda},\, \frac{1}{\lambda^3 n}\right) \quad \text{for large } n$$

```r
lambda = 3 # true value of lambda
par(mfrow = c(2,3))
n_vals = c(2,3,5,10,15,25) # varying sample size
M = 500 # number of replications
for(n in n_vals){
  y_n = numeric(length = M)
  for(i in 1:M){
    x = rpois(n = n, lambda = lambda)
    y_n[i] = 1/mean(x)
  }
  hist(y_n, probability = TRUE, main = paste("n = ", n),
       xlab = expression(bar(Y[n])), breaks = 20,
       cex.lab = 1.4, cex.main = 1.4)
  curve(dnorm(x, mean = 1/lambda, sd = sqrt(1/(lambda^3*n))),
        add = TRUE, col = "red", lwd = 2)
  points(1/lambda,0, pch = 19, col = "blue", cex = 1.5)
}
```

---

**Delta method – sampling distribution of ψ(X̄_n)**

Simulates the sampling distribution of ψ(X̄_n) = 1 − (1 + X̄_n)e^(−X̄_n) for n ∈ {4, 10, 20, 50, 100, 500} and overlays the Delta method normal approximation, with the true ψ(λ) marked as a red dot.

$$\psi(\lambda) = P(X \geq 1) = 1 - (1+\lambda)e^{-\lambda}, \quad \text{Var}(\psi(\bar{X}_n)) \approx \frac{\lambda^3 e^{-2\lambda}}{n}$$

```r
par(mfrow = c(2,3))
n_vals = c(4, 10, 20, 50, 100, 500)
rep = 1000
lambda = 2
psi = 1-(1+lambda)*exp(-lambda)
for(n in n_vals){
  psi_vals = numeric(length = rep)
  for(i in 1:rep){
    x = rpois(n = n, lambda = lambda)
    psi_vals[i] = 1-(1+mean(x))*exp(-mean(x))
  }
  hist(psi_vals, probability = TRUE, cex.lab=1.5,
       xlab = expression(psi), main = paste("n = ", n))
  points(psi, 0, pch = 19, col = "red", cex = 2)
  curve(dnorm(x,mean = psi,sd=sqrt(lambda^3*exp(-2*lambda)/n)), add = TRUE, col = "red", lwd =2)
}
```

---

**CLT and Delta method – Wald test statistics W_λ and W_ψ**

For n ∈ {5, 10, 25} and λ = 4, simulates 1000 replications of the Wald statistics W_λ (based on CLT) and W_ψ (based on Delta method for ψ = λ²) and overlays N(0,1) on each histogram, confirming approximate standard normality for large n.

$$W_\lambda = \frac{\bar{X}_n - \lambda_0}{\widehat{\text{SE}}(\bar{X}_n)} \xrightarrow{d} \mathcal{N}(0,1), \quad W_\psi = \frac{\hat{\psi}_n - \psi_0}{\widehat{\text{SE}}(\hat{\psi}_n)} \xrightarrow{d} \mathcal{N}(0,1), \quad \widehat{\text{SE}}(\hat{\psi}_n) = 2\sqrt{\frac{\bar{X}_n^3}{n}}$$

```r
par(mfrow = c(2,3))
n_vals = c(5, 10, 25) # varied sample size
lambda = 4 # true parameter
psi = lambda^2 # true value of psi

for(n in n_vals){
  rep = 1000 # number of replications
  lambda_hat = numeric(rep) # estimate of lambda
  se_lambda_hat = numeric(rep) # estimated SE(lambda_hat)
  psi_hat = numeric(rep) # estimated SE(psi_hat)
  se_psi_hat = numeric(rep)
  w_lambda = numeric(length = rep) # Wald statistics (lambda)
  w_psi = numeric(length = rep) # Wald statistics (psi)
  
  for(i in 1:rep){
    x = rpois(n = n, lambda = lambda)
    lambda_hat[i] = mean(x)
    se_lambda_hat[i] = sqrt(lambda_hat[i]/n)
    psi_hat[i] = lambda_hat[i]^2
    se_psi_hat[i] = 2*sqrt(lambda_hat[i]^3/n)
    
    w_lambda[i] = (lambda_hat[i] - lambda)/se_lambda_hat[i]
    w_psi[i] = (psi_hat[i]-psi)/se_psi_hat[i]
  }
  
  hist(w_lambda, probability = TRUE, main = paste("n = ", n),
       xlab = expression(w[lambda]), breaks = 30,
       cex.lab = 1.5, cex.main = 1.3)
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  
  hist(w_psi, probability = TRUE, main = paste("n = ", n),
       xlab = expression(w[psi]), breaks = 30,
       cex.lab = 1.5, cex.main = 1.3)
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
}
```

---

**Multivariate Delta method – sampling distribution of the coefficient of variation**

Simulates the sampling distribution of CV = S_n/ $\overline{X}_n$ from N(μ=3, σ²=4) for n ∈ {5, 10, 20, 50, 100, 250} and overlays the Delta method normal approximation on each histogram.

$$\text{Var}\left(\frac{\hat{\sigma}}{\hat{\mu}}\right) \approx \frac{\sigma^4}{\mu^4} \cdot \frac{\sigma^2}{n} + \frac{1}{\mu^2}\,\text{Var}(\hat{\sigma}), \quad \text{Var}(\hat{\sigma}) = \sigma^2\!\left[1 - \left(\sqrt{\frac{2}{n-1}}\cdot\frac{\Gamma(n/2)}{\Gamma((n-1)/2)}\right)^2\right]$$

```r
n_vals = c(5,10,20,50,100,250) # sample size
M = 1000 # number of replications
mu = 3 # true mean
sigma = 2 # true sd
par(mfrow = c(2,3))
for(n in n_vals){
  sample_cv = numeric(length = M)
  for(i in 1:M){
    x = rnorm(n = n, mean = mu, sd = sigma)
    sample_cv[i] = sd(x)/mean(x) # coefficient of variation
  }
  hist(sample_cv, probability = TRUE, 
       xlab = expression(widehat(CV)), main = paste("n = ",n))
  delta_var = (sigma/mu)^4/n + ((1/mu)^2)*sigma^2*(1-
(sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2))^2)
  curve(dnorm(x, mean = sigma/mu, sd = sqrt(delta_var)),
        add = TRUE, col = "red", lwd = 2)
}
```
