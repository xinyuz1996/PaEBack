phi    <- c(0.5, -0.4, 0.3, 0.2, - 0.1)
sigma <- 2

#Compute auto-covariance function up to lag n
AUTOCOV <- (sigma^2)*ts.extend::ARMA.autocov(n = 5, ar = phi)
AUTOCOV

(sigma^2)*ts.extend::ARMA.autocov(n = 5, ar = phi)
# Example
h = 3
a1h = c(1, phi[1], phi[1]^2+phi[2])
names(a1h) = paste("a1(",c(0:2),")", sep="")
a1h
A = 3*a1h[1]^2+2*a1h[2]^2+a1h[3]^2
A

# # sample covariance
# set.seed(300)
# ar_sim = arima.sim(list(order=c(5,0,0), ar=phi, sd = 1), n=100)
# data_input = ar_sim - mean(ar_sim) # Here we consider the ts with zero mean. Reason literature review.
# data_input = ar_sim # Here we consider the ts with zero mean. Reason literature review.
# cor = acf(ar_sim, type="correlation")
# cov = acf(ar_sim, type="covariance")
# cov[1:5]
# cov=as.numeric((acf(ar_sim, type="covariance"))$acf)

# calculate theoretical autocovariance
cov = AUTOCOV
Gamma <- matrix( c(cov[1:5], 
                   cov[2], cov[1:4],
                   cov[3:2], cov[1], cov[2:3],
                   cov[4:2], cov[1], cov[2],
                   cov[5:1]    ), byrow = T, nrow=5)
Gamma

M = matrix(c( rep(0,5), 1, rep(0, 4), 1, 1, rep(0,3) ), 
           byrow=T, nrow=3)
M1 = diag(5)
M2 = matrix(c(2*phi[1], 1, rep(0, 3),
              phi[2], phi[1], 1, 0, 0, 
              phi[3], 0, 0, 1, 0, 
              phi[4], 0, 0, 1, 1,
              phi[5], 0, 0, 0, phi[1]),
            byrow = T, nrow=5
            )
M3 = matrix(c(2*phi[1]^2+2*phi[2], 2*phi[1], 1, 0, 0,
              2*phi[1]*phi[2]+phi[2]+phi[3], phi[1]^2+phi[1], phi[1], 1, 0,
              2*phi[1]*phi[3]+phi[4], phi[3:1], 1,
              2*phi[1]*phi[4], phi[4], 0, phi[1]^2+phi[2], phi[1],
              2*phi[1]*phi[5], phi[5], 0, 0, phi[1]^2+phi[2] ), byrow=T, nrow=5)
M3
InvG = solve(Gamma)
M2 %*% Gamma
t(M2) %*% InvG

mat1 = t(M1) %*% InvG %*% M1 %*% Gamma
tr1 = sum( diag(mat1) )

mat2 = t(M2) %*% InvG %*% M2 %*% Gamma
tr2 = sum( diag(mat2) )


mat3 = t(M3) %*% InvG %*% M3 %*% Gamma
tr3 = sum( diag(mat3) )
A
B = tr1 + tr2 + tr3
A/B

tr1
tr2
tr3

r=0.3
(1/r-1)/(A/B)
r = 0.6
(1/r-1)/(A/B)
