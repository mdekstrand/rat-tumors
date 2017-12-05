data {
    int<lower=0> J;
    int<lower=0> n[J];
    int<lower=0> y[J];
}
parameters {
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0,upper=1> theta[J];
}
model {
    theta ~ beta(alpha, beta);
    y ~ binomial(n, theta);
}
