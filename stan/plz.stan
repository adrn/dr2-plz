data {
    int n_stars;
    real Alambda; // extinction coefficient thing

    vector[n_stars] plx; // mas
    vector[n_stars] plx_err; // mas

    vector[n_stars] mag;
    vector[n_stars] mag_err;

    vector[n_stars] EBV;
    vector[n_stars] EBV_err;

    // vector[n_stars] FeH;
    // vector[n_stars] FeH_err;

    vector[n_stars] log10P; // [log10(day)]
    vector[n_stars] log10P_err; // [log10(day)]

    real log10P_ref; // [log10(day)]
    real FeH_ref;

}

transformed data {
}

parameters {
    real<lower=-10, upper=1> ln_s_M;
    real plx0; // mas
    real ln_sig_plx_add; // [ln(mas)]

    vector<lower=0>[n_stars] r; // true distance [pc]
    real L; // scale length in distance prior [pc]

    real a;
    // real b;
    real M_ref;
}

transformed parameters {
    vector[n_stars] DM;
    vector[n_stars] model_mag;
    vector[n_stars] model_plx_err;
    vector[n_stars] model_mag_err;

    DM = 5 * log10(r) - 5;
    model_mag = M_ref + Alambda * EBV + a * (log10P - log10P_ref) + DM;
    // b * (FeH - FeH_ref) +

    model_mag_err = sqrt(exp(2 * ln_s_M) +
                         square(a * log10P_err) +
                         square(Alambda * EBV_err) +
                         square(mag_err));
    // square(b * FeH_err) +

    model_plx_err = sqrt(square(plx_err) + exp(2 * ln_sig_plx_add));
}

model {
    ln_s_M ~ uniform(-10, 1);
    
    for (n in 1:n_stars) {
        mag[n] ~ normal(model_mag[n], model_mag_err[n]);
        plx[n] ~ normal(1000 / r[n] + plx0, model_plx_err[n]);

        // Distance hyper-prior
        target += -log(2*pow(L, 3)) + 2*log(r[n]) - r[n]/L;
    }
}
