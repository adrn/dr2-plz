data {
    int n_stars;
    real Alambda; // extinction coefficient thing

    vector[n_stars] mag;
    vector[n_stars] mag_err;

    vector[n_stars] EBV;
    vector[n_stars] EBV_err;

    vector[n_stars] FeH;
    vector[n_stars] FeH_err;

    vector[n_stars] log10P;
    vector[n_stars] log10P_err;

}

transformed data {
    real log10P_ref;
    real FeH_ref;

    log10P_ref = log10(0.52854);
    FeH_ref = -1.4;
}

parameters {
    real ln_s_M;
    real varpi0;

    vector<lower=0>[n_stars] r; // true distance
    real L; // scale length in distance prior

    vector<lower=-3, upper=0>[n_stars] FeH_int; // true [Fe/H]
    vector<lower=0, upper=1>[n_stars] EBV_int; // true E(B-V)
    vector<lower=-1, upper=0>[n_stars] logP_int; // true period

    real a;
    real b;
    real M_ref;
}

transformed parameters {
    vector[n_stars] DM;
    vector[n_stars] DMprime;
    vector[n_stars] var_DM;

    DM = 5 * log10(r) - 5;
    DMprime = mag - Alambda * EBV - (a * (log10P - log10P_ref) +
                                     b * (FeH - FeH_ref) +
                                     M_ref);
    var_DM = pow(mag_err, 2) + pow(a * log10P_err, 2) + pow(b * FeH_err, 2) + 
}

model {
    for (n in 1:n_stars) {
        // True parameter uniform priors:
        FeH_int[n] ~ uniform(-3, 0);
        EBV_int[n] ~ uniform(0, 1);
        logP_int[n] ~ uniform(-1, 0);

        // Distance hyper-prior
        target += -log(2*pow(L, 3)) + 2*log(r[n]) - r[n]/L;
    }
}
