function sat_estimate = get_satellite_struct(refined_estimate)

    global mu_Earth;

    sat_estimate.a = refined_estimate(1);
    sat_estimate.e = refined_estimate(2);
    sat_estimate.i = refined_estimate(3);
    sat_estimate.Omega = refined_estimate(4);
    sat_estimate.w = refined_estimate(5);
    sat_estimate.M0 = refined_estimate(6);
    sat_estimate.n = sqrt(mu_Earth/(sat_estimate.a^3));   % satellite mean motion [rad/s]
    sat_estimate.p = sat_estimate.a*(1 - (sat_estimate.e)^2);      % semilatus rectum [m]

end