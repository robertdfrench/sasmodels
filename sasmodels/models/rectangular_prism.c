double form_volume(double length_a, double b2a_ratio, double c2a_ratio);
double Iq(double q, double sld, double solvent_sld, double length_a, 
          double b2a_ratio, double c2a_ratio);

double form_volume(double length_a, double b2a_ratio, double c2a_ratio)
{
    return length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio);
}

double Fq(double q,
        const double length_a,
        const double length_b,
        const double length_c,
        double sin_theta,
        double cos_theta,
        double sin_phi,
        double cos_phi)
{
    //TODO: Some of this evaulation can be done outside the loop.
    //Need to see what is the peformance loss
    const double a_half = 0.5 * length_a;
    const double b_half = 0.5 * length_b;
    const double c_half = 0.5 * length_c;

    double termA = sinc(q * a_half * sin_theta * sin_phi);
    double termB = sinc(q * b_half * sin_theta * cos_phi);
    double termC = sinc(q * c_half * cos_theta);
    double fq = termA * termB * termC;
    return fq;
}
double Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;
   //Integration limits to use in Gaussian quadrature
    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits
    
    double outer_sum = 0.0;
    for(int i=0; i<N_POINTS_76; i++) {
        const double theta = 0.5 * ( Gauss76Z[i]*(v1b-v1a) + v1a + v1b );
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);

        double inner_sum = 0.0;
        for(int j=0; j<N_POINTS_76; j++) {
            double phi = 0.5 * ( Gauss76Z[j]*(v2b-v2a) + v2a + v2b );
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            const double AP = Fq(q,length_a, length_b, length_c,
                                sin_theta, cos_theta, sin_phi, cos_phi);
            inner_sum += Gauss76Wt[j] * AP * AP;
        }
        inner_sum = 0.5 * (v2b-v2a) * inner_sum;
        outer_sum += Gauss76Wt[i] * inner_sum * sin_theta;
    }

    double answer = 0.5*(v1b-v1a)*outer_sum;

    // Normalize by Pi (Eqn. 16). 
    // The term (ABC)^2 does not appear because it was introduced before on 
    // the definitions of termA, termB, termC.
    // The factor 2 appears because the theta integral has been defined between 
    // 0 and pi/2, instead of 0 to pi.
    answer /= M_PI_2; //Form factor P(q)

    // Multiply by contrast^2 and volume^2
    const double volume = length_a * length_b * length_c;
    answer *= square((sld-solvent_sld)*volume);

    // Convert from [1e-12 A-1] to [cm-1] 
    answer *= 1.0e-4;

    return answer;
}
