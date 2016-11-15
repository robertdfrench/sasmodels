double form_volume(double radius, double num_pearls);

double Iq(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld);

double linear_pearls_kernel(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld);


double form_volume(double radius, double num_pearls)
{
    // Pearl volume
    double pearl_vol = M_4PI_3 * cube(radius);
    // Return total volume
    return num_pearls * pearl_vol;;
}

double Fq(double q,
            double radius,
            double pearl_sld,
            double solvent_sld) {

    double contrast_pearl = pearl_sld - solvent_sld;
    //each volume
    double pearl_vol = M_4PI_3 * cube(radius);
    //mass
    double m_s = contrast_pearl * pearl_vol;

    //sine functions of a pearl
    double psi = sph_j1c(q * radius);
    return psi*m_s;

}
double linear_pearls_kernel(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld)
{
    double n_contrib;
    //each volume
    double pearl_vol = M_4PI_3 * cube(radius);
    //total volume
    double tot_vol = num_pearls * pearl_vol;
    //mass
    double separation = edge_sep + 2.0 * radius;


    // N pearls contribution
    int n_max = num_pearls - 1;
    n_contrib = num_pearls;
    for(int num=1; num<=n_max; num++) {
        n_contrib += (2.0*(num_pearls-num)*sinc(q*separation*num));
    }
    double fq = Fq(q,radius,pearl_sld,solvent_sld);
    // form factor for num_pearls
    double form_factor = 1.0e-4 * n_contrib * square(fq) / tot_vol;

    return form_factor;
}

double Iq(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld)
{

	double result = linear_pearls_kernel(q,
                    radius,
                    edge_sep,
                    num_pearls,
                    pearl_sld,
                    solvent_sld);

	return result;
}
