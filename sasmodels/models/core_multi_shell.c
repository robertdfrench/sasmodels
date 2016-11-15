double
form_volume(double core_radius, double n, double thickness[]);
double
Iq(double q, double core_sld, double core_radius,
   double solvent_sld, double num_shells, double sld[], double thickness[]);

double
form_volume(double core_radius, double n, double thickness[]) {
  double r = core_radius;
  for (int i=0; i < n; i++) {
    r += thickness[i];
  }
  return M_4PI_3 * cube(r);
}

Fq(double q, double core_sld, double core_radius,
    double solvent_sld, double num_shells, double sld[], double thickness[]) {

    const int n = (int)ceil(num_shells);
    double fq, r, last_sld;
    r = core_radius;
    last_sld = core_sld;
    fq = 0.;
    for (int i=0; i<n; i++) {
        fq += M_4PI_3 * cube(r) * (sld[i] - last_sld) * sph_j1c(q*r);
        last_sld = sld[i];
        r += thickness[i];
    }
    fq += M_4PI_3 * cube(r) * (solvent_sld - last_sld) * sph_j1c(q*r);
    return fq;
}

double
Iq(double q, double core_sld, double core_radius,
   double solvent_sld, double num_shells, double sld[], double thickness[])
{
  double fq;
  fq = Fq(q,core_sld, core_radius, solvent_sld,num_shells, sld, thickness);
  return fq * fq * 1.0e-4;
}
