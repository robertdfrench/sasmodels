double form_volume(double radius, double thickness);

double Iq(double q, 
          double sld, double sld_solvent, double volfraction,
          double radius, double thickness);

double form_volume(double radius, double thickness)
{
    //note that for the vesicle model, the volume is ONLY the shell volume
    return M_4PI_3*(cube(radius+thickness) - cube(radius));
}

double Fq(double q,
    double sld,
    double sld_solvent,
    double volfraction,
    double radius,
    double thickness) {

    double vol,contrast,fq;

    // core first, then add in shell
    contrast = sld_solvent-sld;
    vol = M_4PI_3*cube(radius);
    fq = vol * sph_j1c(q*radius) * contrast;

    //now the shell. No volume normalization as this is done by the caller
    contrast = sld-sld_solvent;
    vol = M_4PI_3*cube(radius+thickness);
    fq += vol * sph_j1c(q*(radius+thickness)) * contrast;
    return fq;
}

double Iq(double q,
    double sld,
    double sld_solvent,
    double volfraction,
    double radius,
    double thickness)

/*
   scattering from a unilamellar vesicle.
   same functional form as the core-shell sphere, but more intuitive for
   a vesicle
*/

{

    double fq = Fq(q,sld,sld_solvent,volfraction,radius,thickness);
    //rescale to [cm-1].
    double f2 = volfraction * fq*fq*1.0e-4;
    
    return f2;
}
