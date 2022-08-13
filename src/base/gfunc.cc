#include<gfunc.h>

void gfunc(wavefunction &psi_at_RI, wavefunction &d_psi_dr_at_RI,
                               double energy, wavefunction staticpot, fluid V_ee_0, double nuclear_charge, grid g, wavefunction wf, double R_surff)
{
    long ell, i, m_index, m_limit, index, lmindex;
    const long dimens   = g.dimens();
    const long rgridsize = g.ngps_x();
    const double dr     = g.delt_x();
    const long ind_R    = g.rindex(R_surff);
    const cplxd two_over_3delta_r = (2./3./dr,0);
    const cplxd one_over_12delta_r = (1./12./dr,0);
    const cplxd j_eps(0.0,0.0); // adds an effective window exp(-jeps t)

    wavefunction psi(rgridsize);
    wavefunction dir_psi(rgridsize);
    wavefunction vec_pot(rgridsize);
    wavefunction m2_wf(rgridsize);
    wavefunction m2_dir_wf(rgridsize);
    wavefunction a(rgridsize);
    wavefunction b(rgridsize);
    wavefunction c(rgridsize);

    const double m2_b = -10.0/6.0;
    const double d2_b = -2.0/(dr*dr);

    const double m2_a = -1.0/6.0;
    const double d2_a = 1.0/(dr*dr);

    const double m2_c = m2_a;
    const double d2_c = d2_a;
    
    // *****Upper left corrections******
    int modifiedcornerelements=1;
    const double d2_0 = -2.0/(dr*dr)*(1.0 - (nuclear_charge*dr/(12.0 - 10.0*nuclear_charge*dr)));
    const double m2_0 = -2.0*(1.0 + dr*dr/12.0*d2_0);

    psi_at_RI.nullify();
    d_psi_dr_at_RI.nullify();
  
    for (ell=0; ell<g.ngps_y(); ell++) 
    {
        if (dimens==34) 
        {
            m_limit = 0;
        } 
        else if (dimens==44) 
        {
            m_limit = ell;
        } 
        else 
        {
            cerr << "No such option implemented for qprop dim " << dimens << endl;
        };
        for (m_index=-m_limit; m_index<=m_limit; m_index++) 
        {
            for (i=0; i<rgridsize; i++) 
            {
            index=g.index(i, ell, m_index, 0); 
            psi[i]=wf[index];
            };

            //Multiply M2 and wf;
            m2_wf[0] = m2_b*psi[0] + m2_a*psi[1];
            if ((ell==0) && (modifiedcornerelements==1))
                m2_wf[0] = m2_0*psi[0] + m2_a*psi[1];

            for (i=1; i<rgridsize-1; i++)
                m2_wf[i] = m2_c*psi[i-1] + 
                    m2_b*psi[i] + m2_a*psi[i+1];

            m2_wf[rgridsize-1] = m2_c*psi[rgridsize-2]
                                    + m2_b*psi[rgridsize-1];

            // Build up the first matrix
            for (i=0; i<rgridsize; i++) 
            {
                index = g.index(i,ell,m_index,0);
                vec_pot[i] = staticpot[index] + V_ee_0[i] - energy - j_eps;
            };

            for (i=0; i<rgridsize-1; i++) a[i] = d2_a + m2_a*vec_pot[i+1];
            for (i=1; i<rgridsize; i++)   b[i] = d2_b + m2_b*vec_pot[i];
            for (i=1; i<rgridsize; i++)   c[i] = d2_c + m2_c*vec_pot[i-1];

            b[0] = d2_b + m2_b*vec_pot[0];
            if ((ell==0) && (modifiedcornerelements==1))
                b[0] = d2_0 + m2_0*vec_pot[0];

            // Solve the first matrix
            psi.solve_du(a,b,c, m2_wf, rgridsize);
        
            if (dimens==34) 
            {
                lmindex=ell;
            } else if (dimens==44) 
            {
                lmindex=ell*(ell+1)+m_index;
            } else 
            {
                cerr << "No such option implemented for qprop dim " << dimens << endl;
            };
            psi_at_RI[lmindex] = psi[ind_R];
           
            for (i=0; i<rgridsize; i++) 
            {
                a[i] = 1./6.;
                b[i] = 4./6.;
                c[i] = 1./6.;
                dir_psi[i] = (psi[i+1] - psi[i-1])*0.5/dr;
            };
            double y = sqrt(3.0)-2.;
            dir_psi[0] = psi[1]*0.5/dr + psi[0]*0.5*y/dr;
            dir_psi[rgridsize-1] = -psi[rgridsize-1]*0.5*y/dr - psi[rgridsize-2]*0.5/dr;
            b[0] = (4.+y)/6.;
            b[rgridsize-1] = (4.+y)/6.;

            dir_psi.solve_du(a,b,c, dir_psi, rgridsize);
            d_psi_dr_at_RI[lmindex] = dir_psi[ind_R];
        };
    };
};