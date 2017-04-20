// Compute two-phase specific terms for RAS turbulence models: drag dissipation and stratification terms.

// Exponent of the turbulence damping term for k,epsilon and omega equations
tmfexp = -B*rhoa*6.0*Cmu/(max(beta*K,Tpsmall));

// just need the horizontal ones;
rtemp = (Ub-Ua) & g;
Urtemp = rtemp*g/magSqr(g);
Urh = (Ub-Ua) - Urtemp;

// ESD1 for horizontal, ESD4 for vertical density stratification    
// for 1D ESD1 ==0, obviously
ESD1 = (K/rhob)*(fvc::grad(alpha)) & Urh;

//         k was included in the equation implicitly
ESD3 = -2.0*K*alpha/rhob;

// stratification 
ESD4 = -(SUS/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha);
ESD4 = ESD4 + KE1*ESD1;

// the following term is implemented in epsilon equation.
ESD5 = -(SUS/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha);

// if ESD5 > 0, then it should be implemented, else, no
ESD5 = max(ESD5,scalar(0.0)*ESD5);
ESD5 = ESD5 + KE1*ESD1;

ESD = KE3*ESD3;   //need to check here!
