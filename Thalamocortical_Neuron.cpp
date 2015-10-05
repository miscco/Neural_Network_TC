#include "Thalamocortical_Neuron.h"
/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Thalamocortical_Neuron::set_RNG(void) {
    extern const double dt;
    /* Number of independent random variables */
    int N = 2;

    /* Create RNG for each stream */
    for (int i=0; i<N; ++i){
        /* Add the RNG for I_{l}*/
        MTRands.push_back({ENG(rand()), DIST (0.0, dphi*dt)});

        /* Add the RNG for I_{l,0} */
        MTRands.push_back({ENG(rand()), DIST (0.0, dt)});

        /* Get the random number for the first iteration */
        Rand_vars.push_back(MTRands[2*i]());
        Rand_vars.push_back(MTRands[2*i+1]());
    }
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 RK noise inputs 											*/
/****************************************************************************************************/
double Thalamocortical_Neuron::noise_xRK(int N, int M) const{
    return (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Thalamocortical_Neuron::noise_aRK(int M) const{
    return (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current */
double Thalamocortical_Neuron::I_L	(int N) const{
    double I = g_L * (V[N] - E_L);
    return I;
}

/* Leak current */
double Thalamocortical_Neuron::I_LK	(int N) const{
    double I = g_LK * (V[N] - E_K);
    return I;
}

/* Fast sodium current */
double Thalamocortical_Neuron::I_Na	(int N) const{
    double am_Na = 0.1*(V[N]+33)/(1-exp(-(V[N]+33)/10));
    double bm_Na = 4*exp(-(V[N]+53.7)/12);
    double m_Na  = am_Na/(am_Na+bm_Na);
    double I	 = g_Na * m_Na * m_Na * m_Na * h_Na[N] * (V[N] - E_Na);
    return I;
}

/* Fast potassium current */
double Thalamocortical_Neuron::I_K	(int N) const{
    double I = g_K * n_K[N] * n_K[N] * n_K[N] * n_K[N] * (V[N] - E_K);
    return I;
}

/* Calcium current */
double Thalamocortical_Neuron::I_Ca(int N)  const{
    double m_Ca = 1/(1+exp(-(V[N] + 20)/9));
    double I_Ca = g_Ca * m_Ca * m_Ca * (V[N] - E_Ca);
    return I_Ca;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*											Connectivity	 										*/
/****************************************************************************************************/
void Thalamocortical_Neuron::set_Connections(vector<Pyramidal_Neuron*>& PY,
                                             vector<Thalamocortical_Neuron*>& TC,
                                             vector<Reticular_Neuron*>& RE) {
    PY_Con = PY;
    TC_Con = TC;
    RE_Con = RE;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Synaptic currents	 										*/
/****************************************************************************************************/
double Thalamocortical_Neuron::I_AMPA(int N)  const{
    double tot_s_AMPA = 0.0;
    for (unsigned i=0; i<PY_Con.size(); ++i) {
        tot_s_AMPA += PY_Con[i]->s_AMPA[N];
    }
    double I = g_AMPA * tot_s_AMPA * (V[N] - E_AMPA);
    return I;
}

double Thalamocortical_Neuron::I_NMDA(int N)  const{
    double tot_s_NMDA = 0.0;
    for (unsigned i=0; i<PY_Con.size(); ++i) {
        tot_s_NMDA += PY_Con[i]->s_NMDA[N];
    }
    double I = g_NMDA * tot_s_NMDA * (V[N] - E_NMDA);
    return I;
}

double Thalamocortical_Neuron::I_GABA(int N)  const{
    double tot_s_GABA = 0.0;
    for (unsigned i=0; i<RE_Con.size(); ++i) {
        tot_s_GABA += RE_Con[i]->s_GABA[N];
    }
    double I = g_GABA * tot_s_GABA * (V[N] - E_GABA);
    return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Gating functions	 										*/
/****************************************************************************************************/
/* Sodium activation */
double Thalamocortical_Neuron::alpha_h_Na(int N) const{
    double alpha = 0.128*exp((17 - (V[N] + 50))/18);
    return alpha;
}

/* Sodium activation */
double Thalamocortical_Neuron::beta_h_Na(int N) const{
    double beta = 4/(exp((40 - (V[N] + 50))/5) + 1);
    return beta;
}

/* Sodium inactivation */
double Thalamocortical_Neuron::alpha_m_Na(int N) const{
    double alpha = 0.32*(13 - (V[N] + 50))/(exp((13 - (V[N] + 50))/4) - 1);
    return alpha;
}

/* Sodium inactivation */
double Thalamocortical_Neuron::beta_m_Na(int N) const{
    double beta = 0.28*((V[N] + 50) - 40)/(exp(((V[N] + 50) - 40)/5) - 1);
    return beta;
}

/* Potassium activation */
double Thalamocortical_Neuron::alpha_n_K(int N) const{
    double alpha = 0.032*(15 - (V[N] + 50))/(exp((15 - (V[N] + 50))/5) - 1);
    return alpha;
}

/* Potassium activation */
double Thalamocortical_Neuron::beta_n_K(int N) const{
    double beta = 0.5*exp((10 - (V[N] + 50))/40);
    return beta;
}

/* Activation of T-type Ca current after Destexhe 1996 */
double Thalamocortical_Neuron::m_inf_Ca	(int N) const{
    double m = 1.0/(1+exp(-(V[N]+59)/6.2));
    return m;
}

/* Inactivation of T-type Ca current after Destexhe 1996 */
double Thalamocortical_Neuron::h_inf_Ca	(int N) const{
    double h = 1.0/(1+exp((V[N]+83)/4.));
    return h;
}

/* Activation time constant of T-type Ca current after Destexhe 1996 */
double Thalamocortical_Neuron::tau_m_Ca	(int N) const{
    double tau = (1.0/(exp(-(V[N]+131.6)/16.7)+exp((V[N]+16.8)/18.2)) + 0.612)/pow(3.55, 1.2);
    return tau;
}

/* Inactivation time constant of T-type Ca current after Destexhe 1996 */
double Thalamocortical_Neuron::tau_h_Ca	(int N) const{
    double Shift = 2;
    double tau = (30.8 + (211.4 + exp((V[N] + Shift + 113.2)/5))/
                 (1+exp((V[N] + Shift + 84)/3.2)))/pow(3.0, 1.2);
    return tau;
}

/* Activation of A current after Destexhe 1996 */
double Thalamocortical_Neuron::m_inf_A	(int N) const{
    double m = 1.0/(1+exp(-(V[N]+60)/8.5));
    return m;
}

/* Inactivation of A current after Destexhe 1996 */
double Thalamocortical_Neuron::h_inf_A	(int N) const{
    double h = 1.0/(1+exp((V[N]+78)/6));
    return h;
}

/* Activation time constant of A current after Destexhe 1996 */
double Thalamocortical_Neuron::tau_m_A	(int N) const{
    double tau = (1.0/( exp((V[N]+35.82)/19.69)+exp(-(V[N]+79.69)/12.7) ) +0.37)/pow(3., 1.25);
    return tau;
}

/* Inactivation time constant of A current after Destexhe 1996 */
double Thalamocortical_Neuron::tau_h_A	(int N) const{
    double tau = V[N]>=-63 ? 19.0/pow(3.0, 1.25) :
                             1.0/((exp((V[N]+46.05)/5)+exp(-(V[N]+238.4)/37.45)))/pow(3.0, 1.25);
    return tau;
}

/* Activation of h current after Chen2012 */
double Thalamocortical_Neuron::m_inf_h	(int N) const{
    double h = 1/(1+exp( (V[N]+75)/5.5));
    return h;
}

/* Activation time for slow components in TC population after Chen2012 */
double Thalamocortical_Neuron::tau_m_h	(int N) const{
    double tau = (20 + 1000/(exp((V[N]+ 71.5)/14.2) + exp(-(V[N]+ 89)/11.6)));
    return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
