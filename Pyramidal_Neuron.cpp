#include "Pyramidal_Neuron.h"
/****************************************************************************************************/
/*										 Initialization of RNG 										*/
/****************************************************************************************************/
void Pyramidal_Neuron::set_RNG(void) {
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
double Pyramidal_Neuron::noise_xRK(int N, int M) const{
	return (Rand_vars[2*M] + Rand_vars[2*M+1]/std::sqrt(3))*B[N];
}

double Pyramidal_Neuron::noise_aRK(int M) const{
	return (Rand_vars[2*M] - Rand_vars[2*M+1]*std::sqrt(3))/4;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Somatic currents */
/* Leak current */
double Pyramidal_Neuron::I_L	(int N) const{
	double I = g_L * (Vs[N] - E_L);
	return I;
}

/* Fast sodium current */
double Pyramidal_Neuron::I_Na	(int N) const{
	double am_Na = 0.1*(Vs[N]+33)/(1-exp(-(Vs[N]+33)/10));
	double bm_Na = 4*exp(-(Vs[N]+53.7)/12);
	double m_Na  = am_Na/(am_Na+bm_Na);
	double I	 = g_Na * m_Na * m_Na * m_Na * h_Na[N] * (Vs[N] - E_Na);
	return I;
}

/* Fast potassium current */
double Pyramidal_Neuron::I_K	(int N) const{
	double I = g_K * n_K[N] * n_K[N] * n_K[N] * n_K[N] * (Vs[N] - E_K);
	return I;
}

/* A-type current */
double Pyramidal_Neuron::I_A	(int N) const{
	double m_A	= 1/(1+exp(-(Vs[N]+50)/20));
	double I	= g_A * m_A * m_A * m_A * h_A[N] * (Vs[N] - E_K);
	return I;
}

/* KS-type current */
double Pyramidal_Neuron::I_KS	(int N) const{
	double I	= g_KS * m_KS[N] * (Vs[N] - E_K);
	return I;
}

/* Sodium dependent potassium current */
double Pyramidal_Neuron::I_KNa		(int N)  const{
	double w_KNa  = 0.37/(1+pow(38.7/Na[N], 3.5));
	double I_KNa  = g_KNa * w_KNa * (Vs[N] - E_K);
	return I_KNa;
}

/* Somato-dendritic leak */
double Pyramidal_Neuron::I_sd	(int N) const{
	double I	= g_sd * (Vs[N] - Vd[N]);
	return I;
}

/* Dendritic currents */
/* Calcium current */
double Pyramidal_Neuron::I_Ca(int N)  const{
	double m_Ca = 1/(1+exp(-(Vd[N] + 20)/9));
	double I_Ca = g_Ca * m_Ca * m_Ca * (Vd[N] - E_Ca);
	return I_Ca;
}

/* Calcium dependent potassium current */
double Pyramidal_Neuron::I_KCa(int N)  const{
	double m_KCa  = Ca[N]/ (Ca[N] + K_D);
	double I_KCa  = g_KCa * m_KCa *  (Vd[N] - E_K);
	return I_KCa;
}

/* Persistent potassium current */
double Pyramidal_Neuron::I_NaP(int N)  const{
	double m_NaP = 1/(1+exp(-(Vd[N]+55.7)/7.7));
	double I_NaP = g_NaP * m_NaP * m_NaP * m_NaP * (Vd[N] - E_Na);
	return I_NaP;
}

/* Inwardly rectifying potassium current */
double Pyramidal_Neuron::I_AR(int N)  const{
	double h_AR  = 1/(1+exp( (Vd[N]+75)/4));
	double I_AR  = g_AR * h_AR * (Vd[N] - E_K);
	return I_AR;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*											Connectivity	 										*/
/****************************************************************************************************/
void Pyramidal_Neuron::set_Connections(vector<Pyramidal_Neuron*>& PY,
                                       vector<Inhibitory_Neuron*>& IN,
                                       vector<Thalamocortical_Neuron*>& TC) {
    PY_Con = PY;
    IN_Con = IN;
    TC_Con = TC;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Synaptic currents	 										*/
/****************************************************************************************************/
double Pyramidal_Neuron::I_AMPA(int N)  const{
	double tot_s_AMPA = 0.0;
	for (unsigned i=0; i<PY_Con.size(); ++i) {
		tot_s_AMPA += PY_Con[i]->s_AMPA[N];
	}
	double I = g_AMPA * tot_s_AMPA * (Vd[N] - E_AMPA);
	return I;
}

double Pyramidal_Neuron::I_NMDA(int N)  const{
	double tot_s_NMDA = 0.0;
	for (unsigned i=0; i<PY_Con.size(); ++i) {
		tot_s_NMDA += PY_Con[i]->s_NMDA[N];
	}
	double I = g_NMDA * tot_s_NMDA * (Vd[N] - E_NMDA);
	return I;
}

double Pyramidal_Neuron::I_GABA(int N)  const{
	double tot_s_GABA = 0.0;
	for (unsigned i=0; i<IN_Con.size(); ++i) {
		tot_s_GABA += IN_Con[i]->s_GABA[N];
	}
	double I = g_GABA * tot_s_GABA * (Vd[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Gating functions	 										*/
/****************************************************************************************************/
/* Sodium activation */
double Pyramidal_Neuron::alpha_h_Na(int N) const{
	double alpha = 0.28 *exp(-(Vs[N] + 50)/10);
	return alpha;
}

/* Sodium activation */
double Pyramidal_Neuron::beta_h_Na(int N) const{
	double beta = 4./(1+exp(-(Vs[N] + 20)/10));
	return beta;
}
/* Potassium activation */
double Pyramidal_Neuron::alpha_n_K(int N) const{
	double alpha = 0.04*(Vs[N] + 34)/(1-exp(-(Vs[N] + 34)/10));
	return alpha;
}

/* Potassium activation */
double Pyramidal_Neuron::beta_n_K(int N) const{
	double beta = 0.5*exp(-(Vs[N] + 44)/25);
	return beta;
}

/* A_type current inactivation */
double Pyramidal_Neuron::h_A_inf(int N) const{
	double h = 1/(1+exp( (Vs[N]+80)/6));
	return h;
}

/* Non-inactivating potassium activation variable */
double Pyramidal_Neuron::m_KS_inf(int N) const{
	double m = 1/(1+exp(-(Vs[N]+34)/6.5));
	return m;
}

/* Non-inactivating potassium time constant */
double Pyramidal_Neuron::tau_m_KS(int N) const{
	double tau = 8/(exp( (Vs[N]+55)/30) + exp(-(Vs[N]+55)/30));
	return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*									 		Potassium pump	 										*/
/****************************************************************************************************/
double Pyramidal_Neuron::Na_pump		(int N) const{
	double Na_pump = R_pump*( Na[N]*Na[N]*Na[N]/(Na[N]*Na[N]*Na[N]+3375)
							 -Na_0 *Na_0 *Na_0 /(Na_0 *Na_0 *Na_0 +3375));
	return Na_pump;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										RK iteration of ODEs										*/
/****************************************************************************************************/
void Pyramidal_Neuron::set_RK(int N) {
	extern const double dt;
	Vd	  [N+1]=Vd    [0]+A[N]*dt*(1/C_m *( -(I_Ca(N) + I_KCa (N) + I_NaP(N) + I_AR(N))
											-(I_AMPA(N) + I_NMDA(N) - I_sd(N))/A_d));
	Vs	  [N+1]=Vs    [0]+A[N]*dt*(1/C_m *( -(I_L(N) + I_Na(N) + I_K(N) + I_A(N) + I_KS(N)
											  +I_KNa(N))-(I_GABA(N) + I_sd(N))/A_s));
	Ca    [N+1]=Ca    [0]+A[N]*dt*(-alpha_Ca *  A_d * I_Ca(N) -  Ca[N]/tau_Ca);
	Na    [N+1]=Na    [0]+A[N]*dt*(-alpha_Na *( A_s * I_Na(N) + A_d*I_NaP(N)) - Na_pump(N));
	h_Na  [N+1]=h_Na  [0]+A[N]*dt*(alpha_h_Na(N) *(1-h_Na[N]) - beta_h_Na(N) * h_Na[N]);
	n_K   [N+1]=n_K   [0]+A[N]*dt*(alpha_n_K (N) *(1-n_K [N]) - beta_n_K (N) * n_K [N]);
	h_A   [N+1]=h_A   [0]+A[N]*dt*(h_A_inf(N)  - h_A [N])/tau_A;
	m_KS  [N+1]=m_KS  [0]+A[N]*dt*(m_KS_inf(N) - m_KS[N])/tau_m_KS(N);
	s_AMPA[N+1]=s_AMPA[0]+A[N]*dt*(3.48/(1+exp(-(Vs[N]-20)/2))*(1-s_AMPA[N]) - s_AMPA[N]/tau_AMPA);
	s_NMDA[N+1]=s_NMDA[0]+A[N]*dt*(0.5 * x_NMDA[N] 			  *(1-s_NMDA[N]) - s_NMDA[N]/tau_NMDA);
	x_NMDA[N+1]=x_NMDA[0]+A[N]*dt*(3.48/(1+exp(-(Vs[N]-20)/2))			     - x_NMDA[N]/tau_x);
}

void Pyramidal_Neuron::add_RK(void) {
	Vd    [0] = (-3*Vd    [0] + 2*Vd    [1] + 4*Vd    [2] + 2*Vd    [3] + Vd	[4])/6;
	Vs    [0] = (-3*Vs    [0] + 2*Vs    [1] + 4*Vs    [2] + 2*Vs    [3] + Vs	[4])/6;
	Ca    [0] = (-3*Ca    [0] + 2*Ca    [1] + 4*Ca    [2] + 2*Ca    [3] + Ca	[4])/6;
	Na    [0] = (-3*Na    [0] + 2*Na    [1] + 4*Na    [2] + 2*Na    [3] + Na	[4])/6;
	h_Na  [0] = (-3*h_Na  [0] + 2*h_Na  [1] + 4*h_Na  [2] + 2*h_Na  [3] + h_Na  [4])/6;
	h_A   [0] = (-3*h_A   [0] + 2*h_A   [1] + 4*h_A   [2] + 2*h_A   [3] + h_A	[4])/6;
	n_K   [0] = (-3*n_K   [0] + 2*n_K   [1] + 4*n_K   [2] + 2*n_K   [3] + n_K	[4])/6;
	m_KS  [0] = (-3*m_KS  [0] + 2*m_KS  [1] + 4*m_KS  [2] + 2*m_KS  [3] + m_KS  [4])/6;
	s_AMPA[0] = (-3*s_AMPA[0] + 2*s_AMPA[1] + 4*s_AMPA[2] + 2*s_AMPA[3] + s_AMPA[4])/6;
	s_NMDA[0] = (-3*s_NMDA[0] + 2*s_NMDA[1] + 4*s_NMDA[2] + 2*s_NMDA[3] + s_NMDA[4])/6;
	x_NMDA[0] = (-3*x_NMDA[0] + 2*x_NMDA[1] + 4*x_NMDA[2] + 2*x_NMDA[3] + x_NMDA[4])/6;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
