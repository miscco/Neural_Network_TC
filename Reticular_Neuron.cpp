#include "Reticular_Neuron.h"
/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
/* Leak current */
double Reticular_Neuron::I_L	(int N) const{
	double I = g_L * (V[N] - E_L);
	return I;
}

/* Potassium leak current */
double Reticular_Neuron::I_LK	(int N) const{
	double I = g_LK * (V[N] - E_K);
	return I;
}

/* Fast sodium current */
double Reticular_Neuron::I_Na	(int N) const{
	double am_Na = 0.1*(V[N]+33)/(1-exp(-(V[N]+33)/10));
	double bm_Na = 4*exp(-(V[N]+53.7)/12);
	double m_Na  = am_Na/(am_Na+bm_Na);
	double I	 = g_Na * m_Na * m_Na * m_Na * h_Na[N] * (V[N] - E_Na);
	return I;
}

/* Fast potassium current */
double Reticular_Neuron::I_K	(int N) const{
	double I = g_K * n_K[N] * n_K[N] * n_K[N] * n_K[N] * (V[N] - E_K);
	return I;
}

/* Calcium current */
double Reticular_Neuron::I_Ca(int N)  const{
	double m_Ca = 1/(1+exp(-(V[N] + 20)/9));
	double I_Ca = g_Ca * m_Ca * m_Ca * (V[N] - E_Ca);
	return I_Ca;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Synaptic currents	 										*/
/****************************************************************************************************/
double Reticular_Neuron::I_AMPA(int N)  const{
	double tot_s_AMPA = 0.0;
	for (unsigned i=0; i < PY_Con.size(); i++)
		tot_s_AMPA += PY_Con[i]->s_AMPA[N];
	for (unsigned i=0; i < TC_Con.size(); i++)
		tot_s_AMPA += TC_Con[i]->s_AMPA[N];

	double I = g_AMPA * tot_s_AMPA * (V[N] - E_AMPA);
	return I;
}

double Reticular_Neuron::I_NMDA(int N)  const{
	double tot_s_NMDA = 0.0;
	for (unsigned i=0; i < PY_Con.size(); i++)
		tot_s_NMDA += PY_Con[i]->s_NMDA[N];
	for (unsigned i=0; i < TC_Con.size(); i++)
		tot_s_NMDA += TC_Con[i]->s_NMDA[N];

	double I = g_NMDA * tot_s_NMDA * (V[N] - E_NMDA);
	return I;
}

double Reticular_Neuron::I_GABA(int N)  const{
	double tot_s_GABA = 0.0;
	for (unsigned i=0; i < RE_Con.size(); i++)
		tot_s_GABA += RE_Con[i]->s_GABA[N];

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
double Reticular_Neuron::alpha_h_Na(int N) const{
	double alpha = 0.128*exp((17 - (V[N] + 50))/18);
	return alpha;
}

/* Sodium activation */
double Reticular_Neuron::beta_h_Na(int N) const{
	double beta = 4/(exp((40 - (V[N] + 50))/5) + 1);
	return beta;
}

/* Sodium inactivation */
double Reticular_Neuron::alpha_m_Na(int N) const{
	double alpha = 0.32*(13 - (V[N] + 50))/(exp((13 - (V[N] + 50))/4) - 1);
	return alpha;
}

/* Sodium inactivation */
double Reticular_Neuron::beta_m_Na(int N) const{
	double beta = 0.28*((V[N] + 50) - 40)/(exp(((V[N] + 50) - 40)/5) - 1);
	return beta;
}

/* Potassium activation */
double Reticular_Neuron::alpha_n_K(int N) const{
	double alpha = 0.032*(15 - (V[N] + 50))/(exp((15 - (V[N] + 50))/5) - 1);
	return alpha;
}

/* Potassium activation */
double Reticular_Neuron::beta_n_K(int N) const{
	double beta = 0.5*exp((10 - (V[N] + 50))/40);
	return beta;
}

/* Activation of T-type Ca current after Destexhe 1996 */
double Reticular_Neuron::m_inf_Ca	(int N) const{
	double Shift = 2.0;
	double m = 1.0/(1 + exp(-(V[N] + 50 + Shift)/7.4));
	return m;
}

/* Inactivation of T-type Ca current after Destexhe 1996 */
double Reticular_Neuron::h_inf_Ca	(int N) const{
	double Shift = 2.0;
	double h = 1.0/(1+exp((V[N]+78+Shift)/5.));
	return h;
}

/* Activation time constant of T-type Ca current after Destexhe 1996 */
double Reticular_Neuron::tau_m_Ca	(int N) const{
	double tau = (3.0 + 1.0/(exp((V[N] + 27.)/10.) + exp(-(V[N] + 102.)/15.)))/pow(5.0, 1.2);
	return tau;
}

/* Inactivation time constant of T-type Ca current after Destexhe 1996 */
double Reticular_Neuron::tau_h_Ca	(int N) const{
	double tau = (85.0 + 1.0/(exp((V[N] + 48.)/4.) + exp(-(V[N] + 407.)/50.)))/pow(3.0, 1.2);
	return tau;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										RK iteration of ODEs										*/
/****************************************************************************************************/
void Reticular_Neuron::set_RK(int N) {
	extern const double dt;
	V	  [N+1]=V     [0]+A[N]*dt*(1/C_m *( -(I_L(N) + I_LK(N) + I_Na(N) + I_K(N) + I_Ca(N))
											-(I_AMPA(N) + I_NMDA(N) + I_GABA(N))));
	h_Na  [N+1]=h_Na  [0]+A[N]*dt*(alpha_h_Na(N) *(1-h_Na[N]) - beta_h_Na(N) * h_Na[N]);
	m_Na  [N+1]=m_Na  [0]+A[N]*dt*(alpha_m_Na(N) *(1-m_Na[N]) - beta_m_Na(N) * m_Na[N]);
	n_K   [N+1]=n_K   [0]+A[N]*dt*(alpha_n_K (N) *(1-n_K [N]) - beta_n_K (N) * n_K [N]);
	h_Ca  [N+1]=h_Ca  [0]+A[N]*dt*(h_inf_Ca(N) - h_Ca[N])/tau_h_Ca(N);
	m_Ca  [N+1]=m_Ca  [0]+A[N]*dt*(h_inf_Ca(N) - m_Ca[N])/tau_m_Ca(N);
	s_GABA[N+1]=s_GABA[0]+A[N]*dt*(1/(1+exp(-(V[N]-20)/2))*(1-s_GABA[N]) - s_GABA[N]/tau_GABA);
}

void Reticular_Neuron::add_RK(void) {
	V	  [0] = (-3*V     [0] + 2*V     [1] + 4*V     [2] + 2*V     [3] + V 	[4])/6;
	h_Na  [0] = (-3*h_Na  [0] + 2*h_Na  [1] + 4*h_Na  [2] + 2*h_Na  [3] + h_Na  [4])/6;
	m_Na  [0] = (-3*m_Na  [0] + 2*m_Na  [1] + 4*m_Na  [2] + 2*m_Na  [3] + m_Na  [4])/6;
	n_K   [0] = (-3*n_K   [0] + 2*n_K   [1] + 4*n_K   [2] + 2*n_K   [3] + n_K	[4])/6;
	h_Ca  [0] = (-3*h_Ca  [0] + 2*h_Ca  [1] + 4*h_Ca  [2] + 2*h_Ca  [3] + h_Ca  [4])/6;
	m_Ca  [0] = (-3*m_Ca  [0] + 2*m_Ca  [1] + 4*m_Ca  [2] + 2*m_Ca  [3] + m_Ca  [4])/6;
	s_GABA[0] = (-3*s_GABA[0] + 2*s_GABA[1] + 4*s_GABA[2] + 2*s_GABA[3] + s_GABA[4])/6;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
