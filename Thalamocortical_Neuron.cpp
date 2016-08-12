#include "Thalamocortical_Neuron.h"
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
/*										 Synaptic currents	 										*/
/****************************************************************************************************/
double Thalamocortical_Neuron::I_AMPA(int N)  const{
	double tot_s_AMPA = 0.0;
	for (unsigned i=0; i < PY_Con.size(); i++)
		tot_s_AMPA += PY_Con[i]->s_AMPA[N];

	double I = g_AMPA * tot_s_AMPA * (V[N] - E_AMPA);
	return I;
}

double Thalamocortical_Neuron::I_NMDA(int N)  const{
	double tot_s_NMDA = 0.0;
	for (unsigned i=0; i < PY_Con.size(); i++)
		tot_s_NMDA += PY_Con[i]->s_NMDA[N];

	double I = g_NMDA * tot_s_NMDA * (V[N] - E_NMDA);
	return I;
}

double Thalamocortical_Neuron::I_GABA(int N)  const{
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


/****************************************************************************************************/
/*										RK iteration of ODEs										*/
/****************************************************************************************************/
void Thalamocortical_Neuron::set_RK(int N) {
	extern const double dt;
	V	  [N+1]=V     [0]+A[N]*dt*(1/C_m *( -(I_L(N) + I_LK(N) + I_Na(N) + I_K(N) + I_Ca(N))
											-(I_AMPA(N) + I_NMDA(N) + I_GABA(N))));
	h_Na  [N+1]=h_Na  [0]+A[N]*dt*(alpha_h_Na(N) *(1-h_Na[N]) - beta_h_Na(N) * h_Na[N]);
	m_Na  [N+1]=m_Na  [0]+A[N]*dt*(alpha_m_Na(N) *(1-m_Na[N]) - beta_m_Na(N) * m_Na[N]);
	n_K   [N+1]=n_K   [0]+A[N]*dt*(alpha_n_K (N) *(1-n_K [N]) - beta_n_K (N) * n_K [N]);
	h_Ca  [N+1]=h_Ca  [0]+A[N]*dt*(h_inf_Ca(N) - h_Ca[N])/tau_h_Ca(N);
	m_Ca  [N+1]=m_Ca  [0]+A[N]*dt*(h_inf_Ca(N) - m_Ca[N])/tau_m_Ca(N);
	h_A   [N+1]=h_A   [0]+A[N]*dt*(h_inf_A (N) - h_A [N])/tau_h_A(N);
	m_A   [N+1]=m_A   [0]+A[N]*dt*(h_inf_A (N) - m_A [N])/tau_m_A (N);
	m_h   [N+1]=m_h   [0]+A[N]*dt*(m_inf_h (N) - m_h [N])/tau_m_h(N);
	m_h2  [N+1]=m_h2  [0]+A[N]*dt*(m_inf_h (N) - m_h [N])/tau_m_h(N);
	s_AMPA[N+1]=s_AMPA[0]+A[N]*dt*(3.48/(1+exp(-(V[N]-20)/2))*(1-s_AMPA[N]) - s_AMPA[N]/tau_AMPA);
	s_NMDA[N+1]=s_NMDA[0]+A[N]*dt*(0.5 * x_NMDA[N] 			 *(1-s_NMDA[N]) - s_NMDA[N]/tau_NMDA);
	x_NMDA[N+1]=x_NMDA[0]+A[N]*dt*(3.48/(1+exp(-(V[N]-20)/2))			    - x_NMDA[N]/tau_x);
}

void Thalamocortical_Neuron::add_RK(void) {
	V	  [0] = (-3*V     [0] + 2*V     [1] + 4*V     [2] + 2*V     [3] + V 	[4])/6;
	h_Na  [0] = (-3*h_Na  [0] + 2*h_Na  [1] + 4*h_Na  [2] + 2*h_Na  [3] + h_Na  [4])/6;
	m_Na  [0] = (-3*m_Na  [0] + 2*m_Na  [1] + 4*m_Na  [2] + 2*m_Na  [3] + m_Na  [4])/6;
	n_K   [0] = (-3*n_K   [0] + 2*n_K   [1] + 4*n_K   [2] + 2*n_K   [3] + n_K	[4])/6;
	h_Ca  [0] = (-3*h_Ca  [0] + 2*h_Ca  [1] + 4*h_Ca  [2] + 2*h_Ca  [3] + h_Ca  [4])/6;
	m_Ca  [0] = (-3*m_Ca  [0] + 2*m_Ca  [1] + 4*m_Ca  [2] + 2*m_Ca  [3] + m_Ca  [4])/6;
	h_Na  [0] = (-3*h_Na  [0] + 2*h_Na  [1] + 4*h_Na  [2] + 2*h_Na  [3] + h_Na  [4])/6;
	m_Na  [0] = (-3*m_Na  [0] + 2*m_Na  [1] + 4*m_Na  [2] + 2*m_Na  [3] + m_Na  [4])/6;
	m_h   [0] = (-3*m_h   [0] + 2*m_h   [1] + 4*m_h   [2] + 2*m_h   [3] + m_h   [4])/6;
	m_h2  [0] = (-3*m_h2  [0] + 2*m_h2  [1] + 4*m_h2  [2] + 2*m_h2  [3] + m_h2  [4])/6;
	s_AMPA[0] = (-3*s_AMPA[0] + 2*s_AMPA[1] + 4*s_AMPA[2] + 2*s_AMPA[3] + s_AMPA[4])/6;
	s_NMDA[0] = (-3*s_NMDA[0] + 2*s_NMDA[1] + 4*s_NMDA[2] + 2*s_NMDA[3] + s_NMDA[4])/6;
	x_NMDA[0] = (-3*x_NMDA[0] + 2*x_NMDA[1] + 4*x_NMDA[2] + 2*x_NMDA[3] + x_NMDA[4])/6;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/
