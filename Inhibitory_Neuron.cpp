#include "Inhibitory_Neuron.h"
/****************************************************************************************************/
/*										 Current functions 											*/
/****************************************************************************************************/
double Inhibitory_Neuron::I_Na(int N)  const{
	double alpha = 0.5*(V[N] + 35) /(1-exp(-(V[N] + 35)/10));
	double beta  = 20*exp(-(V[N] + 60)/18);
	double m_Na  = alpha/(alpha+beta);
	double I     = g_Na * m_Na * m_Na * m_Na * h_Na[N] * (V[N] - E_Na);
	return I;
}

double Inhibitory_Neuron::I_K(int N)  const{
	double I = g_K * n_K[N] * n_K[N] * n_K[N] * n_K[N] * (V[N] - E_K);
	return I;
}

double Inhibitory_Neuron::I_L(int N)  const{
	double I = g_L * (V[N]- E_L);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Gating functions	 										*/
/****************************************************************************************************/
double Inhibitory_Neuron::alpha_h_Na(int N)  const{
	double alpha  = 0.35*exp(-(V[N] + 58)/20);
	return alpha;
}

double Inhibitory_Neuron::beta_h_Na(int N)  const{
	double beta  = 5/ (1+exp(-(V[N] + 28)/10));
	return beta;
}

double Inhibitory_Neuron::alpha_n_K(int N)  const{
	double alpha  = 0.05*(V[N] + 34)/(1-exp(-(V[N] + 34)/10));
	return alpha;
}

double Inhibitory_Neuron::beta_n_K(int N)  const{
	double beta  = 0.625*exp(-(V[N] + 44)/80);
	return beta;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										 Synaptic currents	 										*/
/****************************************************************************************************/
double Inhibitory_Neuron::I_AMPA(int N)  const{
	double tot_s_AMPA = 0.0;
	for (unsigned i=0; i<PY_Con.size(); i++)
		tot_s_AMPA += PY_Con[i]->s_AMPA[N];
	for (unsigned i=0; i<TC_Con.size(); i++)
		tot_s_AMPA += TC_Con[i]->s_AMPA[N];

	double I = g_AMPA * tot_s_AMPA * (V[N] - E_AMPA);
	return I;
}

double Inhibitory_Neuron::I_NMDA(int N)  const{
	double tot_s_NMDA = 0.0;
	for (unsigned i=0; i<PY_Con.size(); i++)
		tot_s_NMDA += PY_Con[i]->s_NMDA[N];
	for (unsigned i=0; i<TC_Con.size(); i++)
		tot_s_NMDA += TC_Con[i]->s_NMDA[N];

	double I = g_NMDA * tot_s_NMDA * (V[N] - E_NMDA);
	return I;
}

double Inhibitory_Neuron::I_GABA(int N)  const{
	double tot_s_GABA = 0.0;
	for (unsigned i=0; i<IN_Con.size(); i++)
		tot_s_GABA += IN_Con[i]->s_GABA[N];

	double I = g_GABA* tot_s_GABA * (V[N] - E_GABA);
	return I;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*										RK iteration of ODEs										*/
/****************************************************************************************************/
void Inhibitory_Neuron::set_RK(int N) {
	extern const double dt;
	V	  [N+1] =V	   [0]+A[N]*dt*(1/C_m *(-(I_L(N) + I_Na(N) + I_K(N))
											-(I_AMPA(N) + I_NMDA(N) + I_GABA(N))/A_i));
	h_Na  [N+1] =h_Na  [0]+A[N]*dt*(alpha_h_Na(N) *(1-h_Na[N]) - beta_h_Na(N) * h_Na[N]);
	n_K   [N+1] =n_K   [0]+A[N]*dt*(alpha_n_K (N) *(1-n_K [N]) - beta_n_K (N) * n_K [N]);
	s_GABA[N+1] =s_GABA[0]+A[N]*dt*(1/(1+exp(-(V[N]-20)/2))*(1-s_GABA[N]) - s_GABA[N]/tau_GABA);
}

void Inhibitory_Neuron::add_RK(void) {
	V     [0] = (-3*V     [0] + 2*V     [1] + 4*V     [2] + 2*V     [3] + V 	[4])/6;
	h_Na  [0] = (-3*h_Na  [0] + 2*h_Na  [1] + 4*h_Na  [2] + 2*h_Na  [3] + h_Na  [4])/6;
	n_K   [0] = (-3*n_K   [0] + 2*n_K   [1] + 4*n_K   [2] + 2*n_K   [3] + n_K	[4])/6;
	s_GABA[0] = (-3*s_GABA[0] + 2*s_GABA[1] + 4*s_GABA[2] + 2*s_GABA[3] + s_GABA[4])/6;
}
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/

