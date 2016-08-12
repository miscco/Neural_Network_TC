#ifndef RETICULAR_NEURON_H
#define RETICULAR_NEURON_H
#pragma once
#include <cmath>
#include <vector>

#include "Pyramidal_Neuron.h"
#include "Thalamocortical_Neuron.h"

using std::vector;
class Pyramidal_Neuron;
class Inhibitory_Neuron;
class Thalamocortical_Neuron;

/****************************************************************************************************/
/*									Macro for vector initialization									*/
/****************************************************************************************************/
#ifndef _INIT
#define _INIT(x)	{x, 0.0, 0.0, 0.0, 0.0}
#endif
/****************************************************************************************************/
/*										 		end			 										*/
/****************************************************************************************************/


/****************************************************************************************************/
/*						Implementation of the pyramidal neuron after Bazhenov2002 					*/
/****************************************************************************************************/

class Reticular_Neuron {
public:
	Reticular_Neuron(void) {}

	Reticular_Neuron(vector<double> Param)
		: E_L(Param[0]), g_L(Param[1])
	{}

	/* Set strength of input current */
	void	set_Input(double I) {Input = I;}

	/* ODE functions */
	void 	set_RK		(int);
	void 	add_RK	 	(void);

private:
	/* Current functions */
	double 	I_L     (int) const;
	double 	I_LK    (int) const;
	double 	I_Na    (int) const;
	double 	I_K     (int) const;
	double  I_Ca    (int) const;
	double  I_h     (int) const;

	/* Synaptic currents */
	double 	I_AMPA  (int) const;
	double 	I_NMDA  (int) const;
	double 	I_GABA  (int) const;

	/* Gating functions */
	double 	alpha_h_Na(int) const;
	double 	alpha_m_Na(int) const;
	double 	alpha_n_K (int) const;
	double 	beta_h_Na (int) const;
	double 	beta_m_Na (int) const;
	double 	beta_n_K  (int) const;

	double  m_inf_Ca  (int) const;
	double  h_inf_Ca  (int) const;
	double  m_inf_A   (int) const;
	double  h_inf_A   (int) const;
	double  m_inf_h   (int) const;

	double  tau_m_Ca  (int) const;
	double  tau_h_Ca  (int) const;
	double  tau_m_A   (int) const;
	double  tau_h_A   (int) const;
	double  tau_m_h   (int) const;

	/* Connections (Neurons that target THIS Neuron*/
	vector<Pyramidal_Neuron*>       PY_Con;
	vector<Thalamocortical_Neuron*> TC_Con;
	vector<Reticular_Neuron*>       RE_Con;

	/* Paramter constants */
	/* Membrane conductivity */
	const int		C_m		= 1;

	/* Averaged membrane area */
	const double	A_i		= 20E-5;

	/* Reversal potentials */
	const double	E_L		= -63.8;
	const int		E_K		= -90;
	const int		E_Na	= 55;
	const int		E_Ca	= 55;

	const int		E_AMPA	= 0;
	const int		E_NMDA	= 0;
	const int		E_GABA  = -70;

	/* Channel conductivities */
	const double	g_L		= 102.5E-3;
	const double	g_LK	= 102.5E-3;
	const double	g_Na	= 35;
	const double	g_K		= 9;
	const double	g_Ca	= 35;

	const double	g_AMPA	= 2.25E-6;
	const double	g_NMDA	= 0.5E-6;
	const double	g_GABA	= 0.165E-6;

	/* Synapse time constants */
	const int		tau_GABA= 10;

	/* Calcium */
	const double    Ca_0    = 0.1;

	/* Noise parameters */
	const double	dphi	= 5E-3;

	/* Input current */
	double			Input	= 0.0;

	/* Parameters for the RK iteration */
	const vector<double> A = {0.5, 0.5, 1.0, 1.0};
	const vector<double> B = {0.75, 0.75, 0.0, 0.0};

	/* Variables of the neuron */
	vector<double> 	V		= _INIT(E_L),		/* Somatic membrane voltage			*/
					h_Na	= _INIT(0.0),		/* inactivation of Na channel		*/
					m_Na	= _INIT(0.0),		/* activation   of Na channel		*/
					n_K		= _INIT(0.0),   	/* activation 	of K  channel		*/
					h_Ca	= _INIT(0.0),		/* inactivation of Ca channel		*/
					m_Ca	= _INIT(0.0),		/* activation   of Ca channel		*/
					s_GABA	= _INIT(0.0);   	/* Fraction of open AMPA channels	*/

	/* Other neuron types that recieve input from this neuron type */
	friend class Thalamocortical_Neuron;

	friend void get_data(int counter,
						 vector<Pyramidal_Neuron>& PY,
						 vector<Inhibitory_Neuron>& IN,
						 vector<Thalamocortical_Neuron>& TC,
						 vector<Reticular_Neuron>& RE,
						 vector<double*> pData);

	friend void connectNeurons(vector<Pyramidal_Neuron>& PY,
							   vector<Inhibitory_Neuron>& IN,
							   vector<Thalamocortical_Neuron>& TC,
							   vector<Reticular_Neuron>& RE);
};
#endif // RETICULAR_NEURON_H
