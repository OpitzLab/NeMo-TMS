TITLE dual-exponential model of NMDA receptors

COMMENT
Classic double-exponential model of NMDAR 
Mg++ voltage dependency from Spruston95 -> Woodhull, 1973 
Keivan Moradi 2011

--- (and now back to the original exp2syn comments) ---

Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

In the initial block we initialize the factor and total and A and B to starting values. The factor is 
defined in terms of tp, a local variable which defined the time of the peak of the function as 
determined by the tau1 and tau2.  tp is the maximum of the function exp(-t/tau2) – exp(-t/tau1).  To 
verify this for yourself, take the derivative, set it to 0 and solve for t.  The result is tp as defined 
here. Factor is the value of this function at time tp, and 1/factor is the normalization applied so 
that the peak is 1.  Then the synaptic weight determines the maximum synaptic conductance.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method. 

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp2NMDA
	NONSPECIFIC_CURRENT i
	RANGE tau1, tau2, e, i, Mg, K0, delta, wf
	THREADSAFE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
	(S)  = (siemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(J)  = (joules)
}

PARAMETER {
: Parameters Control Neurotransmitter and Voltage-dependent gating of NMDAR
	tau1 = 8.8		(ms)	<1e-9,1e9>	: Spruston95 CA1 dend [at Mg = 0 v=-80]	becarful: Mg can change these values
	tau2 = 500		(ms)
: Parameters Control voltage-dependent gating of NMDAR
	celsius 		(degC)	: actual temperature for simulation, defined in Neuron, usually about 35
: Parameters Control Mg block of NMDAR
	Mg = 1			(mM)	: external magnesium concentration from Spruston95
	K0 = 4.1		(mM)	: IC50 at 0 mV from Spruston95
	delta = 0.8 	(1)		: the electrical distance of the Mg2+ binding site from the outside of the membrane from Spruston95
: Parameter Controls Ohm's law in NMDAR
	e = -0.7		(mV)	: in CA1-CA3 region = -0.7 from Spruston95
}

CONSTANT {
	T = 273.16	(degC)
	F = 9.648e4	(coul)	: Faraday's constant (coulombs/mol)
	R = 8.315	(J/degC): universal gas constant (joules/mol/K)
	z = 2		(1)		: valency of Mg2+
}

ASSIGNED {
	v		(mV)
	dt		(ms)
	i		(nA)
	factor
	wf
}

STATE {
	A
	B
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	
	wf = 1
	Mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	i = (B - A)*Mgblock(v)*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight) {
	wf = weight*factor
	A = A + wf
	B = B + wf
}

FUNCTION Mgblock(v(mV)) {
	: from Spruston95
	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))
}