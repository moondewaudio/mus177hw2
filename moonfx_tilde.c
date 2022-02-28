// moonfx~.cpp : Defines the exported functions for the DLL application.
//
#define PD_LONGINTTYPE long long
#include "m_pd.h"
//#include "stdafx.h"
#include <stdlib.h>
#include <math.h>

static t_class *moonfx_tilde_class;

typedef struct _moonfx_tilde {
	t_object x_obj;

	t_float oosamplerate, samplerate;
	t_float pi;

	// variables for folding
	t_float x_gain;
	t_float x_offset;

	// variables for LPF
	t_float x_cutoff;
	double in1, in2;

	// variables for delay
	t_float x_time;
	t_float x_fb;
	t_float x_mix;
	t_float peak;
	t_float *delayline;
	long delaysize;
	long wp;

	t_float f;

	t_outlet *x_out;
} t_moonfx_tilde;


// peak detector for delay compression
float moonfx_tilde_peak_detect(t_moonfx_tilde *x, float in) {
	float rect;
	rect = in;
	if (rect < 0.0) rect *= -1.0;			// absolute value
	if (x->peak < rect)
		x->peak += (rect - x->peak) * 0.3f;	// smoothing on the peak value tracking to avoid large jumps and better show average level
	else									
		x->peak *= 0.999f;					// decay peak over time to track new ones
	return(x->peak);
}

void *moonfx_tilde_new(t_floatarg f) {
	t_moonfx_tilde *x = (t_moonfx_tilde *)pd_new(moonfx_tilde_class);
	outlet_new(&x->x_obj, gensym("signal"));

	// wavefolder inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_gain"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_offset"));

	// LPF inlet
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_cutoff"));

	// delay inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_time"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_fb"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_mix"));

	// defaults
	x->x_gain = 0.5f;
	x->x_offset = 0.0f;
	x->x_cutoff = 20000.0f;
	x->x_time = 1000.0f;
	x->x_fb = 0.0f;
	x->x_mix = 0.5f;
	x->pi = 4.0 * atan(1.0);
	x->peak = 0.1f;
	x->delaysize = 1048576;
	x->delayline = (t_float *)malloc(sizeof(t_float) * x->delaysize);
	x->wp = 0;
	x->in1 = x->in2 = 0.0;

	return (void *)x;
}

void moonfx_tilde_settime(t_moonfx_tilde *x, float f) {
	if (f > 15000) f = 15000.0f;	// maximum delay time in ms set to 15 seconds
	if (f < 1) f = 1.0f;			// delay time cant be below 1 ms
	x->x_time = f;
}

void moonfx_tilde_setmix(t_moonfx_tilde *x, float f) {
	if (f > 100.0) f = 100.0f;		// signal cannot be more than 100% wet or 0% wet
	if (f < 0.0) f = 0.0f;
	x->x_mix = f / 100.0f;			// scale from 0 to 1
}

void moonfx_tilde_setfb(t_moonfx_tilde *x, float f) {
	if (f > 120.0) f = 120.0f;		// maximum feedback above 100 to allow for self oscillation in the delay
	if (f < 0.0) f = 0.0f;			
	x->x_fb = f / 100.0f;			// scale from 0 to 1.2
}

void moonfx_tilde_setcutoff(t_moonfx_tilde *x, float f) {
	if (f > 20000) f = 20000.0f;	// maximum cutoff for LPF is 20kHz
	if (f < 1) f = 1.0f;
	x->x_cutoff = f;
}

void moonfx_tilde_setgain(t_moonfx_tilde *x, float f) {
	if (f < 0) f = 0.0f;			// signal cant have negative gain, but can scale infinitely
	x->x_gain = f / 10.0f;
}

void moonfx_tilde_setoffset(t_moonfx_tilde *x, float f) {
	x->x_offset = f / 100.0f;		// dc offset on signal can go as high or low as user wants for maximum waveshaping possibilities
}

// future plans: filter after clipping before wavefolder to give more intersting shape for folder to play with
// clip -> lpf -> folder -> BPF -> phase -> delay

static t_int *moonfx_tilde_perform(t_int *w)
{
	t_moonfx_tilde *x = (t_moonfx_tilde *)(w[1]);
	t_float *in = (t_float *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
	float cascade;
	float ingain;
	double tf, c; // tan freq, coefficient for LPF

	tf = tan(x->pi * (x->x_cutoff / x->samplerate));	// convert the cutoff frequency from hz to radians
	c = (tf - 1.0) / (tf + 1.0);						// filter coefficient

	// delay line
	float rpf, frac, x0, x1;
	long rp, rpp1;

	float fbinsum, input, rectif;

	float delaysamples = (x->x_time * 0.001 * x->samplerate);
	if (delaysamples < 1.0) delaysamples = 1.0;								// delay cant be less than 1 sample or longer than the delay line
	if (delaysamples >= x->delaysize) delaysamples = x->delaysize - 1.0;	

	int sample;
	for (sample = 0; sample < n; sample++) {
		// wavefolding based on a triangle waveshape, from http://synthnotes.ucsd.edu/wp4/index.php/2019/10/31/wavefolding/
		ingain = (x->x_gain * *(in + sample)) + x->x_offset;
		cascade = cos(0.5 * x->pi * ingain) 
			- 1.0 / 9.0 * cos(1.5 * x->pi * ingain)
			+ 1.0 / 25.0 * cos(2.5 * x->pi * ingain)
			- 1.0 / 49.0 * cos(3.5 * x->pi * ingain);

		// all pass based LPF, from http://synthnotes.ucsd.edu/wp4/index.php/2019/10/14/first-order-low-pass-and-high-pass-filters/
		x->in2 = (c * cascade) + x->in1 - (c * x->in2);	// output
		x->in1 = cascade;								// current sample
		cascade = (x->in1 + x->in2)*0.5;

		// fractional delay line with mix control

		rpf = x->wp + delaysamples;
		rp = (long)rpf;
		frac = rpf - rp;

		rpp1 = rp + 1;
		if (rp >= x->delaysize) rp -= x->delaysize;		// wrap the read pointers to lie within the delay line
		if (rpp1 >= x->delaysize) rpp1 -= x->delaysize;
		x0 = *(x->delayline + rp);
		x1 = *(x->delayline + rpp1);
		
		input = cascade;
		fbinsum = x0 + (x1 - x0) * frac;								// output from the delay line, using linear interpolation to estimate values between samples
		*(out + sample) = input * (1 - x->x_mix) + fbinsum * x->x_mix;	// mix the input with the output. at 0, signal is fully dry. at 1, signal is fully wet.
		fbinsum = input + x->x_fb * fbinsum;							// new input to the delay line

		// compression
		moonfx_tilde_peak_detect(x, fbinsum);
		if (x->peak > 2.0) fbinsum *= 0.5;
		if (x->peak > 0.5)
			fbinsum *= (1.601539 + x->peak * (-1.605725 + x->peak * (0.8883899 - 0.180484 * x->peak))); // polynomial compression algorithm, from http://synthnotes.ucsd.edu/wp4/index.php/2019/11/16/delay-with-feedback/
		if (fbinsum > 1.0) fbinsum = 1.0; // set feedback gain to be within -1 and 1 to avoid any clipping in the delay
		if (fbinsum < -1.0) fbinsum = -1.0;
		*(x->delayline + x->wp) = fbinsum; // write to the delay line
		(x->wp)--;
		if (x->wp < 0) x->wp += x->delaysize; // wrap the write pointer back into the delay line if necessary
	}
	return(w + 5);
}

void moonfx_tilde_free(t_moonfx_tilde *x) {
	free(x->delayline);
}

void moonfx_tilde_dsp(t_moonfx_tilde *x, t_signal **sp) {
	x->oosamplerate = 1.0 / sp[0]->s_sr;
	x->samplerate = sp[0]->s_sr;
	dsp_add(moonfx_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void moonfx_tilde_setup(void) {
	moonfx_tilde_class = class_new(gensym("moonfx~"),
		(t_newmethod)moonfx_tilde_new,
		(t_method)moonfx_tilde_free,
		sizeof(t_moonfx_tilde),
		CLASS_DEFAULT,
		A_DEFFLOAT, 0);

	CLASS_MAINSIGNALIN(moonfx_tilde_class, t_moonfx_tilde, f);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_dsp, gensym("dsp"), A_CANT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setmix,
		gensym("x_mix"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_settime,
		gensym("x_time"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setfb,
		gensym("x_fb"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setcutoff,
		gensym("x_cutoff"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setgain,
		gensym("x_gain"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setoffset,
		gensym("x_offset"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setcutoff,
		gensym("x_cutoff"), A_DEFFLOAT, 0);


}


