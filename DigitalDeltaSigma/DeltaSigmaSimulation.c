#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "DeltaSigmaSimulation.h"
#include <stdint.h>

#define order 2 //Order of Delta-Sigma (1-5)
static float xv[NZEROS+1], yv[NPOLES+1];
#define SINE_POSITION_MAX 200

#define N 44100 *2

#define SAMP_SIGMA_DELTA 27 //Oversampling rate

int sinPosition = 0;
h_g_coefs *transfer_coefs;
double h_var = 0;

//Initialise the coeficients for Delta-sigma modulator
h_g_coefs *initTransferCoefs(){
	h_g_coefs *transfer_coefs = malloc(sizeof(h_g_coefs));
	transfer_coefs->a_g_coef = malloc(sizeof(double)*6);
	transfer_coefs->b_g_coef = malloc(sizeof(double)*6);
	transfer_coefs->a_h_coef = malloc(sizeof(double)*6);
	transfer_coefs->b_h_coef = malloc(sizeof(double)*6);
	transfer_coefs->u_h = malloc(sizeof(double)*6);
	transfer_coefs->u_g = malloc(sizeof(double)*6);
	int i;
	for(i=0;i<6;i++){
		transfer_coefs->a_g_coef[i]=A_G_COEF[order-1][i];
		transfer_coefs->b_g_coef[i]=B_G_COEF[order-1][i];
		transfer_coefs->a_h_coef[i]=A_H_COEF[order-1][i];
		transfer_coefs->b_h_coef[i]=B_H_COEF[order-1][i];
		transfer_coefs->u_g[i]=0;
		transfer_coefs->u_h[i]=0;
	}
	return transfer_coefs;
}

void tearDown(h_g_coefs *transfer_coefs){
	free(transfer_coefs->a_g_coef);
	free(transfer_coefs->b_g_coef);
	free(transfer_coefs->a_h_coef);
	free(transfer_coefs->b_h_coef);
	free(transfer_coefs->u_g);
	free(transfer_coefs->u_h);
	free(transfer_coefs);
}

//butterworth filter for modelling output
double butterworthFilter(double s){
    xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; 
    xv[3] = s / GAIN;
    yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; 
	yv[3] =   (xv[0] + xv[3]) + 3 * (xv[1] + xv[2])
                     + (  0.7283702835 * yv[0]) + ( -2.4154885130 * yv[1])
                     + (  2.6837110312 * yv[2]);
    return yv[3];
}

double quantize(double inp){
	if(inp<QUAN_THR_LOW){
		return QUAN_MIN;
	}
	if(inp>QUAN_THR_HIGH){
		return QUAN_MAX;
	}
	return 0.0;
}

double iir(double iirin, const double *a_coefs, const double *b_coefs, double *state_coefs){
	int i;
	double out=0;
	/*Circulate buffer*/
	for(i=order; i>0; i--)
		memcpy(&state_coefs[i],&state_coefs[i-1],sizeof(double));
	state_coefs[0] = iirin*a_coefs[0];
	/*Tap with a coefficients*/
	for(i=1;i<=order;i++)
		state_coefs[0]-=a_coefs[i]*state_coefs[i];
	/*Tap with b coefficients*/
	for(i=0;i<=order;i++)
		out+= b_coefs[i]*state_coefs[i];
	return out;
}

double modulate(double sample, double h_var,h_g_coefs *transfer_coefs,double *ret_sig,int j) {
	// double sample = sin(2000 * (2 * M_PI) * sinPosition / 44100);
	// printf("%f\n",sample);
	// if(++sinPosition >= SINE_POSITION_MAX) sinPosition=0;
	double s = sample - h_var;
	/*G filter*/
	double g = iir(s,transfer_coefs->a_g_coef,transfer_coefs->b_g_coef,transfer_coefs->u_g);
	/*Quantize output of G filter*/
	double y = quantize(g);
	//OUTPUT Y
	double filtered_s = butterworthFilter(y);
	printf("%f, %f\n",filtered_s,y);
	ret_sig[j]=filtered_s;
	/*H filter*/
	double c = iir(y,transfer_coefs->a_h_coef,transfer_coefs->b_h_coef,transfer_coefs->u_h);

	return c;
}

int main(int argc, char* argv[]){
	transfer_coefs = initTransferCoefs();
	double h_var = 0;

	int i;
	int j;

	//TEST SINE WAVE:
	int nb = 2000;
	double *buffer = malloc(sizeof(double)*nb);
	for (i = 0; i < nb; i++){
	    buffer[i] = sin(2000 * (2 * M_PI) * i / 44100);
	    //buffer[i]=0;
	}

	// int16_t buf[N] = {0};
	double *P_s = malloc(sizeof(double)*nb);
	double *P_n = malloc(sizeof(double)*nb);

	// FILE *pipein;
 //    pipein = popen("ffmpeg -i testWav.wav -f s16le -ac 1 -", "r");
 //    fread(buf, 2, N, pipein);
 //    pclose(pipein);

    // printf("INPUT:\n");
    // for(i = (44100*2)-(44100/50);i<N;i++){
    // 	printf("%f\n",((double)buf[i])/65536);
    // }

    printf("INPUT:\n");
    for(i = 0;i<nb;i++){
     	printf("%f\n",buffer[i]);
    }

    printf("OUTPUT:\n");

	//for(i = (44100*2)-(44100/50);i<N;i++){
    for(i = 0;i<nb;i++){
    	double *ret_sig = malloc(sizeof(double)*SAMP_SIGMA_DELTA);
		for(j=0;j<SAMP_SIGMA_DELTA;j++){
			//h_var = modulate(buf[i],h_var,transfer_coefs,P_n,((i-((44100*2)-(44100/50)))*SAMP_SIGMA_DELTA)+j);
			h_var = modulate(buffer[i],h_var,transfer_coefs,ret_sig,j);
			//printf("%f\n",P_n[((i-((44100*2)-(44100/50)))*SAMP_SIGMA_DELTA)+j]);
		}
		double sum = 0;
		for(j=0;j<SAMP_SIGMA_DELTA;j++){
			sum+=ret_sig[j];
		}
		P_n[i]=(double)pow(fabs(sum/SAMP_SIGMA_DELTA)-fabs(buffer[i]),2);
		//P_s[i-((44100*2)-(44100/50))]=pow(((double)buf[i])/65536,2);
		P_s[i]=(double)pow(buffer[i],2);
	}
	double TP_s = 0;
	for(i=500;i<nb;i++)
		TP_s+=P_s[i];
	double TP_n = 0;
	for(i=500;i<nb;i++)
		TP_n+=P_n[i];
	printf("S=%lf, N=%lf\n",(TP_s/nb),(TP_n/(nb)));
	printf("SNR = %f\n",(TP_s/nb)/(TP_n/(nb)));
	free(P_s);
	free(P_n);
	tearDown(transfer_coefs);
	return 0;
}