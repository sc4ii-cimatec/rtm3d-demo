#include <cstdlib>
#include <string.h>
#include <iostream>
#include<math.h>
#include<time.h>
#include<Misc.hpp>

#include <chrono>
#include <vector>

using namespace std;

timepoint tic()
{
	return Clock::now();
}
timepoint toc(){
	return Clock::now();
}

float elapsed_s(timepoint t0, timepoint t1)
{
	// seconds sec = std::chrono::duration_cast<seconds>(t1 - t0);
	// return  (sec.count()*1.0);
	float sec = elapsed_ms(t0, t1)/1000.0;
	return sec;
}

float elapsed_ms(timepoint t0, timepoint t1)
{
	milliseconds sec = std::chrono::duration_cast<milliseconds>(t1 - t0);
	return  (sec.count()*1.0);
}


float max(float a, float b){
	if(a>b)
		return a;
	else
		return b;
}

float min(float a, float b){
	if(a<b)
		return a;
	else
		return b;
}


float randbetween(float v, float maxvel, float minvel, int dist, int range){
	//return randparabolicval(v, maxvel, minvel, dist, range);
	return randlinearval(v, dist, range);
 }

float randparabolicval(float v, float maxvel, float minvel, int dist, int range){
	char found = 0;
	float v_ave=0;
	float l_lim = 0.01*maxvel; // lowest allowed velocity value
	float delta = 0.01*v; // random variation limits
	float vtmp=0;
	int frate0;

	float c = v;// - (v - l_lim)*((dist)/(range));
	float a = (-c)/(range*range);
	
	int nrand0 = 0;
	int nrand1 = 0;
	srand(rand()*time(0)); 
	while (found == 0){

		nrand0 = rand()%(((int) delta)*dist + 1);
		nrand1 = rand()%2;

		frate0 = a*(dist*dist) +c;
		if (nrand1==0){
			vtmp = min(v, frate0 + nrand0);
		}else{
			vtmp = max(l_lim, frate0 - nrand0);
		}
		if(vtmp <= (maxvel)){
			found=1;
		}
	}
	if (vtmp < l_lim){
		vtmp+= l_lim;
	}
	return vtmp;
}

float randlinearval(float v, int k, int nkb){
	
	float v_ave = 0, l_lim = 300., delta = 200., vrand=0;
	v_ave = v - (v - l_lim) * (k) / (nkb - 1);
	int d = (int)(v + delta - (v_ave - delta) + 1) + v_ave - delta;
	vrand = rand()%d;
	return vrand;
}


void makeo2 (float *coef,int order){
	float h_beta, alpha1=0.0;
	float alpha2=0.0;
	float  central_term=0.0; 
	float coef_filt=0; 
	float arg=0.0; 
	float  coef_wind=0.0;
	int msign,ix; 
	float alpha = .54;
	float beta = 6.;

	h_beta = 0.5*beta;
	alpha1=2.*alpha-1.0;
	alpha2=2.*(1.0-alpha);
	central_term=0.0;

	msign=-1;

	for (ix=1; ix <= order/2; ix++){      
		msign=-msign ;            
		coef_filt = (2.*msign)/(ix*ix); 
		arg = M_PI*ix/(2.*(order/2+2));
		coef_wind=pow((alpha1+alpha2*cos(arg)*cos(arg)),h_beta);
		coef[order/2+ix] = coef_filt*coef_wind;
		central_term = central_term + coef[order/2+ix]; 
		coef[order/2-ix] = coef[order/2+ix]; 
	}
	
	coef[order/2]  = -2.*central_term;

	return; 
}

void laplacian_coefs(float * coefs, int order){
	 switch(order){
		case 2:
			coefs[0] = 1.;
			coefs[1] = -2.;
			coefs[2] = 1.;
			break;
		case 4:
			coefs[0] = -1./12.;
			coefs[1] = 4./3.;
			coefs[2] = -5./2.;
			coefs[3] = 4./3.;
			coefs[4] = -1./12.;
			break;
		case 6:
			coefs[0] = 1./90.;
			coefs[1] = -3./20.;
			coefs[2] = 3./2.;
			coefs[3] = -49./18.;
			coefs[4] = 3./2.;
			coefs[5] = -3./20.;
			coefs[6] = 1./90.;
			break;
		case 8:
			coefs[0] = -1./560.;
			coefs[1] = 8./315.;
			coefs[2] = -1./5.;
			coefs[3] = 8./5.;
			coefs[4] = -205./72.;
			coefs[5] = 8./5.;
			coefs[6] = -1./5.;
			coefs[7] = 8./315.;
			coefs[8] = -1./560.;
			break;
		default:
			makeo2(coefs,order);
	}
}

void SWAP_PTR(void ** A, void ** B){
    void * swap = *B;
	*B = *A;
	*A = swap;
}

void TO_UPPER(string &data){
	// convert string to upper case
	std::for_each(data.begin(), data.end(), [](char & c){
		c = ::toupper(c);
	});
}

void TO_LOWER(string &data){
	// convert string to upper case
	std::for_each(data.begin(), data.end(), [](char & c){
		c = ::tolower(c);
	});
}



template <>
HostBuffer_t<float> & SLICE(HostBuffer_t<float> &v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;
	HostBuffer_t<float> * vec = new HostBuffer_t<float>(first, last);
	return *vec;
}

