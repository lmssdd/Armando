#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cutil_inline.h>
#include <cudpp.h>

#define MAXP 100000
#define MAXN 21
#define MAXG 1280000

#define THREADS 256
#define ParticlesInSet 12
#define SetsInBlock 20

struct  grid {
    float oX, oY;
    float size;
    int nX, nY;
    int SN;
    int* set;
    int* nump;
    int* cell;
    int* perm;
};

struct  simulation {
    float minX, maxX;
    float minY, maxY;
    float dt;
    int tsn;
    int ssi;
    int nsi;
};


// Host Variables

int *hMaterial;
float *hPosX;
float *hPosY;
float *hVelX;
float *hVelY;
float *hDensity;
float *hEnergy;
float *hPressure;
float *hVelDotX;
float *hVelDotY;
float *hDensityDot;
float *hEnergyDot;
int *hList;
int *hHash;
int *hIndex;
int *hSetStart;
int *hSetStop;

int hPN;
float hSmooth, hMass, hSound;
int hMatType[10];
float hMatProp[10][10];
struct simulation hRun;
struct grid hGrid;

float *hPosX0;
float *hPosY0;
float *hVelX0;
float *hVelY0;
float *hDensity0;
float *hEnergy0;

CUDPPHandle hSortHandle;


// Device Variables

__device__ int *dMaterial;
__device__ float *dPosX;
__device__ float *dPosY;
__device__ float *dVelX;
__device__ float *dVelY;
__device__ float *dDensity;
__device__ float *dEnergy;
__device__ float *dPressure;
__device__ float *dVelDotX;
__device__ float *dVelDotY;
__device__ float *dDensityDot;
__device__ float *dEnergyDot;
__device__ int *dList;
__device__ int *dHash;
__device__ int *dIndex;
__device__ int *dSetStart;
__device__ int *dSetStop;

__device__ int *dIntDummy;
__device__ float *dFloatDummy;

__device__ __constant__ int dPN;
__device__ __constant__ float dSmooth, dMass, dSound;
__device__ __constant__ int dMatType[10];
__device__ __constant__ float dMatProp[10][10];
__device__ __constant__ struct simulation dRun;
__device__ struct grid dGrid;

__device__ float *dPosX0;
__device__ float *dPosY0;
__device__ float *dVelX0;
__device__ float *dVelY0;
__device__ float *dDensity0;
__device__ float *dEnergy0;


// Device code


__device__ float kernelWendland(float r, float h) {

    float q, alpha, w;
    /**
     * \brief Wendland kernel
     *
     * \date Feb 8, 2011
     * \author Luca Massidda
     */
	
	q = r / h;

    // for 3D
	//alpha = 21.0 / (16.0 * M_PI * h * h * h);
	
    // for 2D
	alpha = 7.0 / (4.0 * M_PI * h * h);
	
    w = 0.0;
    if (q < 2) {
        w = 1.0 - 0.5*q;
        w *= w;
        w *= w;
        w *= 1.0 + 2.0*q;
        w *= alpha;
    }

    return w;
}


__device__ float kernelDerivWendland(float r, float h) {

    float q, alpha, dwdr;
    /**
     * \brief Wendland kernel derivative
     *
     * \date Feb 8, 2011
     * \author Luca Massidda
     */

	q = r / h;
	
    // for 3D
	//alpha = 21.0 / (16.0 * M_PI * h * h * h);
	
    // for 2D
	alpha = 7.0 / (4.0 * M_PI * h * h);
	
    dwdr = 0;
    if (q < 2) {
        dwdr = 5.0 / 8.0 * q * pow((q - 2.0), 3) ;
        dwdr *= alpha / h;
    }

    return dwdr;
}


float kernelDerivWendlandHost(float r, float h) {

    float q, alpha, dwdr;
    /**
     * \brief Wendland kernel derivative
     *
     * \date Feb 8, 2011
     * \author Luca Massidda
     */

	q = r / h;
	
    // for 3D
	//alpha = 21.0 / (16.0 * M_PI * h * h * h);
	
    // for 2D
	alpha = 7.0 / (4.0 * M_PI * h * h);
	
    dwdr = 0;
    if (q < 2) {
        dwdr = 5.0 / 8.0 * q * pow((q - 2.0), 3) ;
        dwdr *= alpha / h;
    }

    return dwdr;
}


__device__ float kernelGauss(float r, float h) {

    float r2, q2, h2, alpha, w;//, dwdr;
    /**
     * \brief Gauss kernel
     *
     * \date Dec 21, 2010
     * \author Luca Massidda
     */

    r2 = r * r ;
    h2 = h * h;
    q2 = r2 / h2;


    //alpha = 1.0 / (pow(h, 1) * pow(3.14, 0.5));
    alpha = 1.0 / (3.14 * h2);

    w = 0.0;
    //dwdr = 0.0;

    if (q2 < 4.0) {
        w = alpha * expf(-q2);
        //dwdr = w * (-2.0 * r / h2);
    }

    return w;
}


__device__ float kernelDerivGauss(float r, float h) {

    float r2, q2, h2, alpha, w, dwdr;
    /**
     * \brief Gauss kernel
     *
     * \date Dec 21, 2010
     * \author Luca Massidda
     */

    r2 = r * r ;
    h2 = h * h;
    q2 = r2 / h2;


    //alpha = 1.0 / (pow(h, 1) * pow(3.14, 0.5));
    alpha = 1.0 / (3.14 * h2);

    w = 0.0;
    dwdr = 0.0;

    if (q2 < 4.0) {
        w = alpha * expf(-q2);
        dwdr = w * (-2.0 * r / h2);
    }

    return dwdr;
}


__device__ float pressureGas(int mat ,float rho, float u) {
    /**
     * \brief Ideal gas Equation Of State
     *
     * p = (k -1) rho u
     * c = (k(k -1) u)^0.5
     *
     * k = dMatProp[mat][1]
     * pshift = dMatProp[mat][2]
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;
//    float c;

    p = (dMatProp[mat][1] - 1.0) * rho * u;
    p += dMatProp[mat][2];

//    c = sqrtf(dMatProp[mat][1] * (dMatProp[mat][1] - 1.0) * u);

    return p;
}



__device__ float pressurePoly(int mat , float rho, float u) {
    /**
     * \brief Mie-Gruneisen polynomial Equation Of State
     *
     * p = a1 mu + a2 mu^2 + a3 mu^3 + (b0 + b1 mu) rho0 u  in compression
     * p = t1 mu + t2 mu^2 + b0 rho0 u                      in tension
     *
     * rho0 = dMatProp[mat][0];
     * a1 = dMatProp[mat][1];
     * a2 = dMatProp[mat][2];
     * a3 = dMatProp[mat][3];
     * b0 = dMatProp[mat][4];
     * b1 = dMatProp[mat][5];
     * t1 = dMatProp[mat][6];
     * t2 = dMatProp[mat][7];
     * pmin = dMatProp[mat][8];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p;
//    float c;

    mu = (rho - dMatProp[mat][0]) / dMatProp[mat][0];

    if (mu < 0)
        p = (dMatProp[mat][6] * mu + dMatProp[mat][7] * mu*mu)
            + (dMatProp[mat][4] * dMatProp[mat][0] * u);
    else
        p = (dMatProp[mat][1] * mu + dMatProp[mat][2] * mu*mu
             + dMatProp[mat][3] * mu*mu*mu)
            + ((dMatProp[mat][4] + dMatProp[mat][5] * mu)
               * dMatProp[mat][0] * u);

    if (p < dMatProp[mat][8]) p = dMatProp[mat][8];

//    c = sqrtf(dMatProp[mat][1] / rho);

    return p;
}

__device__ float pressureShock(int mat, float rho, float u) {
    /**
     * \brief Mie-Gruneisen Shock Hugoniot Equation Of State
     *
     * mu = rho / rho0 -1
     * g = g * rho0 / rho
     * ph = (rho0 c0^2 mu (1 + mu)) / (1 - (s0 - 1) * mu)^2
     * uh = 1/2 ph/rho0 * (mu / (1 + mu))
     * p = ph + g * rho * (u - uh)
     *
     * rho0 = dMatProp[mat][0];
     * c0 = dMatProp[mat][1];
     * g0 = dMatProp[mat][2];
     * s0 = dMatProp[mat][3];
     * pmin = dMatProp[mat][4];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p, ph;
//    float c;

    mu = (rho - dMatProp[mat][0]) / dMatProp[mat][0];

    ph = (dMatProp[mat][0] * powf(dMatProp[mat][1], 2) * mu*(1.0 +mu))
         / powf((1.0 - (dMatProp[mat][3] -1.0) * mu), 2);

    p = ph + dMatProp[mat][2] * dMatProp[mat][0]
        * (u - (0.5 * ph / dMatProp[mat][0] * (mu / (1.0 + mu))));

    if (p < dMatProp[mat][4]) p = dMatProp[mat][4];

//    c = dMatProp[mat][1];

    return p;
}


__device__ float pressureTait(int mat, float rho, float u) {
    /**
     * \brief Tait Equation Of State
     *
     * p = rho0 * c0 * c0 / 7.0 * (powf((rho / rho0), 7) - 1.0);
     * c = c0;
     *
     * rho0 = dMatProp[mat][0];
     * c0 = dMatProp[mat][1];
     * pmin = dMatProp[mat][2];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;
//    float c;

    p = dMatProp[mat][0] * powf(dMatProp[mat][1], 2) / 7.0
        * (powf((rho / dMatProp[mat][0]), 7) - 1.0);

    if (p < dMatProp[mat][2]) p = dMatProp[mat][2];

//    c = dMatProp[mat][1];

    return p;
}


float pressureGasHost(int mat ,float rho, float u) {
    /**
     * \brief Ideal gas Equation Of State
     *
     * p = (k -1) rho u
     * c = (k(k -1) u)^0.5
     *
     * k = dMatProp[mat][1]
     * pshift = dMatProp[mat][2]
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;
//    float c;

    p = (hMatProp[mat][1] - 1.0) * rho * u;
    p += hMatProp[mat][2];

//    c = sqrtf(hMatProp[mat][1] * (hMatProp[mat][1] - 1.0) * u);

    return p;
}



float pressurePolyHost(int mat , float rho, float u) {
    /**
     * \brief Mie-Gruneisen polynomial Equation Of State
     *
     * p = a1 mu + a2 mu^2 + a3 mu^3 + (b0 + b1 mu) rho0 u  in compression
     * p = t1 mu + t2 mu^2 + b0 rho0 u                      in tension
     *
     * rho0 = dMatProp[mat][0];
     * a1 = dMatProp[mat][1];
     * a2 = dMatProp[mat][2];
     * a3 = dMatProp[mat][3];
     * b0 = dMatProp[mat][4];
     * b1 = dMatProp[mat][5];
     * t1 = dMatProp[mat][6];
     * t2 = dMatProp[mat][7];
     * pmin = dMatProp[mat][8];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p;
//    float c;

    mu = (rho - hMatProp[mat][0]) / hMatProp[mat][0];

    if (mu < 0)
        p = (hMatProp[mat][6] * mu + hMatProp[mat][7] * mu*mu)
            + (hMatProp[mat][4] * hMatProp[mat][0] * u);
    else
        p = (hMatProp[mat][1] * mu + hMatProp[mat][2] * mu*mu
             + hMatProp[mat][3] * mu*mu*mu)
            + ((hMatProp[mat][4] + hMatProp[mat][5] * mu)
               * hMatProp[mat][0] * u);

    if (p < hMatProp[mat][8]) p = hMatProp[mat][8];

//    c = sqrtf(hMatProp[mat][1] / rho);

    return p;
}

float pressureShockHost(int mat, float rho, float u) {
    /**
     * \brief Mie-Gruneisen Shock Hugoniot Equation Of State
     *
     * mu = rho / rho0 -1
     * g = g * rho0 / rho
     * ph = (rho0 c0^2 mu (1 + mu)) / (1 - (s0 - 1) * mu)^2
     * uh = 1/2 ph/rho0 * (mu / (1 + mu))
     * p = ph + g * rho * (u - uh)
     *
     * rho0 = dMatProp[mat][0];
     * c0 = dMatProp[mat][1];
     * g0 = dMatProp[mat][2];
     * s0 = dMatProp[mat][3];
     * pmin = dMatProp[mat][4];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p, ph;
//    float c;

    mu = (rho - hMatProp[mat][0]) / hMatProp[mat][0];

    ph = (hMatProp[mat][0] * powf(hMatProp[mat][1], 2) * mu*(1.0 +mu))
         / powf((1.0 - (hMatProp[mat][3] -1.0) * mu), 2);

    p = ph + hMatProp[mat][2] * hMatProp[mat][0]
        * (u - (0.5 * ph / hMatProp[mat][0] * (mu / (1.0 + mu))));

    if (p < hMatProp[mat][4]) p = hMatProp[mat][4];

//    c = hMatProp[mat][1];

    return p;
}


float pressureTaitHost(int mat, float rho, float u) {
    /**
     * \brief Tait Equation Of State
     *
     * p = rho0 * c0 * c0 / 7.0 * (powf((rho / rho0), 7) - 1.0);
     * c = c0;
     *
     * rho0 = dMatProp[mat][0];
     * c0 = dMatProp[mat][1];
     * pmin = dMatProp[mat][2];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;
//    float c;

    p = hMatProp[mat][0] * powf(hMatProp[mat][1], 2) / 7.0
        * (powf((rho / hMatProp[mat][0]), 7) - 1.0);

    if (p < hMatProp[mat][2]) p = hMatProp[mat][2];

//    c = hMatProp[mat][1];

    return p;
}


// Global code

__global__ void balanceMassMomentumDevice(const int* dList,
	const float* dPosX, const float* dPosY,
	const float* dVelX, const float* dVelY, 
	const float* dDensity, const float* dPressure, 
	float* dDensityDot, float* dVelDotX, float* dVelDotY) {
	
    /**
     * \brief Interate particles
     *
     * \date Jan 6, 2011
     * \author Luca Massidda
     */

	int ip, il, jp;
    float iDensityDot;
    float iVelDotX, iVelDotY;
    volatile float dx, dy, dz, dr, dvr, dwdr, f;
	
    ip = threadIdx.x + blockDim.x * blockIdx.x;
	
    if (ip < dPN) {
		iDensityDot = 0.0;
		iVelDotX = 0.0;
		iVelDotY = 0.0;
		
		for (il = 0; il < MAXN; il++) {
			jp = dList[ip * MAXN + il];
			
			dx = dPosX[ip] - dPosX[jp];
			dy = dPosY[ip] - dPosY[jp];
			dz = 0.0;
			dr = sqrtf(dx * dx + dy * dy + dz * dz);
			
			if (dr < 0.1 * dSmooth) dr = 100.0 * dSmooth;
			
			//dwdr = kernelDerivGauss(dr, dSmooth);
			dwdr = kernelDerivWendland(dr, dSmooth);
			
			dvr = 0.0;
			dvr += (dPosX[ip] - dPosX[jp]) * (dVelX[ip] - dVelX[jp]);
			dvr += (dPosY[ip] - dPosY[jp]) * (dVelY[ip] - dVelY[jp]);
			
			iDensityDot += dMass * dvr * dwdr / dr;
			
			// Calculate interparticle pressure action
			f = -(dPressure[ip] + dPressure[jp])
				/ (dDensity[ip] * dDensity[jp]);
			
			iVelDotX += dMass * f * dwdr * (dPosX[ip] - dPosX[jp]) / dr;
			iVelDotY += dMass * f * dwdr * (dPosY[ip] - dPosY[jp]) / dr;
			
			// Calculate shock correction for mass
			f = dDensity[ip] - dDensity[jp];
			f *= 2.0 * dSound / (dDensity[ip] + dDensity[jp]);
			
			iDensityDot += dMass * f * dwdr;
			
			// Calculate shock correction for momentum
			if (dvr < 0) f = dvr;
			else f = 0.0;
			
			f *= dSmooth / (dr * dr + 0.01 * dSmooth * dSmooth);
			f *= 2. * dSound / (dDensity[ip] + dDensity[jp]);
			f *= 0.03;
			
			iVelDotX += dMass * f * dwdr * (dPosX[ip] - dPosX[jp]) / dr;
			iVelDotY += dMass * f * dwdr * (dPosY[ip] - dPosY[jp]) / dr;
		}
		
		dDensityDot[ip] += iDensityDot;
		dVelDotX[ip] += iVelDotX;
		dVelDotY[ip] += iVelDotY;

    }
}


void balanceMassMomentumHost(void) {
	
    /**
     * \brief Interate particles
     *
     * \date Jan 6, 2011
     * \author Luca Massidda
     */

	int ip, il, jp;
    float iDensityDot;
    float iVelDotX, iVelDotY;
    float dx, dy, dz, dr, dvr, dwdr, f;
	
    for (ip = 0; ip < hPN; ip++) {
		iDensityDot = 0.0;
		iVelDotX = 0.0;
		iVelDotY = 0.0;
		
		for (il = 0; il < MAXN; il++) {
			jp = hList[ip * MAXN + il];
			
			dx = hPosX[ip] - hPosX[jp];
			dy = hPosY[ip] - hPosY[jp];
			dz = 0.0;
			dr = sqrtf(dx * dx + dy * dy + dz * dz);
			
			if (dr < 0.1 * hSmooth) dr = 100.0 * hSmooth;
			
			//dwdr = kernelDerivGauss(dr, hSmooth);
			dwdr = kernelDerivWendlandHost(dr, hSmooth);
			
			dvr = 0.0;
			dvr += (hPosX[ip] - hPosX[jp]) * (hVelX[ip] - hVelX[jp]);
			dvr += (hPosY[ip] - hPosY[jp]) * (hVelY[ip] - hVelY[jp]);
			
			iDensityDot += hMass * dvr * dwdr / dr;
			
			// Calculate interparticle pressure action
			f = -(hPressure[ip] + hPressure[jp])
				/ (hDensity[ip] * hDensity[jp]);
			
			iVelDotX += hMass * f * dwdr * (hPosX[ip] - hPosX[jp]) / dr;
			iVelDotY += hMass * f * dwdr * (hPosY[ip] - hPosY[jp]) / dr;
			
			// Calculate shock correction for mass
			f = hDensity[ip] - hDensity[jp];
			f *= 2.0 * hSound / (hDensity[ip] + hDensity[jp]);
			
			iDensityDot += hMass * f * dwdr;
			
			// Calculate shock correction for momentum
			if (dvr < 0) f = dvr;
			else f = 0.0;
			
			f *= hSmooth / (dr * dr + 0.01 * hSmooth * hSmooth);
			f *= 2. * hSound / (hDensity[ip] + hDensity[jp]);
			f *= 0.03;
			
			iVelDotX += hMass * f * dwdr * (hPosX[ip] - hPosX[jp]) / dr;
			iVelDotY += hMass * f * dwdr * (hPosY[ip] - hPosY[jp]) / dr;
		}
		
		hDensityDot[ip] += iDensityDot;
		hVelDotX[ip] += iVelDotX;
		hVelDotY[ip] += iVelDotY;

    }
}

__global__ void balanceEnergyDevice(const float* dPressure,
		const float* dDensity, const float* dDensityDot,
		float* dEnergyDot) {

    /**
     * \brief Interate particles
     *
     * \date Jan 9, 2011
     * \author Luca Massidda
     */

    volatile int ip;
    float iPressure, iDensity, iDensityDot;
    float iEnergyDot;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < dPN) {
        iPressure = dPressure[ip];
        iDensity = dDensity[ip];
        iDensityDot = dDensityDot[ip];

        iEnergyDot = (iPressure * iDensityDot) / (iDensity * iDensity);

        dEnergyDot[ip] += iEnergyDot;
    }
}


void balanceEnergyHost(void) {

    /**
     * \brief Interate particles
     *
     * \date Jan 9, 2011
     * \author Luca Massidda
     */

    int ip;
    float iPressure, iDensity, iDensityDot;
    float iEnergyDot;

    for (ip = 0; ip < hPN; ip++) {
        iPressure = hPressure[ip];
        iDensity = hDensity[ip];
        iDensityDot = hDensityDot[ip];

        iEnergyDot = (iPressure * iDensityDot) / (iDensity * iDensity);

        hEnergyDot[ip] += iEnergyDot;
    }
}


__global__ void updateParticlesDevice(const int* dMaterial,
		const float* dVelDotX, const float* dVelDotY,
		const float* dDensityDot, const float* dEnergyDot, const float alpha, 
		const float* dPosX0, const float* dPosY0,
		const float* dVelX0, const float* dVelY0,
		const float* dDensity0, const float* dEnergy0,
		float* dPosX, float* dPosY,
		float* dVelX, float* dVelY, 
		float* dDensity, float* dEnergy, float* dPressure) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip;
    float f;
    int iMaterial;
    float iDensity, iEnergy;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < dPN) {
        f = dPosX0[ip] + alpha * (dPosX[ip] + dRun.dt * dVelX[ip] - dPosX0[ip]);
        /*
		if (f < dRun.minX)
			f += dRun.maxX - dRun.minX;
		if (f > dRun.maxX)
			f -= dRun.maxX - dRun.minX;
		*/
        dPosX[ip] = f;

        f = dPosY0[ip] + alpha * (dPosY[ip] + dRun.dt * dVelY[ip] - dPosY0[ip]);
        /*
		if (f < dRun.minY)
			f += dRun.maxY - dRun.minY;
		if (f > dRun.maxY)
			f -= dRun.maxY - dRun.minY;
		*/
        dPosY[ip] = f;

        f = dVelX0[ip] + alpha * (dVelX[ip] + dRun.dt * dVelDotX[ip] - dVelX0[ip]);
        dVelX[ip] = f;
        
        f = dVelY0[ip] + alpha * (dVelY[ip] + dRun.dt * dVelDotY[ip] - dVelY0[ip]);
        dVelY[ip] = f;

        f = dDensity0[ip] + alpha * (dDensity[ip] + dRun.dt * dDensityDot[ip] - dDensity0[ip]);
        dDensity[ip] = f;

        f = dEnergy0[ip] + alpha * (dEnergy[ip] + dRun.dt * dEnergyDot[ip] - dEnergy0[ip]);
        dEnergy[ip] = f;
        
        iMaterial = dMaterial[ip];
        
        if (iMaterial < 0) {
			dVelX[ip] = dVelX0[ip];
			dVelY[ip] = dVelY0[ip];
        }

        iMaterial = abs(iMaterial);
        iDensity = dDensity[ip];
        iEnergy = dEnergy[ip];

        switch (dMatType[iMaterial]) {
        case (1) : // IDEAL GAS EOS
            dPressure[ip] = pressureGas(iMaterial, iDensity, iEnergy);
            break;
        case (2) : // MIE-GRUNEISEN POLYNOMIAL EOS
            dPressure[ip] = pressurePoly(iMaterial, iDensity, iEnergy);
            break;
        case (3) : // MIE-GRUNEISEN SHOCK EOS
            dPressure[ip] = pressureShock(iMaterial, iDensity, iEnergy);
            break;
        case (4) : // TAIT EOS
            dPressure[ip] = pressureTait(iMaterial, iDensity, iEnergy);
            break;
        default :
            dPressure[ip] = 0.0;
        }
        
        
	}
}


void updateParticlesHost(const float alpha) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip;
    float f;
    int iMaterial;
    float iDensity, iEnergy;

    for (ip = 0; ip < hPN; ip++) {
        f = hPosX0[ip] + alpha * (hPosX[ip] + hRun.dt * hVelX[ip] - hPosX0[ip]);
        /*
		if (f < hRun.minX)
			f += hRun.maxX - hRun.minX;
		if (f > hRun.maxX)
			f -= hRun.maxX - hRun.minX;
		*/
        hPosX[ip] = f;

        f = hPosY0[ip] + alpha * (hPosY[ip] + hRun.dt * hVelY[ip] - hPosY0[ip]);
        /*
		if (f < hRun.minY)
			f += hRun.maxY - hRun.minY;
		if (f > hRun.maxY)
			f -= hRun.maxY - hRun.minY;
		*/
        hPosY[ip] = f;

        f = hVelX0[ip] + alpha * (hVelX[ip] + hRun.dt * hVelDotX[ip] - hVelX0[ip]);
        hVelX[ip] = f;
        
        f = hVelY0[ip] + alpha * (hVelY[ip] + hRun.dt * hVelDotY[ip] - hVelY0[ip]);
        hVelY[ip] = f;

        f = hDensity0[ip] + alpha * (hDensity[ip] + hRun.dt * hDensityDot[ip] - hDensity0[ip]);
        hDensity[ip] = f;

        f = hEnergy0[ip] + alpha * (hEnergy[ip] + hRun.dt * hEnergyDot[ip] - hEnergy0[ip]);
        hEnergy[ip] = f;
        
        iMaterial = hMaterial[ip];
        
        if (iMaterial < 0) {
			hVelX[ip] = hVelX0[ip];
			hVelY[ip] = hVelY0[ip];
        }

        iMaterial = abs(iMaterial);
        iDensity = hDensity[ip];
        iEnergy = hEnergy[ip];

        switch (hMatType[iMaterial]) {
        case (1) : // IDEAL GAS EOS
            hPressure[ip] = pressureGasHost(iMaterial, iDensity, iEnergy);
            break;
        case (2) : // MIE-GRUNEISEN POLYNOMIAL EOS
            hPressure[ip] = pressurePolyHost(iMaterial, iDensity, iEnergy);
            break;
        case (3) : // MIE-GRUNEISEN SHOCK EOS
            hPressure[ip] = pressureShockHost(iMaterial, iDensity, iEnergy);
            break;
        case (4) : // TAIT EOS
            hPressure[ip] = pressureTaitHost(iMaterial, iDensity, iEnergy);
            break;
        default :
            hPressure[ip] = 0.0;
        }
        
        
	}
}


__global__ void updateForcesDevice(const int* dMaterial,
		float* dVelDotX, float* dVelDotY, float* dDensityDot,
		float* dEnergyDot) {
	
    int ip;
    int iMaterial;
    float iVelDotX, iVelDotY, iDensityDot, iEnergyDot;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < dPN) {
        iVelDotX = 0.0;
        iVelDotY = 0.0;
        iDensityDot = 0.0;
        iEnergyDot = 0.0;

        iMaterial = dMaterial[ip];

        if (iMaterial > 0) iVelDotY = -9.81;

        dVelDotX[ip] = iVelDotX;
        dVelDotY[ip] = iVelDotY;
        dDensityDot[ip] = iDensityDot;
        dEnergyDot[ip] = iEnergyDot;
    }
}

void updateForcesHost(void) {
	
    int ip;
    int iMaterial;
    float iVelDotX, iVelDotY, iDensityDot, iEnergyDot;

    for (ip = 0; ip < hPN; ip++) {
        iVelDotX = 0.0;
        iVelDotY = 0.0;
        iDensityDot = 0.0;
        iEnergyDot = 0.0;

        iMaterial = hMaterial[ip];

        if (iMaterial > 0) iVelDotY = -9.81;

        hVelDotX[ip] = iVelDotX;
        hVelDotY[ip] = iVelDotY;
        hDensityDot[ip] = iDensityDot;
        hEnergyDot[ip] = iEnergyDot;
    }
}


__global__ void kerSortFloat(float* dArrayOut, const float* dArrayIn, 
	const int* dIndex) {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

	int ip;
	
	ip = threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= dPN) return;
	
	dArrayOut[ip] = dArrayIn[dIndex[ip]];
}


__global__ void kerSortInt(int* dArrayOut, const int* dArrayIn, 
	const int* dIndex) {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

	int ip;
	
	ip = threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= dPN) return;
	
	dArrayOut[ip] = dArrayIn[dIndex[ip]];
}


__global__ void updateHashDevice(const struct grid dGrid, 
	const float* dPosX, const float* dPosY, 
	int* dHash, int* dIndex) {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

	int ip, ix, iy, ic;
	
    ip = threadIdx.x + blockDim.x * blockIdx.x;
	
    if (ip < dPN) {
        ix = (int) truncf((dPosX[ip] - dGrid.oX) / dGrid.size);
        iy = (int) truncf((dPosY[ip] - dGrid.oY) / dGrid.size);
		ic = ix + iy * dGrid.nX;
		
		dHash[ip] = ic;
		dIndex[ip] = ip;
	}
}


__global__ void updateSetsDevice(int *dSetStart, int *dSetStop, 
	const int* dHash) {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */
	
	__shared__ int prevHash[THREADS];
	__shared__ int nextHash[THREADS];
	
	int ip;
	int hash;
	
    ip = threadIdx.x + blockDim.x * blockIdx.x;
    if (ip >= dPN) return;
	
	hash = dHash[ip];
	
	if (threadIdx.x < THREADS -1) prevHash[threadIdx.x +1] = hash;
	if (threadIdx.x == 0) {
		if (ip == 0) prevHash[threadIdx.x] = -1;
		else prevHash[threadIdx.x] = dHash[ip -1];
	}
	
	if (threadIdx.x > 0) nextHash[threadIdx.x -1] = hash;
	if (threadIdx.x == THREADS -1) {
		if (ip == dPN -1) nextHash[threadIdx.x] = -1;
		else nextHash[threadIdx.x] = dHash[ip +1];
	}
	
	__syncthreads();
	
	if (hash != prevHash[threadIdx.x]) dSetStart[hash] = ip;
	
	if (hash != nextHash[threadIdx.x]) dSetStop[hash] = ip +1;
	
}


__global__ void updateListDevice(int *dList, 
	const int* dSetStart, const int* dSetStop, 
	const struct grid dGrid, 
	const float* dPosX, const float* dPosY) {
	
	int ip, ic, ix, iy, i, j, jp, jc, np;
	float dx, dy, dr;
	
    // Particles list is filled
    ip = threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= dPN) return;
	
	ix = (int) ((dPosX[ip] - dGrid.oX) / dGrid.size);
	iy = (int) ((dPosY[ip] - dGrid.oY) / dGrid.size);
	ic = ix + iy * dGrid.nX;
	np = 0;
	
	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			jc = ic + i + j * dGrid.nX;
			
			for (jp = dSetStart[jc]; jp < dSetStop[jc]; jp++) {
				dx = dPosX[ip] - dPosX[jp];
				dy = dPosY[ip] - dPosY[jp];
				dr = sqrtf(dx * dx + dy * dy);
				
				if ((dr < 2.0 * dSmooth) && (np < MAXN)) {
					dList[ip * MAXN + np] = jp;
					np++;
				}
					
			}
				
		}
	}
	
	while (np < MAXN) {
		dList[ip * MAXN + np] = ip;
		np++;
	}
	
}



// Host code

int initHost() {

    hMaterial = (int *) malloc(MAXP * sizeof(int));
    hPosX = (float *) malloc(MAXP * sizeof(float));
    hPosY = (float *) malloc(MAXP * sizeof(float));
    hVelX = (float *) malloc(MAXP * sizeof(float));
    hVelY = (float *) malloc(MAXP * sizeof(float));
    hDensity = (float *) malloc(MAXP * sizeof(float));
    hEnergy = (float *) malloc(MAXP * sizeof(float));
    hPressure = (float *) malloc(MAXP * sizeof(float));
    hVelDotX = (float *) malloc(MAXP * sizeof(float));
    hVelDotY = (float *) malloc(MAXP * sizeof(float));
    hDensityDot = (float *) malloc(MAXP * sizeof(float));
    hEnergyDot = (float *) malloc(MAXP * sizeof(float));
	
    hList = (int *) malloc(MAXP * MAXN * sizeof(int));
    hHash = (int *) malloc(MAXP * sizeof(int));
    hIndex = (int *) malloc(MAXP * sizeof(int));
    hSetStart = (int *) malloc(MAXG * sizeof(int));
    hSetStop = (int *) malloc(MAXG * sizeof(int));
	
    hGrid.set = (int *) malloc(MAXG * sizeof(int));
    hGrid.nump = (int *) malloc(MAXG * sizeof(int));
    hGrid.cell = (int *) malloc(MAXG * sizeof(int));
    hGrid.perm = (int *) malloc(MAXP * sizeof(int));

    hPosX0 = (float *) malloc(MAXP * sizeof(float));
    hPosY0 = (float *) malloc(MAXP * sizeof(float));
    hVelX0 = (float *) malloc(MAXP * sizeof(float));
    hVelY0 = (float *) malloc(MAXP * sizeof(float));
    hDensity0 = (float *) malloc(MAXP * sizeof(float));
    hEnergy0 = (float *) malloc(MAXP * sizeof(float));
    
    return 0;
}


int initDevice() {

    cutilSafeCall( cudaMalloc((void**) &(dMaterial), (MAXP * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dPosX), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dPosY), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelX), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelY), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dDensity), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dEnergy), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dPressure), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelDotX), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelDotY), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dDensityDot), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dEnergyDot), (MAXP * sizeof(float))) );
    
    cutilSafeCall( cudaMalloc((void**) &(dList), (MAXP * MAXN * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dHash), (MAXP * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dIndex), (MAXP * sizeof(int))) );
    
    cutilSafeCall( cudaMalloc((void**) &(dSetStart), (MAXG * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dSetStop), (MAXG * sizeof(int))) );
    
    cutilSafeCall( cudaMalloc((void**) &(dIntDummy), (MAXP * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dFloatDummy), (MAXP * sizeof(float))) );

    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.size = hGrid.size;
    dGrid.SN = hGrid.SN;
    cutilSafeCall( cudaMalloc((void**) &(dGrid.set), (MAXG * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dGrid.nump), (MAXG * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dGrid.cell), (MAXG * sizeof(int))) );
    cutilSafeCall( cudaMalloc((void**) &(dGrid.perm), (MAXP * sizeof(int))) );

    cutilSafeCall( cudaMalloc((void**) &(dPosX0), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dPosY0), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelX0), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dVelY0), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dDensity0), (MAXP * sizeof(float))) );
    cutilSafeCall( cudaMalloc((void**) &(dEnergy0), (MAXP * sizeof(float))) );

    cutilSafeCall( cudaMemcpyToSymbol("dPN", &hPN, sizeof(int)) );
    cutilSafeCall( cudaMemcpyToSymbol("dSmooth", &hSmooth, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMass", &hMass, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dSound", &hSound, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dRun", &hRun, sizeof(struct simulation)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMatType", hMatType, 10 * sizeof(int)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMatProp", hMatProp, 100 * sizeof(float)) );

    return 0;
}


int initCUDPP() {

    CUDPPConfiguration sortConfig;
    sortConfig.algorithm = CUDPP_SORT_RADIX;
    sortConfig.datatype = CUDPP_UINT;
    sortConfig.op = CUDPP_ADD;
    sortConfig.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
    cudppPlan(&hSortHandle, sortConfig, hPN, 1, 0);
    
    return 0;
}


int copyHostToDevice() {

    cutilSafeCall( cudaMemcpy(dMaterial, hMaterial,
                              (MAXP * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dPosX, hPosX,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dPosY, hPosY,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dVelX, hVelX,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dVelY, hVelY,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dDensity, hDensity,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dEnergy, hEnergy,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dPressure, hPressure,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dVelDotX, hVelDotX,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dVelDotY, hVelDotY,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dDensityDot, hDensityDot,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dEnergyDot, hEnergyDot,
                              (MAXP * sizeof(float)), cudaMemcpyHostToDevice) );

    cutilSafeCall( cudaMemcpy(dList, hList,
                              (MAXP * MAXN * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dHash, hHash,
                              (MAXP * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dIndex, hIndex,
                              (MAXP * sizeof(int)), cudaMemcpyHostToDevice) );
                              
    cutilSafeCall( cudaMemcpy(dSetStart, hSetStart,
                              (MAXG * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dSetStop, hSetStop,
                              (MAXG * sizeof(int)), cudaMemcpyHostToDevice) );
	
    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.size = hGrid.size;
    dGrid.SN = hGrid.SN;
    cutilSafeCall( cudaMemcpy(dGrid.set, hGrid.set,
                              (MAXG * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dGrid.nump, hGrid.nump,
                              (MAXG * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dGrid.cell, hGrid.cell,
                              (MAXG * sizeof(int)), cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(dGrid.perm, hGrid.perm,
                              (MAXP * sizeof(int)), cudaMemcpyHostToDevice) );

    return 0;
}


int copyDeviceToHost() {

    cutilSafeCall( cudaMemcpy(hMaterial, dMaterial,
                              (MAXP * sizeof(int)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hPosX, dPosX,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hPosY, dPosY,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hVelX, dVelX,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hVelY, dVelY,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hDensity, dDensity,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hEnergy, dEnergy,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hPressure, dPressure,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hVelDotX, dVelDotX,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hVelDotY, dVelDotY,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hDensityDot, dDensityDot,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hEnergyDot, dEnergyDot,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToHost) );

    cutilSafeCall( cudaMemcpy(hList, dList,
                              (MAXP * MAXN * sizeof(int)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hHash, dHash,
                              (MAXP * sizeof(int)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hIndex, dIndex,
                              (MAXP * sizeof(int)), cudaMemcpyDeviceToHost) );
                              
    cutilSafeCall( cudaMemcpy(hSetStart, dSetStart,
                              (MAXG * sizeof(int)), cudaMemcpyDeviceToHost) );
    cutilSafeCall( cudaMemcpy(hSetStop, dSetStop,
                              (MAXG * sizeof(int)), cudaMemcpyDeviceToHost) );
	
    return 0;
}

int backupData() {

    cutilSafeCall( cudaMemcpy(dPosX0, dPosX,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
    cutilSafeCall( cudaMemcpy(dPosY0, dPosY,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
    cutilSafeCall( cudaMemcpy(dVelX0, dVelX,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
    cutilSafeCall( cudaMemcpy(dVelY0, dVelY,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
    cutilSafeCall( cudaMemcpy(dDensity0, dDensity,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
    cutilSafeCall( cudaMemcpy(dEnergy0, dEnergy,
                              (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );

    return 0;
}

int backupDataHost() {

    memcpy(hPosX0, hPosX, MAXP * sizeof(float));
    memcpy(hPosY0, hPosY, MAXP * sizeof(float));
    memcpy(hVelX0, hVelX, MAXP * sizeof(float));
    memcpy(hVelY0, hVelY, MAXP * sizeof(float));
    memcpy(hDensity0, hDensity, MAXP * sizeof(float));
    memcpy(hEnergy0, hEnergy, MAXP * sizeof(float));

    return 0;
}


int initRun() {

    /**
     * \brief Input run data
     *
     * Reads the input file for run data
     *
     * \date Oct 21, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    char tok[10];
    int i, m, p, pn;
    int iv;
    float fv;
    int mpn, mpp[10];

    // Open stream file
    stream = fopen("armando.run", "r");

    while (!feof(stream)) {
        sprintf(tok, " ");
        fscanf(stream, "%s", tok);

        if (strcmp(tok, "MAT") == 0) {
            fscanf(stream, "%i", &iv);
            if ((iv > 0) && (iv <= 50))
                m = iv;

            for (p = 0; p < 10; p++)
                hMatProp[m][p] = 0.0;

            if ((m > 0) && (m <= 10))
                pn = 3;
            if ((m > 10) && (m <= 20))
                pn = 9;
            if ((m > 20) && (m <= 30))
                pn = 10;
            if ((m > 30) && (m <= 40))
                pn = 5;
            if ((m > 40) && (m <= 50))
                pn = 3;

            for (p = 0; p < pn; p++) {
                fscanf(stream, "%f", &fv);
                hMatProp[m][p] = fv;
            }

            printf("Material %d\n", m);
            printf("hMatProp: \n");
            for (p = 0; p < pn; p++)
                printf(" %f\n", hMatProp[m][p]);
            printf("\n");
        }

        if (strcmp(tok, "TIME") == 0) {
            fscanf(stream, "%f", &fv);
            if (fv > 0.0)
                hRun.dt = fv;

            fscanf(stream, "%i", &iv);
            if (iv > 0)
                hRun.tsn = iv;

            fscanf(stream, "%i", &iv);
            if (iv > 0)
                hRun.ssi = iv;

            printf("Time step: %f\n", hRun.dt);
            printf("Steps: %i\n", hRun.tsn);
            printf("Save step: %i\n", hRun.ssi);
            printf("\n");
        }

        if (strcmp(tok, "LIMITS") == 0) {
            fscanf(stream, "%f", &fv);
            hRun.minX = fv;

            fscanf(stream, "%f", &fv);
            hRun.maxX = fv;

            fscanf(stream, "%f", &fv);
            hRun.minY = fv;

            fscanf(stream, "%f", &fv);
            hRun.maxY = fv;

            printf("Domain limits: \n");
            printf("X: %+e - %+e \n", hRun.minX, hRun.maxX);
            printf("Y: %+e - %+e \n", hRun.minY, hRun.maxY);
            printf("\n");
        }

        if (strcmp(tok, "MONITORS") == 0) {
            fscanf(stream, "%i", &iv);
            mpn = iv;

            for (i = 0; i < mpn; i++) {
                fscanf(stream, "%i", &iv);
                mpp[i] = iv;
            }

            printf("Monitored particles: %i \n", mpn);
            if (mpn > 0) {
                printf("Index:");
                for (i = 0; i < mpn; i++)
                    printf(" %i", mpp[i]);
                printf("\n");
                printf("\n");
            }
        }
    }

    fclose(stream);

    hSound = hSmooth / hRun.dt;

    return 0;
}

int scanData() {
    /**
     * \brief Input particle data file
     *
     * Reads particle data from a disk file
     *
     * \date Oct 20, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    int i;
    float fv1, fv2, fv3;
    int iv;

    // Stream file position
    stream = fopen("in_pos.txt", "r");
    for (i = 0; !feof(stream); i++) {
        fscanf(stream, "%e %e ", &fv1, &fv2);
        hPosX[i] = fv1;
        hPosY[i] = fv2;
    }
    fclose(stream);
    hPN = i;

    // Stream file velocity
    stream = fopen("in_vel.txt", "r");
    for (i = 0; i < hPN; i++) {
        fscanf(stream, "%e %e", &fv1, &fv2);
        hVelX[i] = fv1;
        hVelY[i] = fv2;
    }
    fclose(stream);

    // Stream file info
    stream = fopen("in_info.txt", "r");
    for (i = 0; i < hPN; i++) {
        fscanf(stream, "%i %e %e ", &iv, &fv1, &fv2);
        hMaterial[i] = iv;
        hMass = fv1;
        hSmooth = fv2;
    }
    fclose(stream);

    // Stream file field
    stream = fopen("in_field.txt", "r");
    for (i = 0; i < hPN; i++) {
        fscanf(stream, "%e %e %e ", &fv1, &fv2, &fv3);
        hDensity[i] = fv1;
        hPressure[i] = fv2;
        hEnergy[i] = fv3;
    }
    fclose(stream);

    return 0;
}

int printData() {
    /**
     * \brief Particle data file output
     *
     * Saves particle data on a disk file
     *
     * \date Oct 21, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    int i;

    // Stream file position
    stream = fopen("new_pos.txt", "w");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e\n", hPosX[i], hPosY[i]);
    fclose(stream);

    // Stream file velocity
    stream = fopen("new_vel.txt", "w");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e \n", hVelX[i], hVelY[i]);
    fclose(stream);

    // Stream file info
    stream = fopen("new_info.txt", "w");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%i %+14.8e %+14.8e \n", hMaterial[i], hMass, hSmooth);
    fclose(stream);

    // Stream file field
    stream = fopen("new_field.txt", "w");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e \n", hDensity[i], hPressure[i], hEnergy[i]);
    fclose(stream);

    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e %+14.8e \n", hDensityDot[i],
                hVelDotX[i], hVelDotY[i], hEnergyDot[i]);
    fclose(stream);

    return 0;
}


int outputCase() {
    /**
     * \brief Output Case file
     *
     * Saves ensight case file
     *
     * \date Jul 5, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    int ts;

    // Open stream file
    stream = fopen("armando.case", "w");

    fprintf(stream, "# Ensight formatted case file for Armando\n");
    fprintf(stream, "\n");
    fprintf(stream, "FORMAT\n");
    fprintf(stream, "type: ensight gold\n");
    fprintf(stream, "\n");
    fprintf(stream, "GEOMETRY\n");
    fprintf(stream, "model:    1           armando_pos_*****.geo\n");
    fprintf(stream, "\n");
    fprintf(stream, "VARIABLE\n");
    fprintf(stream, "vector per node:    1 velocity armando_vel_*****.dat\n");
    fprintf(stream, "scalar per node:    1 density  armando_rho_*****.dat\n");
    fprintf(stream, "scalar per node:    1 pressure armando_pre_*****.dat\n");
    fprintf(stream, "scalar per node:    1 energy   armando_ene_*****.dat\n");
    fprintf(stream, "\n");
    fprintf(stream, "TIME\n");
    fprintf(stream, "time set: %i\n", 1);
    fprintf(stream, "number of steps: %i\n", (hRun.tsn / hRun.ssi + 1));
    fprintf(stream, "filename start number: %i\n", 0);
    fprintf(stream, "filename increment: %i\n", 1);
    fprintf(stream, "time values:\n");

    for (ts = 0; ts <= hRun.tsn; ts++)
        if ((ts % hRun.ssi) == 0)
            fprintf(stream, "%14.8e\n", (ts * hRun.dt));

    // Close stream file
    fclose(stream);

    return 0;
}

int outputData(int ss) {
    /**
     * \brief Output Data file
     *
     * Saves ensight data file
     *
     * \date Oct 21, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    char filename[80];
    int i;

    // Stream position file
    sprintf(filename, "armando_pos_%05d.geo", ss);
    stream = fopen(filename, "w");

    fprintf(stream, "Armando output in EnSight Gold format\n");
    fprintf(stream, "EnSight 8.0.7\n");
    fprintf(stream, "node id assign\n");
    fprintf(stream, "element id assign\n");
    fprintf(stream, "extents\n");
    fprintf(stream, " 1.00000e+38-1.00000e+38\n");
    fprintf(stream, " 1.00000e+38-1.00000e+38\n");
    fprintf(stream, " 1.00000e+38-1.00000e+38\n");
    fprintf(stream, "part\n");
    fprintf(stream, "%10i\n", 1);
    fprintf(stream, "SPH particles\n");
    fprintf(stream, "coordinates\n");
    fprintf(stream, "%10i\n", hPN);

    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hPosX[i]);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hPosY[i]);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", 0.0);

    fclose(stream);

    // Stream velocity file
    sprintf(filename, "armando_vel_%05d.dat", ss);
    stream = fopen(filename, "w");

    fprintf(stream, "particle velocity in EnSight Gold format\n");
    fprintf(stream, "part\n");
    fprintf(stream, "%10i\n", 1);
    fprintf(stream, "coordinates\n");

    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hVelX[i]);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hVelY[i]);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", 0.0);

    fclose(stream);

    // Stream density file
    sprintf(filename, "armando_rho_%05d.dat", ss);
    stream = fopen(filename, "w");

    fprintf(stream, "particle density in EnSight Gold format\n");
    fprintf(stream, "part\n");
    fprintf(stream, "%10i\n", 1);
    fprintf(stream, "coordinates\n");

    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hDensity[i]);

    fclose(stream);

    // Stream pressure file
    sprintf(filename, "armando_pre_%05d.dat", ss);
    stream = fopen(filename, "w");

    fprintf(stream, "particle pressure in EnSight Gold format\n");
    fprintf(stream, "part\n");
    fprintf(stream, "%10i\n", 1);
    fprintf(stream, "coordinates\n");

    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hPressure[i]);

    fclose(stream);

    // Stream energy file
    sprintf(filename, "armando_ene_%05d.dat", ss);
    stream = fopen(filename, "w");

    fprintf(stream, "particle energy in EnSight Gold format\n");
    fprintf(stream, "part\n");
    fprintf(stream, "%10i\n", 1);
    fprintf(stream, "coordinates\n");

    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e\n", hEnergy[i]);

    fclose(stream);

    return 0;
}


int outputVTK(int ss) {
    /**
     * \brief Output Data file
     *
     * Saves vtk data file
     *
     * \date Oct 21, 2010
     * \author Luca Massidda
     */

    FILE *stream;
    char filename[80];
    int i;

    // Stream position file
    sprintf(filename, "out%05d.vtk", ss);
    stream = fopen(filename, "w");
	
	fprintf(stream, "# vtk DataFile Version 2.0\n");
    fprintf(stream, "Unstructured Grid Example\n");
    fprintf(stream, "ASCII\n");
    fprintf(stream, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(stream, "POINTS %i float\n", hPN);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e %+e %+e \n", hPosX[i], hPosY[i], 0.0);

    fprintf(stream, "CELLS %i %i \n", hPN, 2*hPN);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%i %i \n", 1, i);

    fprintf(stream, "CELL_TYPES %i \n", hPN);
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%i \n", 1);
    
    fprintf(stream, "POINT_DATA %i \n", hPN);
    
    fprintf(stream, "SCALARS density float 1 \n", hPN);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e \n", hDensity[i]);
    
    fprintf(stream, "SCALARS pressure float 1 \n", hPN);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e \n", hPressure[i]);
    
    fprintf(stream, "SCALARS energy float 1 \n", hPN);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e \n", hEnergy[i]);
    
    fprintf(stream, "VECTORS velocity float\n");
    for (i = 0; i < hPN; i++)
        fprintf(stream, "%+e %+e %+e \n", hVelX[i], hVelY[i], 0.0);
    
    fclose(stream);

    return 0;
}


void initDamBreak() {

    int i, j, m, pi;
    double rho, c0, pmin;
    double dr;

    m = 1;
    rho = 1000.;
    c0 = 50.;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.02 / 4; // x4
    pi = 0;

    for (j = 0; j <= 50 * 4 ; j++) {
        for (i = 0; i <= 100 * 4; i++) {
            hPosX[pi] = i * dr + 0.5 * dr;
            hPosY[pi] = j * dr + 0.5 * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = m;
            hDensity[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0;
            hPressure[pi] = 0.0;
            pi++;
        }
    }
    // 0 - 268   0 - 150
    /*
        for (j = 151; j <= 153; j++) {
            for (i = -3; i <= 271; i++) {
                hPosX[pi] = i * dr;
                hPosY[pi] = j * dr;

                hVelX[pi] = 0.0;
                hVelY[pi] = 0.0;
                hMaterial[pi] = -m;
                hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
                hEnergy[pi] = 0.0;
            hPressure[pi] = 0.0;
                pi++;
            }
        }
    */
    for (j = -3; j <= -1; j++) {
        for (i = -3; i <= 269 * 4 + 2; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0;
            hPressure[pi] = 0.0;
            pi++;
        }
    }

    for (j = -0; j <= 80 * 4; j++) {
        for (i = -3; i <= -1; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0;
            hPressure[pi] = 0.0;
            pi++;
        }
    }

    for (j = -0; j <= 80 * 4; j++) {
        for (i = 269 * 4; i <= 269 * 4 +2; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0;
            hPressure[pi] = 0.0;
            pi++;
        }
    }

    hPN = pi;
    hSmooth = 1.2 * dr;
    hMass = rho * dr * dr;
    hSound = c0;

    hRun.minX = -1.0;
    hRun.maxX =  6.0;
    hRun.minY = -1.0;
    hRun.maxY =  4.0;

    hRun.dt = 4.0e-4 / 4; //1.0e-3;
    hRun.tsn = 10000 * 4; //1000;
    hRun.ssi = 200 * 4;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.size = 2.0 * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;


    printf("Dam break in a box \n");
    printf("Particles: %i \n", hPN);
}


void initFree() {

    int i, j, m, pi;
    double rho, c0, pmin;
    double dr;

    m = 1;
    rho = 1000.;
    c0 = 50.;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.01; // x4
    pi = 0;

    for (j = 0; j < 100; j++) {
        for (i = 0; i < 100; i++) {
            hPosX[pi] = i * dr + 0.0 * dr;
            hPosY[pi] = j * dr + 0.0 * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = m;
            hDensity[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0;
            hPressure[pi] = 1.0;
            pi++;
        }
    }
    
    hPN = pi;
    hSmooth = 1.2 * dr;
    hMass = rho * dr * dr;
    hSound = c0;

    hRun.minX = -0.5;
    hRun.maxX =  1.5;
    hRun.minY = -0.5;
    hRun.maxY =  1.5;

    hRun.dt = 0.5e-2; //1.0e-3;
    hRun.tsn = 3; //1000;
    hRun.ssi = 1;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.size = 2.0 * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;


    printf("Freefall\n");
    printf("Particles: %i \n", hPN);
}



int iSort(int *array, int *perm, int n) {
    int i;
    static int* dummy = NULL;

    if (!dummy) dummy = (int *) malloc(MAXP * sizeof(int));

    for (i = 0; i < n; i++) dummy[i] = array[i];
    for (i = 0; i < n; i++) array[i] = dummy[perm[i]];

    return 0;
}

int fSort(float *array, int *perm, int n) {
    int i;
    static float* dummy = NULL;

    if (!dummy) dummy = (float *) malloc(MAXP * sizeof(float));

    for (i = 0; i < n; i++) dummy[i] = array[i];
    for (i = 0; i < n; i++) array[i] = dummy[perm[i]];

    return 0;
}


int sortArraysDevice(void) {
	int blocks, threads;
	
	threads = THREADS;
	blocks = (hPN + threads - 1) / threads;

    // Particles are re ordered
    
	kerSortInt <<< blocks, threads >>>
	(dIntDummy, dMaterial, dIndex);
    cutilSafeCall( cudaMemcpy(dMaterial, dIntDummy, 
		(MAXP * sizeof(int)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dPosX, dIndex);
    cutilSafeCall( cudaMemcpy(dPosX, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dPosY, dIndex);
    cutilSafeCall( cudaMemcpy(dPosY, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dVelX, dIndex);
    cutilSafeCall( cudaMemcpy(dVelX, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dVelY, dIndex);
    cutilSafeCall( cudaMemcpy(dVelY, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dDensity, dIndex);
    cutilSafeCall( cudaMemcpy(dDensity, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dEnergy, dIndex);
    cutilSafeCall( cudaMemcpy(dEnergy, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(dFloatDummy, dPressure, dIndex);
    cutilSafeCall( cudaMemcpy(dPressure, dFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
    return 0;
}


int sortArraysHost(void) {

    // Particles are re ordered
    
    iSort(hHash, hIndex, hPN);
    iSort(hMaterial, hIndex, hPN);
    fSort(hPosX, hIndex, hPN);
    fSort(hPosY, hIndex, hPN);
    fSort(hVelX, hIndex, hPN);
    fSort(hVelY, hIndex, hPN);
    fSort(hDensity, hIndex, hPN);
    fSort(hEnergy, hIndex, hPN);
    fSort(hPressure, hIndex, hPN);
    
    return 0;
}


int indexCompare(const void *a, const void *b)
{
	int c, i1, i2;
	c = 0;
	i1 = *(int*)a;
	i2 = *(int*)b;
	if (hHash[i1] < hHash[i2]) c = -1;
	if (hHash[i1] > hHash[i2]) c = 1;
  return c;
}



int updateHashHost() {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

	int ip, ix, iy, ic;
	
    for (ip = 0; ip < hPN; ip++) {
        ix = (int) ((hPosX[ip] - hGrid.oX) / hGrid.size);
        iy = (int) ((hPosY[ip] - hGrid.oY) / hGrid.size);
		ic = ix + iy * hGrid.nX;
		
		hHash[ip] = ic;
		hIndex[ip] = ip;
	}
	
	return 0;
}

int updateSetsHost() {
	
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */
	
	int ip;
	
	hSetStart[hHash[0]] = 0;
	if (hHash[0] != hHash[1]) hSetStop[hHash[0]] = 1;

    for (ip = 1; ip < hPN -1; ip++) {
		if (hHash[ip] != hHash[ip -1]) hSetStart[hHash[ip]] = ip;
		if (hHash[ip] != hHash[ip +1]) hSetStop[hHash[ip]] = ip +1;
	}
	
	if (hHash[hPN -1] != hHash[hPN -2]) hSetStart[hHash[hPN -1]] = hPN -1;
	hSetStop[hHash[hPN -1]] = hPN;
	
	return 0;
}


int updateListHost(void) {
	int ip, ic, ix, iy, il, i, j, jp, jc, np;
	float dx, dy, dr;
	
    // Particles list is filled
    for (ip = 0; ip < hPN; ip++) {
		for (il = 0; il < MAXN; il++) {
			hList[ip * MAXN + il] = ip;
		}
		
        ix = (int) ((hPosX[ip] - hGrid.oX) / hGrid.size);
        iy = (int) ((hPosY[ip] - hGrid.oY) / hGrid.size);
		ic = ix + iy * hGrid.nX;
		
		np = 0;
        for (j = -1; j <= 1; j++) {
            for (i = -1; i <= 1; i++) {
				jc = ic + i + j * hGrid.nX;
				
				for (jp = hSetStart[jc]; jp < hSetStop[jc]; jp++) {
					dx = hPosX[ip] - hPosX[jp];
                    dy = hPosY[ip] - hPosY[jp];
                    dr = sqrtf(dx * dx + dy * dy);
					
					if ((dr < 2.0 * hSmooth) && (np < MAXN)) {
						hList[ip * MAXN + np] = jp;
						np++;
					}
				}
			}
		}
		
		
	}
	
	return 0;
}


int neighbourListDevice() {
	int blocks, threads;
	
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
	updateHashDevice <<< blocks, threads >>>
	(dGrid, dPosX, dPosY, dHash, dIndex);
	
	cudppSort(hSortHandle, dHash, dIndex, 24, hPN);
	
	sortArraysDevice();
	
	cutilSafeCall( cudaMemset(dSetStart, 0, dGrid.nX*dGrid.nY * sizeof(int)));
	cutilSafeCall( cudaMemset(dSetStop, 0, dGrid.nX*dGrid.nY * sizeof(int)));
	
	updateSetsDevice <<< blocks, threads >>>
	(dSetStart, dSetStop, dHash);
	
	updateListDevice <<< blocks, threads >>>
	(dList, dSetStart, dSetStop, dGrid, dPosX, dPosY);
	
	return 0;
}

int neighbourListHost() {
	
	updateHashHost();
	
	qsort(hIndex, hPN, sizeof(int), indexCompare);
	
	sortArraysHost();
	
	updateSetsHost();
	
	updateListHost();
	
	return 0;
}

int RKstepDevice(float alpha) {
	int blocks, threads;
	
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
	// External forces
	updateForcesDevice <<< blocks, threads >>>
	(dMaterial, dVelDotX, dVelDotY, dDensityDot, dEnergyDot);
	
	// Calculate particle interactions
	balanceMassMomentumDevice <<< blocks, threads >>>
	(dList, dPosX, dPosY, dVelX, dVelY, dDensity, dPressure, 
	dDensityDot, dVelDotX, dVelDotY);
	
	balanceEnergyDevice <<< blocks, threads >>>
	(dPressure, dDensity, dDensityDot, dEnergyDot);
	
	// Update particles
	updateParticlesDevice  <<< blocks, threads >>> 
	(dMaterial, dVelDotX, dVelDotY, dDensityDot, dEnergyDot, alpha,
	dPosX0, dPosY0, dVelX0, dVelY0, dDensity0, dEnergy0,
	dPosX, dPosY, dVelX, dVelY, dDensity, dEnergy, dPressure);
	
	return 0;
}

int RKstepHost(float alpha) {
	
	// External forces
	updateForcesHost();
		
	// Calculate particle interactions
	balanceMassMomentumHost();
	
	balanceEnergyHost();
	
	// Update particles
	updateParticlesHost(alpha);
	
	return 0;
}

int RKintegrateDevice(void) {

    /**
     * \brief Runge Kutta 3rd order time integration
     *
     * Integrate the Navier Stokes equations in time with the
     * Total Variation Diminishing Runge-Kutta algorithm of the 3rd order
     *
     * \date Dec 20, 2010
     * \author Luca Massidda
     */

    int ts;

    // Output
    //outputCase();
	
	copyHostToDevice();

    // TIME CYCLE
    for (ts = 0; ts <= hRun.tsn; ts++) {
		
        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
			copyDeviceToHost();
            printData();
            //outputData(ts / hRun.ssi);
            outputVTK(ts / hRun.ssi);
        }
		
		// Calculate neighbour list
		neighbourListDevice();
		
		/*
		copyDeviceToHost();
		updateHash();
		sortArraysHost();
		for (int i = 0; i < hGrid.nX*hGrid.nY; i++) {
			hSetStart[i] = 0;
			hSetStop[i] = 0;
		}
		updateGrid();
		updateList();
		copyHostToDevice();
		*/
		
		// Save initial condition
		backupData();
		
        // Step 1
		RKstepDevice(1.0);
		
        // Step 2
		RKstepDevice(1.0 / 4.0);
		
        // Step 3
		RKstepDevice(2.0 / 3.0);
	}
	
	cutilSafeCall( cudaThreadExit() );
	
	return 0;
}

int RKintegrateHost(void) {

    /**
     * \brief Runge Kutta 3rd order time integration
     *
     * Integrate the Navier Stokes equations in time with the
     * Total Variation Diminishing Runge-Kutta algorithm of the 3rd order
     *
     * \date Dec 20, 2010
     * \author Luca Massidda
     */

    int ts;
	
    // Output
    //outputCase();
	
    // TIME CYCLE
    for (ts = 0; ts <= hRun.tsn; ts++) {
		
        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
            printData();
            //outputData(ts / hRun.ssi);
            outputVTK(ts / hRun.ssi);
        }
        
		// Calculate neighbour list
		neighbourListHost();
		
		// Save initial condition
		backupDataHost();
		
        // Step 1
		RKstepHost(1.0);
		
        // Step 2
		RKstepHost(1.0 / 4.0);
		
        // Step 3
		RKstepHost(2.0 / 3.0);
	}
	
	return 0;
}

int check() {
    float *rVelDotX, *rVelDotY, *rDensityDot;
    float sumVelDotX, maxVelDotX;
    float sumVelDotY, maxVelDotY;
    float sumDensityDot, maxDensityDot;
    FILE *stream;
    
    rVelDotX = (float *) malloc(MAXP * sizeof(float));
    rVelDotY = (float *) malloc(MAXP * sizeof(float));
    rDensityDot = (float *) malloc(MAXP * sizeof(float));
	
	for (int i = 0; i < hPN; i++) {
		hPressure[i] = 1000.0 * (sin(10*M_PI * hPosX[i]) + hPosY[i]*hPosY[i]);
		hVelX[i] = sin(10*M_PI * hPosX[i]);
		hVelY[i] = hPosY[i]*hPosY[i];
	}
	
	int blocks = (hPN + THREADS - 1) / THREADS;
	copyHostToDevice();
	
	neighbourListDevice();
	
	// Calculate particle interactions
	balanceMassMomentumDevice <<< blocks, THREADS >>>
	(dList, dPosX, dPosY, dVelX, dVelY, dDensity, dPressure, 
	dDensityDot, dVelDotX, dVelDotY);
	
	copyDeviceToHost();
	
	for (int i = 0; i < hPN; i++) {
		rVelDotX[i] = -1 * 10*M_PI * cos(10*M_PI * hPosX[i]);
		rVelDotY[i] = -1 * 2.0 *hPosY[i];
		rDensityDot[i] = -1000.0 * (10*M_PI * cos(10*M_PI * hPosX[i]) + 2.0 *hPosY[i]);
	}
	
	for (int i = 0; i < hPN; i++) {
		if ((hPosX[i] < 1*hGrid.size) ||
			(hPosY[i] < 1*hGrid.size) ||
			(hPosX[i] > 1.0 - 1*hGrid.size) ||
			(hPosY[i] > 1.0 - 1*hGrid.size)) {
				hVelDotX[i] = rVelDotX[i];
				hVelDotY[i] = rVelDotY[i];
				hDensityDot[i] = rDensityDot[i];
		}
	}
	
	sumVelDotX = 0.0;
	maxVelDotX = 0.0;
	sumVelDotY = 0.0;
	maxVelDotY = 0.0;
	sumDensityDot = 0.0;
	maxDensityDot = 0.0;
	
	for (int i = 0; i < hPN; i++) {
		sumVelDotX += abs(rVelDotX[i] - hVelDotX[i]);
		sumVelDotY += abs(rVelDotY[i] - hVelDotY[i]);
		sumDensityDot += abs(rDensityDot[i] - hDensityDot[i]);
		
		if (abs(rVelDotX[i]) > maxVelDotX) 
			maxVelDotX = abs(rVelDotX[i]);
		if (abs(rVelDotY[i]) > maxVelDotY) 
			maxVelDotY = abs(rVelDotY[i]);
		if (abs(rDensityDot[i]) > maxDensityDot) 
			maxDensityDot = abs(rDensityDot[i]);		
	}
	
	printf("VelDotX errors %f%% \n", 100 * sumVelDotX/hPN/maxVelDotX);
	printf("VelDotY errors %f%% \n", 100 * sumVelDotY/hPN/maxVelDotY);
	printf("DensityDot errors %f%% \n", 100 * sumDensityDot/hPN/maxDensityDot);
	
    stream = fopen("debug.txt", "w");
    for (int i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e\n", hDensityDot[i],
                hVelDotX[i], hVelDotY[i]);
    fclose(stream);

    stream = fopen("reference.txt", "w");
    for (int i = 0; i < hPN; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e\n", rDensityDot[i],
                rVelDotX[i], rVelDotY[i]);
    fclose(stream);
    
	for (int i = 0; i < hPN; i++) {
		hVelX[i] = hVelDotX[i];
		hVelY[i] = hVelDotY[i];
		hDensity[i] = hDensityDot[i];
	}
	
	hRun.tsn = 1;
	outputCase();
	outputData(0);
	
	return 0;
}


int main() {
    /**
     * \brief armando2D v2.0
     *
     * An SPH code for non stationary fluid dynamics.
     * This is the reviewed and improved C version of Armando v1.0
     * developed at CERN in 2008
     *
     * \date Oct 20, 2010
     * \author Luca Massidda
     */

    initHost();

    initDamBreak();
    //initFree();

    //initDevice();
    //initCUDPP();
    
	//RKintegrateDevice();
	RKintegrateHost();
	//check();
	
    return 0;
}
