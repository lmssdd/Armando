#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/for_each.h>
#include <thrust/system_error.h>


#include <cutil_inline.h>
#include <cudpp.h>

#define MAXP 20000
#define MAXN 21
#define MAXG 200000
#define PI 3.14159f

#define THREADS 256
#define ParticlesInSet 12
#define SetsInBlock 20

struct grid {
    float oX, oY;
    float size;
    int nX, nY;
};

struct simulation {
    float minX, maxX;
    float minY, maxY;
    float dt;
    int tsn;
    int ssi;
    int nsi;
};

struct load {
    float minX, maxX;
    float minY, maxY;
    float gx;
    float gy;
    float w;
};


// Host Variables

thrust::host_vector<int> hMaterial;
thrust::host_vector<float> hPosX;
thrust::host_vector<float> hPosY;
thrust::host_vector<float> hVelX;
thrust::host_vector<float> hVelY;
thrust::host_vector<float> hDensity;
thrust::host_vector<float> hEnergy;
thrust::host_vector<float> hPressure;
thrust::host_vector<float> hVelDotX;
thrust::host_vector<float> hVelDotY;
thrust::host_vector<float> hDensityDot;
thrust::host_vector<float> hEnergyDot;
thrust::host_vector<int> hList;
thrust::host_vector<int> hHash;
thrust::host_vector<int> hIndex;
thrust::host_vector<int> hSetStart;
thrust::host_vector<int> hSetStop;

int hPN;
float hSmooth, hMass, hSound;
int hMatType[10];
float hMatProp[10][10];
struct simulation hRun;
struct grid hGrid;
struct load hLoad[10];

// Device Variables

thrust::device_vector<int> dMaterial;
thrust::device_vector<float> dPosX;
thrust::device_vector<float> dPosY;
thrust::device_vector<float> dVelX;
thrust::device_vector<float> dVelY;
thrust::device_vector<float> dDensity;
thrust::device_vector<float> dEnergy;
thrust::device_vector<float> dPressure;
thrust::device_vector<float> dVelDotX;
thrust::device_vector<float> dVelDotY;
thrust::device_vector<float> dDensityDot;
thrust::device_vector<float> dEnergyDot;
thrust::device_vector<int> dList;
thrust::device_vector<int> dHash;
thrust::device_vector<int> dIndex;
thrust::device_vector<int> dSetStart;
thrust::device_vector<int> dSetStop;

__device__ __constant__ int dPN;
__device__ __constant__ float dSmooth, dMass, dSound;
__device__ __constant__ int dMatType[10];
__device__ __constant__ float dMatProp[10][10];
__device__ __constant__ struct simulation dRun;
__device__ struct grid dGrid;
__device__ struct load dLoad[10];

thrust::device_vector<float> dPosX0;
thrust::device_vector<float> dPosY0;
thrust::device_vector<float> dVelX0;
thrust::device_vector<float> dVelY0;
thrust::device_vector<float> dDensity0;
thrust::device_vector<float> dEnergy0;

thrust::device_vector<int> dIntDummy;
thrust::device_vector<float> dFloatDummy;

int *pMaterial;
float *pPosX, *pPosY;
float *pVelX, *pVelY;
float *pDensity, *pEnergy, *pPressure;
float *pVelDotX, *pVelDotY;
float *pDensityDot, *pEnergyDot;
float *pPosX0, *pPosY0;
float *pVelX0, *pVelY0;
float *pDensity0, *pEnergy0;
int *pHash, *pIndex, *pSetStart, *pSetStop, *pList;
int *pIntDummy;
float *pFloatDummy;


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
	//alpha = 21.0f / (16.0f * M_PI * h * h * h);
	
    // for 2D
	alpha = 7.0f / (4.0f * PI * h * h);
	
    w = 0.0f;
    if (q < 2) {
        w = 1.0f - 0.5f*q;
        w *= w;
        w *= w;
        w *= 1.0f + 2.0f*q;
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
	//alpha = 21.0f / (16.0f * M_PI * h * h * h);
	
    // for 2D
	alpha = 7.0f / (4.0f * PI * h * h);
	
    dwdr = 0.0f;
    if (q < 2) {
        dwdr = 5.0f / 8.0f * q * powf((q - 2.0f), 3) ;
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
    alpha = 1.0f / (3.14f * h2);

    w = 0.0f;
    //dwdr = 0.0;

    if (q2 < 4.0f) {
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
    alpha = 1.0f / (3.14f * h2);

    w = 0.0f;
    dwdr = 0.0f;

    if (q2 < 4.0f) {
        w = alpha * expf(-q2);
        dwdr = w * (-2.0f * r / h2);
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

    p = (dMatProp[mat][1] - 1.0f) * rho * u;
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

    ph = (dMatProp[mat][0] * powf(dMatProp[mat][1], 2) * mu*(1.0f +mu))
         / powf((1.0f - (dMatProp[mat][3] -1.0f) * mu), 2);

    p = ph + dMatProp[mat][2] * dMatProp[mat][0]
        * (u - (0.5f * ph / dMatProp[mat][0] * (mu / (1.0f + mu))));

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

    p = dMatProp[mat][0] * powf(dMatProp[mat][1], 2) / 7.0f
        * (powf((rho / dMatProp[mat][0]), 7) - 1.0f);

    if (p < dMatProp[mat][2]) p = dMatProp[mat][2];

//    c = dMatProp[mat][1];

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
		iDensityDot = 0.0f;
		iVelDotX = 0.0f;
		iVelDotY = 0.0f;
		
		for (il = 0; il < MAXN; il++) {
			jp = dList[ip * MAXN + il];
			
			dx = dPosX[ip] - dPosX[jp];
			dy = dPosY[ip] - dPosY[jp];
			dz = 0.0f;
			dr = sqrtf(dx * dx + dy * dy + dz * dz);
			
			if (dr < 0.1f * dSmooth) dr = 100.0f * dSmooth;
			
			//dwdr = kernelDerivGauss(dr, dSmooth);
			dwdr = kernelDerivWendland(dr, dSmooth);
			
			dvr = 0.0f;
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
			f *= 2.0f * dSound / (dDensity[ip] + dDensity[jp]);
			
			iDensityDot += dMass * f * dwdr;
			
			// Calculate shock correction for momentum
			if (dvr < 0.0f) f = dvr;
			else f = 0.0f;
			
			f *= dSmooth / (dr * dr + 0.01f * dSmooth * dSmooth);
			f *= 2.0f * dSound / (dDensity[ip] + dDensity[jp]);
			f *= 0.03f;
			
			iVelDotX += dMass * f * dwdr * (dPosX[ip] - dPosX[jp]) / dr;
			iVelDotY += dMass * f * dwdr * (dPosY[ip] - dPosY[jp]) / dr;
		}
		
		dDensityDot[ip] += iDensityDot;
		dVelDotX[ip] += iVelDotX;
		dVelDotY[ip] += iVelDotY;

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
            dPressure[ip] = 0.0f;
        }
        
        
	}
}


struct updateHash {
    const float ox, oy, d;
    const int nx, ny;
    updateHash(float _ox, float _oy, float _d, int _nx, int _ny) : 
    ox(_ox), oy(_oy), d(_d), nx(_nx), ny(_ny) {}
    
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t) {
        int ix, iy, ic;
        float px, py;
        
        px = thrust::get<0>(t);
        py = thrust::get<1>(t);
        
        ix = (int) truncf((px - ox) / d);
        iy = (int) truncf((py - oy) / d);
		ic = ix + iy * nx;
		
        thrust::get<2>(t) = ic;
    }
};


struct updateForces {
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t) {
        
        int iMaterial;
        float iVelDotX, iVelDotY, iDensityDot, iEnergyDot;
        
        iVelDotX = 0.0f;
        iVelDotY = 0.0f;
        iDensityDot = 0.0f;
        iEnergyDot = 0.0f;

        iMaterial = thrust::get<0>(t);
        
        if (iMaterial > 0) iVelDotY = -9.81f;
        
        thrust::get<1>(t) = iVelDotX;
        thrust::get<2>(t) = iVelDotY;
        thrust::get<3>(t) = iDensityDot;
        thrust::get<4>(t) = iEnergyDot;
    }
};


struct updateBoundary {
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t) {
        int iMaterial;
        
        iMaterial = thrust::get<0>(t);
        
        if (iMaterial < 0) {
            thrust::get<3>(t) = thrust::get<1>(t);
            thrust::get<4>(t) = thrust::get<2>(t);
        }
    }
};


struct updatePressure {
    template <typename Tuple>
    __device__
    void operator()(Tuple t) {
        int iMaterial;
        float iDensity, iEnergy, iPressure;
        
        iMaterial = abs(thrust::get<0>(t));
        iDensity = thrust::get<1>(t);
        iEnergy = thrust::get<2>(t);
        
        switch (dMatType[iMaterial]) {
        case (1) : // IDEAL GAS EOS
            iPressure = pressureGas(iMaterial, iDensity, iEnergy);
            break;
        case (2) : // MIE-GRUNEISEN POLYNOMIAL EOS
            iPressure = pressurePoly(iMaterial, iDensity, iEnergy);
            break;
        case (3) : // MIE-GRUNEISEN SHOCK EOS
            iPressure = pressureShock(iMaterial, iDensity, iEnergy);
            break;
        case (4) : // TAIT EOS
            iPressure = pressureTait(iMaterial, iDensity, iEnergy);
            break;
        default :
            iPressure = 0.0f;
        }
        
        thrust::get<3>(t) = iPressure;
    }
};


struct updateEnergyDot {
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t) {
        float iDensity, iPressure, iDensityDot, iEnergyDot;
        
        iDensity = thrust::get<0>(t);
        iPressure = thrust::get<1>(t);
        iDensityDot = thrust::get<2>(t);
        
        iEnergyDot = (iPressure * iDensityDot) / (iDensity * iDensity);

        thrust::get<3>(t) += iEnergyDot;
    }
};


struct updateVector {
    const float alpha;
    updateVector(float _alpha) : alpha(_alpha) {}
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t) {
        float u0, ud, u1, u2;
        
        u0 = thrust::get<0>(t);
        ud = thrust::get<1>(t);
        u1 = thrust::get<2>(t);
        u2 = u0 + alpha * (u1 + dRun.dt * ud - u0);
        thrust::get<2>(t) = u2;
    }
};


__global__ void updateForcesDevice(const int* dMaterial,
		float* dVelDotX, float* dVelDotY, float* dDensityDot,
		float* dEnergyDot) {
	
    int ip;
    int iMaterial;
    float iVelDotX, iVelDotY, iDensityDot, iEnergyDot;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < dPN) {
        iVelDotX = 0.0f;
        iVelDotY = 0.0f;
        iDensityDot = 0.0f;
        iEnergyDot = 0.0f;

        iMaterial = dMaterial[ip];

        if (iMaterial > 0) iVelDotY = -9.81f;

        dVelDotX[ip] = iVelDotX;
        dVelDotY[ip] = iVelDotY;
        dDensityDot[ip] = iDensityDot;
        dEnergyDot[ip] = iEnergyDot;
    }
}


__global__ void updateLoadsDevice(const int* dMaterial,
		float* dPosX, float* dPosY, 
		float* dVelDotX, float* dVelDotY, float* dEnergyDot) {
	
    int ip, i;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if ((ip < dPN) && (dMaterial[ip] > 0)) {
        for (i = 0; i < 10; i++) {
            if ((dPosX[ip] > dLoad[i].minX) && 
                (dPosX[ip] < dLoad[i].maxX) && 
                (dPosY[ip] > dLoad[i].minY) && 
                (dPosY[ip] < dLoad[i].maxY)) {
                    dVelDotX[ip] += dLoad[i].gx;
                    dVelDotY[ip] += dLoad[i].gy;
                    dEnergyDot[ip] += dLoad[i].w;
            }
        }
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
				
				if ((dr < 2.0f * dSmooth) && (np < MAXN)) {
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


int initDevice() {
    
    try {
        dMaterial = hMaterial;
        dPosX = hPosX;
        dPosY = hPosY;
        dVelX = hVelX;
        dVelY = hVelY;
        dDensity = hDensity;
        dEnergy = hEnergy;
        dPressure = hPressure;
        
        dPosX0 = hPosX;
        dPosY0 = hPosY;
        dVelX0 = hVelX;
        dVelY0 = hVelY;
        dDensity0 = hDensity;
        dEnergy0 = hEnergy;
        
        dVelDotX = hVelDotX;
        dVelDotY = hVelDotY;
        dDensityDot = hDensityDot;
        dEnergyDot = hEnergyDot;
        dList = hList;
        dHash = hHash;
        dIndex = hIndex;
        dSetStart = hSetStart;
        dSetStop = hSetStop;
        
        dIntDummy.resize(hPN);
        dFloatDummy.resize(hPN);
    } catch(std::bad_alloc &e) {
        printf("Ran out of memory while copying to device\n");
        exit(-1);
    } catch(thrust::system_error &e) {
        printf("Error while copying vector on device: %s \n", e.what());
        exit(-1);
    }
    
    try {
        pMaterial = thrust::raw_pointer_cast(&dMaterial[0]);
        pPosX = thrust::raw_pointer_cast(&dPosX[0]);
        pPosY = thrust::raw_pointer_cast(&dPosY[0]);
        pVelX = thrust::raw_pointer_cast(&dVelX[0]);
        pVelY = thrust::raw_pointer_cast(&dVelY[0]);
        pDensity = thrust::raw_pointer_cast(&dDensity[0]);
        pEnergy = thrust::raw_pointer_cast(&dEnergy[0]);
        pPressure = thrust::raw_pointer_cast(&dPressure[0]);
        pVelDotX = thrust::raw_pointer_cast(&dVelDotX[0]);
        pVelDotY = thrust::raw_pointer_cast(&dVelDotY[0]);
        pDensityDot = thrust::raw_pointer_cast(&dDensityDot[0]);
        pEnergyDot = thrust::raw_pointer_cast(&dEnergyDot[0]);
        pPosX0 = thrust::raw_pointer_cast(&dPosX0[0]);
        pPosY0 = thrust::raw_pointer_cast(&dPosY0[0]);
        pVelX0 = thrust::raw_pointer_cast(&dVelX0[0]);
        pVelY0 = thrust::raw_pointer_cast(&dVelY0[0]);
        pDensity0 = thrust::raw_pointer_cast(&dDensity0[0]);
        pEnergy0 = thrust::raw_pointer_cast(&dEnergy0[0]);
        pHash = thrust::raw_pointer_cast(&dHash[0]);
        pIndex = thrust::raw_pointer_cast(&dIndex[0]);
        pSetStart = thrust::raw_pointer_cast(&dSetStart[0]);
        pSetStop = thrust::raw_pointer_cast(&dSetStop[0]);
        pList = thrust::raw_pointer_cast(&dList[0]);
        pIntDummy = thrust::raw_pointer_cast(&dIntDummy[0]);
        pFloatDummy = thrust::raw_pointer_cast(&dFloatDummy[0]);
    } catch(std::bad_alloc &e) {
        printf("Ran out of memory while castin pointers\n");
        exit(-1);
    } catch(thrust::system_error &e) {
        printf("Error while casting pointers on vectors: %s \n", e.what());
        exit(-1);
    }
    
    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.size = hGrid.size;
    
    cutilSafeCall( cudaMemcpyToSymbol("dPN", &hPN, sizeof(int)) );
    cutilSafeCall( cudaMemcpyToSymbol("dSmooth", &hSmooth, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMass", &hMass, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dSound", &hSound, sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dRun", &hRun, sizeof(struct simulation)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMatType", hMatType, 10 * sizeof(int)) );
    cutilSafeCall( cudaMemcpyToSymbol("dMatProp", hMatProp, 100 * sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol("dLoad", &hLoad, 10 * sizeof(struct load)) );

    return 0;
}


int copyHostToDevice() {

    thrust::copy(hMaterial.begin(), hMaterial.end(), dMaterial.begin());
    thrust::copy(hPosX.begin(), hPosX.end(), dPosX.begin());
    thrust::copy(hPosY.begin(), hPosY.end(), dPosY.begin());
    thrust::copy(hVelX.begin(), hVelX.end(), dVelX.begin());
    thrust::copy(hVelY.begin(), hVelY.end(), dVelY.begin());
    thrust::copy(hDensity.begin(), hDensity.end(), dDensity.begin());
    thrust::copy(hEnergy.begin(), hEnergy.end(), dEnergy.begin());
    thrust::copy(hPressure.begin(), hPressure.end(), dPressure.begin());
    thrust::copy(hVelDotX.begin(), hVelDotX.end(), dVelDotX.begin());
    thrust::copy(hVelDotY.begin(), hVelDotY.end(), dVelDotY.begin());
    thrust::copy(hDensityDot.begin(), hDensityDot.end(), dDensityDot.begin());
    thrust::copy(hEnergyDot.begin(), hEnergyDot.end(), dEnergyDot.begin());
    
    thrust::copy(hList.begin(), hList.end(), dList.begin());
    thrust::copy(hHash.begin(), hHash.end(), dHash.begin());
    thrust::copy(hIndex.begin(), hIndex.end(), dIndex.begin());
                              
    thrust::copy(hSetStart.begin(), hSetStart.end(), dSetStart.begin());
    thrust::copy(hSetStop.begin(), hSetStop.end(), dSetStop.begin());
	
    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.size = hGrid.size;

    return 0;
}


int copyDeviceToHost() {

    thrust::copy(dMaterial.begin(), dMaterial.end(), hMaterial.begin());
    thrust::copy(dPosX.begin(), dPosX.end(), hPosX.begin());
    thrust::copy(dPosY.begin(), dPosY.end(), hPosY.begin());
    thrust::copy(dVelX.begin(), dVelX.end(), hVelX.begin());
    thrust::copy(dVelY.begin(), dVelY.end(), hVelY.begin());
    thrust::copy(dDensity.begin(), dDensity.end(), hDensity.begin());
    thrust::copy(dEnergy.begin(), dEnergy.end(), hEnergy.begin());
    thrust::copy(dPressure.begin(), dPressure.end(), hPressure.begin());
    thrust::copy(dVelDotX.begin(), dVelDotX.end(), hVelDotX.begin());
    thrust::copy(dVelDotY.begin(), dVelDotY.end(), hVelDotY.begin());
    thrust::copy(dDensityDot.begin(), dDensityDot.end(), hDensityDot.begin());
    thrust::copy(dEnergyDot.begin(), dEnergyDot.end(), hEnergyDot.begin());
    
    thrust::copy(dList.begin(), dList.end(), hList.begin());
    thrust::copy(dHash.begin(), dHash.end(), hHash.begin());
    thrust::copy(dIndex.begin(), dIndex.end(), hIndex.begin());
                              
    thrust::copy(dSetStart.begin(), dSetStart.end(), hSetStart.begin());
    thrust::copy(dSetStop.begin(), dSetStop.end(), hSetStop.begin());
	
    return 0;
}

int backupData() {
    
    thrust::copy(dPosX.begin(), dPosX.end(), dPosX0.begin());
    thrust::copy(dPosY.begin(), dPosY.end(), dPosY0.begin());
    thrust::copy(dVelX.begin(), dVelX.end(), dVelX0.begin());
    thrust::copy(dVelY.begin(), dVelY.end(), dVelY0.begin());
    thrust::copy(dDensity.begin(), dDensity.end(), dDensity0.begin());
    thrust::copy(dEnergy.begin(), dEnergy.end(), dEnergy0.begin());
    
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
                hMatProp[m][p] = 0.0f;

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
            if (fv > 0.0f)
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
        fprintf(stream, "%+e\n", 0.0f);

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
        fprintf(stream, "%+e\n", 0.0f);

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
        fprintf(stream, "%+e %+e %+e \n", hPosX[i], hPosY[i], 0.0f);

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
        fprintf(stream, "%+e %+e %+e \n", hVelX[i], hVelY[i], 0.0f);
    
    fclose(stream);

    return 0;
}


void initDamBreak() {

    int i, j, m, pi, k;
    double rho, c0, pmin;
    double dr;
    
    k = 1;
    
    m = 1;
    rho = 1000.f;
    c0 = 50.f;
    pmin = -1.e12;
    
    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.02f / k; // x4
    pi = 0;

    for (j = 0; j <= 50 * k ; j++) {
        for (i = 0; i <= 100 * k; i++) {
            hPosX[pi] = i * dr + 0.5f * dr;
            hPosY[pi] = j * dr + 0.5f * dr;

            hVelX[pi] = 0.0f;
            hVelY[pi] = 0.0f;
            hMaterial[pi] = m;
            hDensity[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0f;
            hPressure[pi] = 0.0f;
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
        for (i = -3; i <= 269 * k + 2; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0f;
            hVelY[pi] = 0.0f;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0f;
            hPressure[pi] = 0.0f;
            pi++;
        }
    }

    for (j = -0; j <= 80 * k; j++) {
        for (i = -3; i <= -1; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0f;
            hVelY[pi] = 0.0f;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0f;
            hPressure[pi] = 0.0f;
            pi++;
        }
    }

    for (j = -0; j <= 80 * k; j++) {
        for (i = 269 * k; i <= 269 * k +2; i++) {
            hPosX[pi] = i * dr;
            hPosY[pi] = j * dr;

            hVelX[pi] = 0.0f;
            hVelY[pi] = 0.0f;
            hMaterial[pi] = -m;
            hDensity[pi] = rho; // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0f;
            hPressure[pi] = 0.0f;
            pi++;
        }
    }

    hPN = pi;
    hSmooth = 1.2f * dr;
    hMass = rho * dr * dr;
    hSound = c0;

    hRun.minX = -1.0f;
    hRun.maxX =  6.0f;
    hRun.minY = -1.0f;
    hRun.maxY =  4.0f;

    hRun.dt = 4.0e-4 / k; //1.0e-3;
    hRun.tsn = 10000 * k; //1000;
    hRun.ssi = 200 * k;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.size = 2.0f * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    
    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].gy = -9.81f;
    
    printf("Dam break in a box \n");
    printf("Particles: %i \n", hPN);
}

void initPump() {

    int i, j, m, k;
    double rho, c0, pmin;
    double dr;
    
    k = 1;
    
    m = 1;
    rho = 1000.f;
    c0 = 50.f;
    pmin = -1.e12;
    
    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.02f / k;

    for (j = 0; j < 20 * k ; j++) {
        for (i = 0; i < 40 * k; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(m);
            hDensity.push_back(rho); //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }
    
    for (j = 20 * k; j < 50 * k ; j++) {
        for (i = 0; i < 20 * k -1; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(m);
            hDensity.push_back(rho); //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }
    
    for (j = 20 * k; j < 50 * k ; j++) {
        for (i = 20 * k +1; i < 40 * k; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(m);
            hDensity.push_back(rho); //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }
    
    for (j = -3; j < 0; j++) {
        for (i = -3; i < 40 * k + 3; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(-m);
            hDensity.push_back(rho); // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }

    for (j = 0; j < 80 * k; j++) {
        for (i = -3; i < 0; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(-m);
            hDensity.push_back(rho); // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }

    for (j = 0; j < 80 * k; j++) {
        for (i = 40 * k; i < 40 * k +3; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(-m);
            hDensity.push_back(rho); // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }
    
    for (j = 20 * k +1; j < 60 * k; j++) {
        for (i = 20 * k -1; i < 20 * k +1; i++) {
            hPosX.push_back(i * dr);
            hPosY.push_back(j * dr);

            hVelX.push_back(0.0f);
            hVelY.push_back(0.0f);
            hMaterial.push_back(-m);
            hDensity.push_back(rho); // + (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy.push_back(0.0f);
            hPressure.push_back(0.0f);
        }
    }

    hPN = hMaterial.size();
    if ((hPosX.size() != hPN) ||
	    (hPosY.size() != hPN) ||
	    (hVelX.size() != hPN) ||
	    (hVelY.size() != hPN) ||
	    (hDensity.size() != hPN) ||
	    (hEnergy.size() != hPN) ||
	    (hPressure.size() != hPN)) {
			printf("Wrong vector size\n");
			exit(-1);
		}
	
    hVelDotX.resize(hPN);
    hVelDotY.resize(hPN);
    hDensityDot.resize(hPN);
    hEnergyDot.resize(hPN);
    hList.resize(hPN*MAXN);
    hHash.resize(hPN);
    hIndex.resize(hPN);
    
    hSmooth = 1.2f * dr;
    hMass = rho * dr * dr;
    hSound = c0;

    hRun.minX = -1.0f;
    hRun.maxX =  4.0f;
    hRun.minY = -1.0f;
    hRun.maxY =  4.0f;

    hRun.dt = 4.0e-4 / k; //1.0e-3;
    hRun.tsn = 10000 * k; //1000;
    hRun.ssi = 200 * k;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.size = 2.0f * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    
    hSetStart.resize(hGrid.nX * hGrid.nY);
    hSetStop.resize(hGrid.nX * hGrid.nY);
    
    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].gy = -9.81f;
    
    hLoad[1].minX = dr * (20 * k);
    hLoad[1].maxX = dr * (40 * k);
    hLoad[1].minY = dr * (20 * k);
    hLoad[1].maxY = dr * (30 * k);
    hLoad[1].gy = 3.0f*9.81f;
    
    printf("Pump in a box \n");
    printf("Particles: %i \n", hPN);
}


void initFree() {

    int i, j, m, pi;
    double rho, c0, pmin;
    double dr;

    m = 1;
    rho = 1000.f;
    c0 = 50.f;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.01f; // x4
    pi = 0;

    for (j = 0; j < 100; j++) {
        for (i = 0; i < 100; i++) {
            hPosX[pi] = i * dr + 0.0 * dr;
            hPosY[pi] = j * dr + 0.0 * dr;

            hVelX[pi] = 0.0f;
            hVelY[pi] = 0.0f;
            hMaterial[pi] = m;
            hDensity[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
            hEnergy[pi] = 0.0f;
            hPressure[pi] = 1.0f;
            pi++;
        }
    }
    
    hPN = pi;
    hSmooth = 1.2f * dr;
    hMass = rho * dr * dr;
    hSound = c0;

    hRun.minX = -0.5f;
    hRun.maxX =  1.5f;
    hRun.minY = -0.5f;
    hRun.maxY =  1.5f;

    hRun.dt = 0.5e-2; //1.0e-3;
    hRun.tsn = 3; //1000;
    hRun.ssi = 1;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.size = 2.0f * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].gy = -9.81;
    
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


int sortArraysHost(void) {
	int *pHash, *pIndex;
	int *pMaterial;
	float *pPosX, *pPosY;
	float *pVelX, *pVelY;
	float *pDensity, *pEnergy, *pPressure;

    pHash = thrust::raw_pointer_cast(&hHash[0]);
    pIndex = thrust::raw_pointer_cast(&hIndex[0]);
    pMaterial = thrust::raw_pointer_cast(&hMaterial[0]);
    pPosX = thrust::raw_pointer_cast(&hPosX[0]);
    pPosY = thrust::raw_pointer_cast(&hPosY[0]);
    pVelX = thrust::raw_pointer_cast(&hVelX[0]);
    pVelY = thrust::raw_pointer_cast(&hVelY[0]);
    pDensity = thrust::raw_pointer_cast(&hDensity[0]);
    pEnergy = thrust::raw_pointer_cast(&hEnergy[0]);
    pPressure = thrust::raw_pointer_cast(&hPressure[0]);
    
    // Particles are re ordered
    
    iSort(pHash, pIndex, hPN);
    iSort(pMaterial, pIndex, hPN);
    fSort(pPosX, pIndex, hPN);
    fSort(pPosY, pIndex, hPN);
    fSort(pVelX, pIndex, hPN);
    fSort(pVelY, pIndex, hPN);
    fSort(pDensity, pIndex, hPN);
    fSort(pEnergy, pIndex, hPN);
    fSort(pPressure, pIndex, hPN);
    
    return 0;
}


int sortArraysDevice(void) {
	int blocks, threads;
    
	threads = THREADS;
	blocks = (hPN + threads - 1) / threads;

    // Particles are re ordered
    
	kerSortInt <<< blocks, threads >>>
	(pIntDummy, pMaterial, pIndex);
    cutilSafeCall( cudaMemcpy(pMaterial, pIntDummy, 
		(MAXP * sizeof(int)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pPosX, pIndex);
    cutilSafeCall( cudaMemcpy(pPosX, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pPosY, pIndex);
    cutilSafeCall( cudaMemcpy(pPosY, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pVelX, pIndex);
    cutilSafeCall( cudaMemcpy(pVelX, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pVelY, pIndex);
    cutilSafeCall( cudaMemcpy(pVelY, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pDensity, pIndex);
    cutilSafeCall( cudaMemcpy(pDensity, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pEnergy, pIndex);
    cutilSafeCall( cudaMemcpy(pEnergy, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
	kerSortFloat <<< blocks, threads >>>
	(pFloatDummy, pPressure, pIndex);
    cutilSafeCall( cudaMemcpy(pPressure, pFloatDummy, 
		(MAXP * sizeof(float)), cudaMemcpyDeviceToDevice) );
	
    return 0;
}


int updateGrid() {
	
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


int updateList(void) {
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

		
	/*
	for (ip = 0; ip < hPN; ip++) {
		for (il = 0; il < MAXN; il++) {
			printf("%d ", hList[ip *MAXN + il]);
		}
		printf("\n");
	}
	printf("\n");
*/
	return 0;
}


int neighbourListThrust() {
	int blocks, threads;
	float *pPosX, *pPosY;
	int *pSetStart, *pSetStop, *pList;
    
    pPosX = thrust::raw_pointer_cast(&dPosX[0]);
    pPosY = thrust::raw_pointer_cast(&dPosY[0]);
    pSetStart = thrust::raw_pointer_cast(&dSetStart[0]);
    pSetStop = thrust::raw_pointer_cast(&dSetStop[0]);
    pList = thrust::raw_pointer_cast(&dList[0]);
    
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
    // transform the points to their bucket indices
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
                     dPosX.begin(), dPosY.begin(), dHash.begin())), 
                     thrust::make_zip_iterator(thrust::make_tuple(
                     dPosX.end(), dPosY.end(), dHash.end())), 
                     updateHash(dGrid.oX, dGrid.oY, dGrid.size, 
                     dGrid.nX, dGrid.nY));

    thrust::sort_by_key(dHash.begin(), dHash.end(), 
                        thrust::make_zip_iterator(thrust::make_tuple(
                        dIndex.begin(), dMaterial.begin(), 
                        dPosX.begin(), dPosY.begin(), 
                        dVelX.begin(), dVelY.begin(), 
                        dDensity.begin(), dEnergy.begin(), 
                        dPressure.begin())));
    
    thrust::fill(dSetStart.begin(), dSetStart.end(), 0);
    thrust::fill(dSetStop.begin(), dSetStop.end(), 0);
    
    thrust::counting_iterator<unsigned int> search_begin(0);
    thrust::lower_bound(dHash.begin(),
                        dHash.end(),
                        search_begin,
                        search_begin + dGrid.nX*dGrid.nY,
                        dSetStart.begin());

    // find the end of each bucket's list of points
    thrust::upper_bound(dHash.begin(),
                        dHash.end(),
                        search_begin,
                        search_begin + dGrid.nX*dGrid.nY,
                        dSetStop.begin());
	
	updateListDevice <<< blocks, threads >>>
	(pList, pSetStart, pSetStop, dGrid, pPosX, pPosY);
	
	return 0;
}


int neighbourListDevice() {
	int blocks, threads;
    
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
	updateHashDevice <<< blocks, threads >>>
	(dGrid, pPosX, pPosY, pHash, pIndex);
	
    try {
        thrust::sort_by_key(dHash.begin(), dHash.end(), dIndex.begin());
        
        thrust::copy(dMaterial.begin(), dMaterial.end(), dIntDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dIntDummy.begin(), dMaterial.begin());
        
        thrust::copy(dPosX.begin(), dPosX.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dPosX.begin());
        
        thrust::copy(dPosY.begin(), dPosY.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dPosY.begin());
        
        thrust::copy(dVelX.begin(), dVelX.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dVelX.begin());
        
        thrust::copy(dVelY.begin(), dVelY.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dVelY.begin());
        
        thrust::copy(dDensity.begin(), dDensity.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dDensity.begin());
        
        thrust::copy(dEnergy.begin(), dEnergy.end(), dFloatDummy.begin());
        thrust::gather(dIndex.begin(), dIndex.end(), dFloatDummy.begin(), dEnergy.begin());
        
    } catch(std::bad_alloc &e) {
        printf("Ran out of memory while sorting\n");
        exit(-1);
    } catch(thrust::system_error &e) {
        printf("Error accessing vector element: %s \n", e.what());
        exit(-1);
    }
    
    thrust::fill(dSetStart.begin(), dSetStart.end(), 0);
    thrust::fill(dSetStop.begin(), dSetStop.end(), 0);
    
    updateSetsDevice <<< blocks, threads >>>
	(pSetStart, pSetStop, pHash);
    
	updateListDevice <<< blocks, threads >>>
	(pList, pSetStart, pSetStop, dGrid, pPosX, pPosY);
	
	return 0;
}

int RKstepDevice(float alpha) {
	int blocks, threads;
    
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
    thrust::fill(dVelDotX.begin(), dVelDotX.end(), 0.0f);
    thrust::fill(dVelDotY.begin(), dVelDotY.end(), 0.0f);
    thrust::fill(dDensityDot.begin(), dDensityDot.end(), 0.0f);
    thrust::fill(dEnergyDot.begin(), dEnergyDot.end(), 0.0f);
    
	// External loads
	updateLoadsDevice <<< blocks, threads >>>
	(pMaterial, pPosX, pPosY, pVelDotX, pVelDotY, pEnergyDot);
	
	// Calculate particle interactions
	balanceMassMomentumDevice <<< blocks, threads >>>
	(pList, pPosX, pPosY, pVelX, pVelY, pDensity, pPressure, 
	pDensityDot, pVelDotX, pVelDotY);
	
	balanceEnergyDevice <<< blocks, threads >>>
	(pPressure, pDensity, pDensityDot, pEnergyDot);
	
	// Update particles
	updateParticlesDevice  <<< blocks, threads >>> 
	(pMaterial, pVelDotX, pVelDotY, pDensityDot, pEnergyDot, alpha,
	pPosX0, pPosY0, pVelX0, pVelY0, pDensity0, pEnergy0,
	pPosX, pPosY, pVelX, pVelY, pDensity, pEnergy, pPressure);
	
	return 0;
}


int RKstepThrust(float alpha) {
	int blocks, threads;
	float *pPosX, *pPosY;
	float *pVelX, *pVelY;
	float *pDensity, *pPressure;
	float *pVelDotX, *pVelDotY;
	float *pDensityDot;
	int *pList;
    
    pPosX = thrust::raw_pointer_cast(&dPosX[0]);
    pPosY = thrust::raw_pointer_cast(&dPosY[0]);
    pVelX = thrust::raw_pointer_cast(&dVelX[0]);
    pVelY = thrust::raw_pointer_cast(&dVelY[0]);
    pDensity = thrust::raw_pointer_cast(&dDensity[0]);
    pPressure = thrust::raw_pointer_cast(&dPressure[0]);
    pVelDotX = thrust::raw_pointer_cast(&dVelDotX[0]);
    pVelDotY = thrust::raw_pointer_cast(&dVelDotY[0]);
    pDensityDot = thrust::raw_pointer_cast(&dDensityDot[0]);
    pList = thrust::raw_pointer_cast(&dList[0]);
    
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
	// External forces
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.begin(), dVelDotX.begin(), dVelDotY.begin(), 
            dDensityDot.begin(), dEnergyDot.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.end(), dVelDotX.end(), dVelDotY.end(), 
            dDensityDot.end(), dEnergyDot.end())),
            updateForces());
	
	// Calculate particle interactions
    
	balanceMassMomentumDevice <<< blocks, threads >>>
	(pList, pPosX, pPosY, pVelX, pVelY, pDensity, pPressure, 
	pDensityDot, pVelDotX, pVelDotY);
	
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dDensity.begin(), dPressure.begin(), dDensityDot.begin(), dEnergyDot.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dDensity.end(), dPressure.end(), dDensityDot.end(), dEnergyDot.end())), 
            updateEnergyDot());
    
    
	// Update particle position
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dPosX0.begin(), dVelX.begin(), dPosX.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dPosX0.end(), dVelX.end(), dPosX.end())), 
            updateVector(alpha));
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dPosY0.begin(), dVelY.begin(), dPosY.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dPosY0.end(), dVelY.end(), dPosY.end())), 
            updateVector(alpha));
    
    
	// Update particle velocity
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dVelX0.begin(), dVelDotX.begin(), dVelX.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dVelX0.end(), dVelDotX.end(), dVelX.end())), 
            updateVector(alpha));
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dVelY0.begin(), dVelDotY.begin(), dVelY.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dVelY0.end(), dVelDotY.end(), dVelY.end())), 
            updateVector(alpha));
    
    
	// Update particle density
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dDensity0.begin(), dDensityDot.begin(), dDensity.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dDensity0.end(), dDensityDot.end(), dDensity.end())), 
            updateVector(alpha));
    
    
    // Update particle energy
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dEnergy0.begin(), dEnergyDot.begin(), dEnergy.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dEnergy0.end(), dEnergyDot.end(), dEnergy.end())), 
            updateVector(alpha));
    
    
    // Update particles of boundary
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.begin(), dVelX0.begin(), dVelY0.begin(), dVelX.begin(), dVelY.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.end(), dVelX0.begin(), dVelY0.begin(), dVelX.end(), dVelY.end())), 
            updateBoundary());
    
    
    // Update particle pressure
    
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.begin(), dDensity.begin(), dEnergy.begin(), dPressure.begin())), 
            thrust::make_zip_iterator(thrust::make_tuple(
            dMaterial.end(), dDensity.end(), dEnergy.end(), dPressure.end())), 
            updatePressure());
    
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
					
					if ((dr < 2.0f * hSmooth) && (np < MAXN)) {
						hList[ip * MAXN + np] = jp;
						np++;
					}
				}
			}
		}
		
		
	}
	
	return 0;
}


int neighbourListHost() {
	int *pIndex;
    
    pIndex = thrust::raw_pointer_cast(&hIndex[0]);
    
	updateHashHost();
	
	qsort(pIndex, hPN, sizeof(int), indexCompare);
	
	sortArraysHost();
	
	updateSetsHost();
	
	updateListHost();
	
	return 0;
}

int neighbourList() {
	int *phIndex, *pdIndex;
	int blocks, threads;
	float *pPosX, *pPosY;
	int *pHash, *pSetStart, *pSetStop, *pList;
    
    phIndex = thrust::raw_pointer_cast(&hIndex[0]);
    pdIndex = thrust::raw_pointer_cast(&dIndex[0]);
    pPosX = thrust::raw_pointer_cast(&dPosX[0]);
    pPosY = thrust::raw_pointer_cast(&dPosY[0]);
    pHash = thrust::raw_pointer_cast(&dHash[0]);
    pSetStart = thrust::raw_pointer_cast(&dSetStart[0]);
    pSetStop = thrust::raw_pointer_cast(&dSetStop[0]);
    pList = thrust::raw_pointer_cast(&dList[0]);
    
    
	blocks = (hPN + THREADS - 1) / THREADS;
	threads = THREADS;
	
	//updateHashHost();
    updateHashDevice <<< blocks, threads >>>
	(dGrid, pPosX, pPosY, pHash, pdIndex);
	
    copyDeviceToHost();
    
	qsort(phIndex, hPN, sizeof(int), indexCompare);
    //thrust::sort_by_key(dHash.begin(), dHash.end(), dIndex.begin());
    
    copyHostToDevice();
    
	//sortArraysHost();
	sortArraysDevice();
    
	//updateSetsHost();
    thrust::fill(dSetStart.begin(), dSetStart.end(), 0);
    thrust::fill(dSetStop.begin(), dSetStop.end(), 0);
	updateSetsDevice <<< blocks, threads >>>
	(pSetStart, pSetStop, pHash);
	
	//updateListHost();
	updateListDevice <<< blocks, threads >>>
	(pList, pSetStart, pSetStop, dGrid, pPosX, pPosY);
	
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

    // TIME CYCLE
    for (ts = 0; ts <= hRun.tsn; ts++) {
		
        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
			copyDeviceToHost();
            printData();
            outputVTK(ts / hRun.ssi);
        }
        
		// Calculate neighbour list
        neighbourListDevice();
		
		// Save initial condition
        backupData();
		
        // Step 1
		RKstepDevice(1.0f);
		
        // Step 2
		RKstepDevice(1.0f / 4.0f);
		
        // Step 3
		RKstepDevice(2.0f / 3.0f);
	}
	
	return 0;
}

int check() {
    float *rVelDotX, *rVelDotY, *rDensityDot;
    float sumVelDotX, maxVelDotX;
    float sumVelDotY, maxVelDotY;
    float sumDensityDot, maxDensityDot;
    FILE *stream;
    
	float *pPosX, *pPosY;
	float *pVelX, *pVelY;
	float *pDensity, *pPressure;
	float *pVelDotX, *pVelDotY;
	float *pDensityDot;
	int *pList;
    
    pPosX = thrust::raw_pointer_cast(&dPosX[0]);
    pPosY = thrust::raw_pointer_cast(&dPosY[0]);
    pVelX = thrust::raw_pointer_cast(&dVelX[0]);
    pVelY = thrust::raw_pointer_cast(&dVelY[0]);
    pDensity = thrust::raw_pointer_cast(&dDensity[0]);
    pPressure = thrust::raw_pointer_cast(&dPressure[0]);
    pVelDotX = thrust::raw_pointer_cast(&dVelDotX[0]);
    pVelDotY = thrust::raw_pointer_cast(&dVelDotY[0]);
    pDensityDot = thrust::raw_pointer_cast(&dDensityDot[0]);
    pList = thrust::raw_pointer_cast(&dList[0]);
    
    rVelDotX = (float *) malloc(MAXP * sizeof(float));
    rVelDotY = (float *) malloc(MAXP * sizeof(float));
    rDensityDot = (float *) malloc(MAXP * sizeof(float));
	
	for (int i = 0; i < hPN; i++) {
		hPressure[i] = 1000.0f * (sinf(10*PI * hPosX[i]) + hPosY[i]*hPosY[i]);
		hVelX[i] = sinf(10*PI * hPosX[i]);
		hVelY[i] = hPosY[i]*hPosY[i];
	}
	
	int blocks = (hPN + THREADS - 1) / THREADS;
	copyHostToDevice();
	
	neighbourListDevice();
	
	// Calculate particle interactions
	balanceMassMomentumDevice <<< blocks, THREADS >>>
	(pList, pPosX, pPosY, pVelX, pVelY, pDensity, pPressure, 
	pDensityDot, pVelDotX, pVelDotY);
	
	copyDeviceToHost();
	
	for (int i = 0; i < hPN; i++) {
		rVelDotX[i] = -1 * 10*PI * cosf(10*PI * hPosX[i]);
		rVelDotY[i] = -1.0f * 2.0f *hPosY[i];
		rDensityDot[i] = -1000.0f * (10*PI * cosf(10*PI * hPosX[i]) + 2.0 *hPosY[i]);
	}
	
	for (int i = 0; i < hPN; i++) {
		if ((hPosX[i] < 1*hGrid.size) ||
			(hPosY[i] < 1*hGrid.size) ||
			(hPosX[i] > 1.0f - 1*hGrid.size) ||
			(hPosY[i] > 1.0f - 1*hGrid.size)) {
				hVelDotX[i] = rVelDotX[i];
				hVelDotY[i] = rVelDotY[i];
				hDensityDot[i] = rDensityDot[i];
		}
	}
	
	sumVelDotX = 0.0f;
	maxVelDotX = 0.0f;
	sumVelDotY = 0.0f;
	maxVelDotY = 0.0f;
	sumDensityDot = 0.0f;
	maxDensityDot = 0.0f;
	
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
     * \date May 2, 2012
     * \author Luca Massidda
     */
    
    //initHost();
    
    initPump();
    //initDamBreak();
    //initFree();
    
    initDevice();
    
    hHash.resize(hPN);
    dHash.resize(hPN);
    hIndex.resize(hPN);
    dIndex.resize(hPN);
    
	RKintegrateDevice();
	//check();

    return 0;
}
