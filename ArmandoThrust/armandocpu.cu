#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/count.h>

#define PI 3.14159f
#define MAXP 8000
#define MAXN 96
#define MAXG 1280000

#define THREADS 256

struct pair {
	int key;
	int value;
};

struct  grid {
    float oX, oY, oZ;
    float size;
    int nX, nY, nZ;
};

struct  simulation {
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
    float dt;
    int tsn;
    int ssi;
    int nsi;
};

struct load {
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
    float gx;
    float gy;
    float gz;
    float w;
};

struct fix {
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
    float velX, velY, velZ;
};

struct outlet {
    float oX, oY, oZ;
    float nX, nY, nZ;
};

struct inlet {
    int Material;
    float Mass, Smooth;
    float oX, oY, oZ;
    float uX, uY, uZ;
    float vX, vY, vZ;
    float nX, nY, nZ;
    int nu, nv;
    float Velocity;
    float Density, Energy;
    float Distance;
};


struct model {
    int pn;
    int* Material;
    float* Mass;
    float* Smooth;
    float* PosX;
    float* PosY;
    float* PosZ;
    float* VelX;
    float* VelY;
    float* VelZ;
    float* Density;
    float* Energy;
    float* Pressure;
    float* Sound;
    float* VelDotX;
    float* VelDotY;
    float* VelDotZ;
    float* DensityDot;
    float* EnergyDot;
    float* PosX0;
    float* PosY0;
    float* PosZ0;
    float* VelX0;
    float* VelY0;
    float* VelZ0;
    float* Density0;
    float* Energy0;
    int* List;
    int* Hash;
    int* Index;
    int* SetStart;
    int* SetStop;
    int* IntDummy;
    float* FloatDummy;
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
struct load hLoad[10];
struct fix hFix[10];
struct outlet hOut[10];
struct inlet hIn[10];

float *hPosX0;
float *hPosY0;
float *hVelX0;
float *hVelY0;
float *hDensity0;
float *hEnergy0;

// Device Variables
__device__ __constant__ int dMatType[10];
__device__ __constant__ float dMatProp[10][10];
__device__ __constant__ struct simulation dRun;
__device__ struct grid dGrid;
__device__ __constant__ struct load dLoad[10];
__device__ __constant__ struct fix dFix[10];
__device__ __constant__ struct outlet dOut[10];
__device__ struct inlet dIn[10];


__host__ __device__ float kernelWendland(float r, float h) {

    float q, alpha, w;
    /**
     * \brief Wendland kernel
     *
     * \date Feb 8, 2011
     * \author Luca Massidda
     */

    q = r / h;

    // for 3D
    alpha = 15.0f / (16.0f * PI * h * h * h);

    // for 2D
    //alpha = 7.0f / (4.0f * PI * h * h);

    w = 0.0f;
    if (q < 2) {
        w = powf((1.0f - 0.5f*q),4);
        w *= 1.0f + 2.0f*q;
        w *= alpha;
    }

    return w;
}


__host__ __device__ float kernelDerivWendland(float r, float h) {

    float q, alpha, dwdr;
    /**
     * \brief Wendland kernel derivative
     *
     * \date Feb 8, 2011
     * \author Luca Massidda
     */

    q = r / h;

    // for 3D
    alpha = 15.0f / (16.0f * PI * h * h * h);

    // for 2D
    //alpha = 7.0f / (4.0f * PI * h * h);

    dwdr = 0.0f;
    if (q < 2) {
        dwdr = 5.0f / 8.0f * q * powf((q - 2.0f), 3) ;
        dwdr *= alpha / h;
    }

    return dwdr;
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

void balanceMassMomentumHost(const int pn, const int* List,
		const int* Material, const float* Mass, const float* Smooth,
		const float* PosX, const float* PosY, const float* PosZ,
		const float* VelX, const float* VelY, const float* VelZ,
		const float* Density, const float* Pressure, const float* Sound,
		float* DensityDot, float* VelDotX, float* VelDotY, float* VelDotZ) {

    /**
     * \brief Interate particles
     *
     * \date Jan 6, 2011
     * \author Luca Massidda
     */

    int ip, il, jp;
    float iDensityDot;
    float iVelDotX, iVelDotY, iVelDotZ;
    float iSmooth, jMass;
    float dx, dy, dz, dr, dvr, dwdr, f, w, w0;

    for (ip = 0; ip < pn; ip++) {
        iDensityDot = 0.0f;
        iVelDotX = 0.0f;
        iVelDotY = 0.0f;
        iVelDotZ = 0.0f;
        iSmooth = Smooth[ip];
		
        for (il = 0; il < MAXN; il++) {
            jp = List[ip * MAXN + il];
            
            jMass = Mass[jp];

            dx = PosX[ip] - PosX[jp];
            dy = PosY[ip] - PosY[jp];
			dz = PosZ[ip] - PosZ[jp];
            dr = sqrtf(dx * dx + dy * dy + dz * dz);

            if (dr < (0.01f * iSmooth)) dr = 100.0f * iSmooth;

			w = kernelWendland(dr, iSmooth);
			w0 = kernelWendland(0.0f, iSmooth);
            dwdr = kernelDerivWendland(dr, iSmooth);

			dvr = 0.0f;
			dvr += (PosX[ip] - PosX[jp]) * (VelX[ip] - VelX[jp]);
			dvr += (PosY[ip] - PosY[jp]) * (VelY[ip] - VelY[jp]);
			dvr += (PosZ[ip] - PosZ[jp]) * (VelZ[ip] - VelZ[jp]);
			
			iDensityDot += jMass * dvr * dwdr / dr;
			
			// Calculate interparticle pressure action
			//f = -(Pressure[ip] + Pressure[jp])
			//	/ (Density[ip] * Density[jp]);
			f = -(Pressure[ip] / powf(Density[ip], 2) + Pressure[jp] / powf(Density[jp], 2));

			iVelDotX += jMass * f * dwdr * (PosX[ip] - PosX[jp]) / dr;
			iVelDotY += jMass * f * dwdr * (PosY[ip] - PosY[jp]) / dr;
			iVelDotZ += jMass * f * dwdr * (PosZ[ip] - PosZ[jp]) / dr;
	
			// Calculate shock correction for mass
			f = Density[ip] - Density[jp];
			f *= 2.0f * Sound[ip] / (Density[ip] + Density[jp]);

			iDensityDot += jMass * f * dwdr;

			// Calculate shock correction for momentum
			if (dvr < 0.0f) f = dvr;
            else f = 0.0f;

            f *= iSmooth / (dr * dr + 0.01f * iSmooth * iSmooth);
            f *= 2.0f * Sound[ip] / (Density[ip] + Density[jp]);
            f *= 0.03f;

            iVelDotX += jMass * f * dwdr * (PosX[ip] - PosX[jp]) / dr;
            iVelDotY += jMass * f * dwdr * (PosY[ip] - PosY[jp]) / dr;
            iVelDotZ += jMass * f * dwdr * (PosZ[ip] - PosZ[jp]) / dr;
			
			// Calculate boundary repulsion
            if (Material[ip] != Material[jp]) {
				f = 0.02f * w / w0 * Sound[ip] * Sound[jp] / iSmooth;
				iVelDotX += jMass  / (Mass[ip] + jMass) * f * (PosX[ip] - PosX[jp]) / dr;
				iVelDotY += jMass  / (Mass[ip] + jMass) * f * (PosY[ip] - PosY[jp]) / dr;
				iVelDotZ += jMass  / (Mass[ip] + jMass) * f * (PosZ[ip] - PosZ[jp]) / dr;
			}
        }

        DensityDot[ip] += iDensityDot;
        VelDotX[ip] += iVelDotX;
        VelDotY[ip] += iVelDotY;
        VelDotZ[ip] += iVelDotZ;
    }
}

__global__ void balanceMassMomentumDevice(const int pn, const int* List,
        const int* Material, const float* Mass, const float* Smooth,
        const float* PosX, const float* PosY, const float* PosZ,
        const float* VelX, const float* VelY, const float* VelZ,
        const float* Density, const float* Pressure, const float* Sound,
        float* DensityDot, float* VelDotX, float* VelDotY, float* VelDotZ) {

    /**
     * \brief Interate particles
     *
     * \date Jan 6, 2011
     * \author Luca Massidda
     */

    int ip, il, jp;
    float iDensityDot;
    float iVelDotX, iVelDotY, iVelDotZ;
    float iSmooth, jMass;
    volatile float dx, dy, dz, dr, dvr, dwdr, f, w, w0, q;

    ip = threadIdx.x + blockDim.x * blockIdx.x;
	
    if (ip < pn) {
        iDensityDot = 0.0f;
        iVelDotX = 0.0f;
        iVelDotY = 0.0f;
        iVelDotZ = 0.0f;
        iSmooth = Smooth[ip];

        for (il = 0; il < MAXN; il++) {
            jp = List[ip * MAXN + il];

            jMass = Mass[jp];

            dx = PosX[ip] - PosX[jp];
            dy = PosY[ip] - PosY[jp];
            dz = PosZ[ip] - PosZ[jp];
            dr = sqrtf(dx * dx + dy * dy + dz * dz);

            if (dr < (0.01f * iSmooth)) dr = 100.0f * iSmooth;
			
			w = kernelWendland(dr, iSmooth);
			dwdr = kernelDerivWendland(dr, iSmooth);
			
            if (Material[ip] == Material[jp]) {
				dvr = 0.0f;
				dvr += (PosX[ip] - PosX[jp]) * (VelX[ip] - VelX[jp]);
				dvr += (PosY[ip] - PosY[jp]) * (VelY[ip] - VelY[jp]);
				dvr += (PosZ[ip] - PosZ[jp]) * (VelZ[ip] - VelZ[jp]);
				
				iDensityDot += jMass * dvr * dwdr / dr;
				
				// Calculate interparticle pressure action
				f = -(Pressure[ip] / powf(Density[ip], 2) + Pressure[jp] / powf(Density[jp], 2));
				f *= jMass * dwdr;
				iVelDotX += f * (PosX[ip] - PosX[jp]) / dr;
				iVelDotY += f * (PosY[ip] - PosY[jp]) / dr;
				iVelDotZ += f * (PosZ[ip] - PosZ[jp]) / dr;
				
				// Calculate shock correction for mass
				f = Density[ip] - Density[jp];
				f *= 2.0f * Sound[ip] / (Density[ip] + Density[jp]);
				iDensityDot += jMass * f * dwdr;
				
				// Calculate shock correction for momentum
				if (dvr < 0.0f) f = dvr;
				else f = 0.0f;
				
				f *= iSmooth / (dr * dr + 0.01f * iSmooth * iSmooth);
				f *= 2.0f * Sound[ip] / (Density[ip] + Density[jp]);
				f *= 0.03f;
				f *= jMass * dwdr;
				
				iVelDotX += f * (PosX[ip] - PosX[jp]) / dr;
				iVelDotY += f * (PosY[ip] - PosY[jp]) / dr;
				iVelDotZ += f * (PosZ[ip] - PosZ[jp]) / dr;
			}
			
			// Calculate boundary repulsion
            if (Material[ip] != Material[jp]) {
				f = 0.25f * w * Mass[jp] / Density[jp] / Smooth[jp] * powf(Sound[jp], 2);
				iVelDotX += f * (PosX[ip] - PosX[jp]) / dr;
				iVelDotY += f * (PosY[ip] - PosY[jp]) / dr;
				iVelDotZ += f * (PosZ[ip] - PosZ[jp]) / dr;
			}
        }
		
        DensityDot[ip] += iDensityDot;
        VelDotX[ip] += iVelDotX;
        VelDotY[ip] += iVelDotY;
        VelDotZ[ip] += iVelDotZ;
    }
}

void balanceMassMomentumHostOld(void) {
	
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
			dwdr = kernelDerivWendland(dr, hSmooth);
			
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


void balanceEnergyHost(const int pn,
                                const float* Pressure, const float* Density,
                                const float* DensityDot, float* EnergyDot) {

    /**
     * \brief Interate particles
     *
     * \date Jan 9, 2011
     * \author Luca Massidda
     */

    int ip;
    float iPressure, iDensity, iDensityDot;
    float iEnergyDot;

    for (ip = 0; ip < pn; ip++) {
        iPressure = Pressure[ip];
        iDensity = Density[ip];
        iDensityDot = DensityDot[ip];

        iEnergyDot = (iPressure * iDensityDot) / (iDensity * iDensity);

        EnergyDot[ip] += iEnergyDot;
    }
}

__global__ void balanceEnergyDevice(const int pn,
                                    const float* Pressure, const float* Density,
                                    const float* DensityDot, float* EnergyDot) {

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

    if (ip < pn) {
        iPressure = Pressure[ip];
        iDensity = Density[ip];
        iDensityDot = DensityDot[ip];

        iEnergyDot = (iPressure * iDensityDot) / (iDensity * iDensity);

        EnergyDot[ip] += iEnergyDot;
    }
}

void balanceEnergyHostOld(void) {

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



__host__ __device__ float pressureGas(float* properties, float rho, float u) {
    /**
     * \brief Ideal gas Equation Of State
     *
     * p = (k -1) rho u
     * c = (k(k -1) u)^0.5
     *
     * k = properties[1]
     * pshift = properties[2]
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;

    p = (properties[1] - 1.0f) * rho * u;
    p += properties[2];

    return p;
}



__host__ __device__ float pressurePoly(float* properties, float rho, float u) {
    /**
     * \brief Mie-Gruneisen polynomial Equation Of State
     *
     * p = a1 mu + a2 mu^2 + a3 mu^3 + (b0 + b1 mu) rho0 u  in compression
     * p = t1 mu + t2 mu^2 + b0 rho0 u                      in tension
     *
     * rho0 = properties[0];
     * a1 = properties[1];
     * a2 = properties[2];
     * a3 = properties[3];
     * b0 = properties[4];
     * b1 = properties[5];
     * t1 = properties[6];
     * t2 = properties[7];
     * pmin = properties[8];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p;

    mu = (rho - properties[0]) / properties[0];

    if (mu < 0)
        p = (properties[6] * mu + properties[7] * mu*mu)
            + (properties[4] * properties[0] * u);
    else
        p = (properties[1] * mu + properties[2] * mu*mu
             + properties[3] * mu*mu*mu)
            + ((properties[4] + properties[5] * mu)
               * properties[0] * u);

    //if (p < properties[8]) p = properties[8];

    return p;
}

__host__ __device__ float pressureShock(float* properties, float rho, float u) {
    /**
     * \brief Mie-Gruneisen Shock Hugoniot Equation Of State
     *
     * mu = rho / rho0 -1
     * g = g * rho0 / rho
     * ph = (rho0 c0^2 mu (1 + mu)) / (1 - (s0 - 1) * mu)^2
     * uh = 1/2 ph/rho0 * (mu / (1 + mu))
     * p = ph + g * rho * (u - uh)
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * g0 = properties[2];
     * s0 = properties[3];
     * pmin = properties[4];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float mu;
    float p, ph;

    mu = (rho - properties[0]) / properties[0];

    ph = (properties[0] * powf(properties[1], 2) * mu*(1.0f +mu))
         / powf((1.0f - (properties[3] -1.0f) * mu), 2);

    p = ph + properties[2] * properties[0]
        * (u - (0.5f * ph / properties[0] * (mu / (1.0f + mu))));

    //if (p < properties[4]) p = properties[4];

    return p;
}


__host__ __device__ float pressureTait(float* properties, float rho, float u) {
    /**
     * \brief Tait Equation Of State
     *
     * p = rho0 * c0 * c0 / 7.0 * (powf((rho / rho0), 7) - 1.0);
     * c = c0;
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * pmin = properties[2];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float p;

    p = properties[0] * powf(properties[1], 2) / 7.0f
        * (powf((rho / properties[0]), 7) - 1.0f);

    //if (p < properties[2]) p = properties[2];

    return p;
}


__host__ __device__ float soundGas(float* properties ,float rho, float u) {
    /**
     * \brief Ideal gas Equation Of State
     *
     * p = (k -1) rho u
     * c = (k(k -1) u)^0.5
     *
     * k = properties[1]
     * pshift = properties[2]
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float c;

    c = sqrtf(properties[1] * (properties[1] - 1.0f) * u);

    return c;
}



__host__ __device__ float soundPoly(float* properties , float rho, float u) {
    /**
     * \brief Mie-Gruneisen polynomial Equation Of State
     *
     * p = a1 mu + a2 mu^2 + a3 mu^3 + (b0 + b1 mu) rho0 u  in compression
     * p = t1 mu + t2 mu^2 + b0 rho0 u                      in tension
     *
     * rho0 = properties[0];
     * a1 = properties[1];
     * a2 = properties[2];
     * a3 = properties[3];
     * b0 = properties[4];
     * b1 = properties[5];
     * t1 = properties[6];
     * t2 = properties[7];
     * pmin = properties[8];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float c;

    c = sqrtf(properties[1] / rho);

    return c;
}

__host__ __device__ float soundShock(float* properties, float rho, float u) {
    /**
     * \brief Mie-Gruneisen Shock Hugoniot Equation Of State
     *
     * mu = rho / rho0 -1
     * g = g * rho0 / rho
     * ph = (rho0 c0^2 mu (1 + mu)) / (1 - (s0 - 1) * mu)^2
     * uh = 1/2 ph/rho0 * (mu / (1 + mu))
     * p = ph + g * rho * (u - uh)
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * g0 = properties[2];
     * s0 = properties[3];
     * pmin = properties[4];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float c;

    c = properties[1];

    return c;
}


__host__ __device__ float soundTait(float* properties, float rho, float u) {
    /**
     * \brief Tait Equation Of State
     *
     * p = rho0 * c0 * c0 / 7.0 * (powf((rho / rho0), 7) - 1.0);
     * c = c0;
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * pmin = properties[2];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float c;

    c = properties[1];

    return c;
}


__host__ __device__ float densityPoly(float* properties , float rho) {
    /**
     * \brief Mie-Gruneisen polynomial Equation Of State
     *
     * p = a1 mu + a2 mu^2 + a3 mu^3 + (b0 + b1 mu) rho0 u  in compression
     * p = t1 mu + t2 mu^2 + b0 rho0 u                      in tension
     *
     * rho0 = properties[0];
     * a1 = properties[1];
     * a2 = properties[2];
     * a3 = properties[3];
     * b0 = properties[4];
     * b1 = properties[5];
     * t1 = properties[6];
     * t2 = properties[7];
     * pmin = properties[8];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float rho0;
    
    rho0 = properties[0];
    if (rho < 0.9f * rho0) rho = 0.9f*rho0;

    return rho;
}

__host__ __device__ float densityShock(float* properties, float rho) {
    /**
     * \brief Mie-Gruneisen Shock Hugoniot Equation Of State
     *
     * mu = rho / rho0 -1
     * g = g * rho0 / rho
     * ph = (rho0 c0^2 mu (1 + mu)) / (1 - (s0 - 1) * mu)^2
     * uh = 1/2 ph/rho0 * (mu / (1 + mu))
     * p = ph + g * rho * (u - uh)
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * g0 = properties[2];
     * s0 = properties[3];
     * pmin = properties[4];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float rho0;
    
    rho0 = properties[0];
    if (rho < 0.9f * rho0) rho = 0.9f*rho0;

    return rho;
}


__host__ __device__ float densityTait(float* properties, float rho) {
    /**
     * \brief Tait Equation Of State
     *
     * p = rho0 * c0 * c0 / 7.0 * (powf((rho / rho0), 7) - 1.0);
     * c = c0;
     *
     * rho0 = properties[0];
     * c0 = properties[1];
     * pmin = properties[2];
     *
     * \date Jun 10, 2010
     * \author Luca Massidda
     */

    float rho0;
    
    rho0 = properties[0];
    if (rho < 0.9f * rho0) rho = 0.9f*rho0;

    return rho;
}


void updateParticlesHost(const int pn, const float alpha,
                                  const int* Material,
                                  const float* VelDotX, const float* VelDotY, const float* VelDotZ,
                                  const float* DensityDot, const float* EnergyDot,
                                  const float* PosX0, const float* PosY0, const float* PosZ0,
                                  const float* VelX0, const float* VelY0, const float* VelZ0,
                                  const float* Density0, const float* Energy0,
                                  float* PosX, float* PosY, float* PosZ,
                                  float* VelX, float* VelY, float* VelZ,
                                  float* Density, float* Energy, float* Pressure, float* Sound) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip, i;
    int iMaterial;

    for (ip = 0; ip < pn; ip++) if (Material[ip] != 0) {
        PosX[ip] = PosX0[ip] + alpha * (PosX[ip] + hRun.dt * VelX[ip] - PosX0[ip]);
        PosY[ip] = PosY0[ip] + alpha * (PosY[ip] + hRun.dt * VelY[ip] - PosY0[ip]);
        PosZ[ip] = PosZ0[ip] + alpha * (PosZ[ip] + hRun.dt * VelZ[ip] - PosZ0[ip]);

        VelX[ip] = VelX0[ip] + alpha * (VelX[ip] + hRun.dt * VelDotX[ip] - VelX0[ip]);
        VelY[ip] = VelY0[ip] + alpha * (VelY[ip] + hRun.dt * VelDotY[ip] - VelY0[ip]);
        VelZ[ip] = VelZ0[ip] + alpha * (VelZ[ip] + hRun.dt * VelDotZ[ip] - VelZ0[ip]);
		//VelZ[ip] = 0.0f;

        Density[ip] = Density0[ip] + alpha * (Density[ip] + hRun.dt * DensityDot[ip] - Density0[ip]);

        Energy[ip] = Energy0[ip] + alpha * (Energy[ip] + hRun.dt * EnergyDot[ip] - Energy0[ip]);

        iMaterial = Material[ip];

        if (iMaterial <= 0) {
            VelX[ip] = VelX0[ip];
            VelY[ip] = VelY0[ip];
            VelZ[ip] = VelZ0[ip];
        }
		
		for (i = 0; i < 10; i++) 
			if ((PosX[ip] > hFix[i].minX) && 
				(PosX[ip] < hFix[i].maxX) && 
				(PosY[ip] > hFix[i].minY) && 
				(PosY[ip] < hFix[i].maxY) && 
				(PosZ[ip] > hFix[i].minZ) && 
				(PosZ[ip] < hFix[i].maxZ)) {
					VelX[ip] = hFix[i].velX;
					VelY[ip] = hFix[i].velY;
					VelZ[ip] = hFix[i].velZ;
		}
		
        iMaterial = abs(iMaterial);

        if (hMatType[iMaterial] == 0) {
            VelX[ip] = VelX0[ip];
            VelY[ip] = VelY0[ip];
            VelZ[ip] = VelZ0[ip];
            Density[ip] = Density0[ip];
            Energy[ip] = Energy0[ip];
        }
		
        switch (hMatType[iMaterial]) {
        case (0) : // BOUNDARY
            Density[ip] = densityTait(hMatProp[iMaterial], Density[ip]);
            Pressure[ip] = 0.0f*pressureTait(hMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundTait(hMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (1) : // IDEAL GAS EOS
            Pressure[ip] = pressureGas(hMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundGas(hMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (2) : // MIE-GRUNEISEN POLYNOMIAL EOS
            Density[ip] = densityPoly(hMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressurePoly(hMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundPoly(hMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (3) : // MIE-GRUNEISEN SHOCK EOS
            Density[ip] = densityShock(hMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressureShock(hMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundShock(hMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (4) : // TAIT EOS
            Density[ip] = densityTait(hMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressureTait(hMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundTait(hMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        default :
            Pressure[ip] = 0.0f;
        }

    }
}

__global__ void updateParticlesDevice(const int pn, const float alpha,
                                      const int* Material,
                                      const float* VelDotX, const float* VelDotY, const float* VelDotZ,
                                      const float* DensityDot, const float* EnergyDot,
                                      const float* PosX0, const float* PosY0, const float* PosZ0,
                                      const float* VelX0, const float* VelY0, const float* VelZ0,
                                      const float* Density0, const float* Energy0,
                                      float* PosX, float* PosY, float* PosZ,
                                      float* VelX, float* VelY, float* VelZ,
                                      float* Density, float* Energy, float* Pressure, float* Sound) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip, i;
    int iMaterial;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < pn) {
        PosX[ip] = PosX0[ip] + alpha * (PosX[ip] + dRun.dt * VelX[ip] - PosX0[ip]);
        PosY[ip] = PosY0[ip] + alpha * (PosY[ip] + dRun.dt * VelY[ip] - PosY0[ip]);
        PosZ[ip] = PosZ0[ip] + alpha * (PosZ[ip] + dRun.dt * VelZ[ip] - PosZ0[ip]);

        VelX[ip] = VelX0[ip] + alpha * (VelX[ip] + dRun.dt * VelDotX[ip] - VelX0[ip]);
        VelY[ip] = VelY0[ip] + alpha * (VelY[ip] + dRun.dt * VelDotY[ip] - VelY0[ip]);
        VelZ[ip] = VelZ0[ip] + alpha * (VelZ[ip] + dRun.dt * VelDotZ[ip] - VelZ0[ip]);
        //VelZ[ip] = 0.0f;

        Density[ip] = Density0[ip] + alpha * (Density[ip] + dRun.dt * DensityDot[ip] - Density0[ip]);

        Energy[ip] = Energy0[ip] + alpha * (Energy[ip] + dRun.dt * EnergyDot[ip] - Energy0[ip]);

        iMaterial = Material[ip];

		for (i = 0; i < 10; i++) 
			if ((PosX[ip] > dFix[i].minX) && 
				(PosX[ip] < dFix[i].maxX) && 
				(PosY[ip] > dFix[i].minY) && 
				(PosY[ip] < dFix[i].maxY) && 
				(PosZ[ip] > dFix[i].minZ) && 
				(PosZ[ip] < dFix[i].maxZ)) {
					VelX[ip] = dFix[i].velX;
					VelY[ip] = dFix[i].velY;
					VelZ[ip] = dFix[i].velZ;
		}
		
        if (dMatType[iMaterial] == 0) {
            VelX[ip] = VelX0[ip];
            VelY[ip] = VelY0[ip];
            VelZ[ip] = VelZ0[ip];
        }

        switch (dMatType[iMaterial]) {
        case (0) : // BOUNDARY
            Density[ip] = densityTait(dMatProp[iMaterial], Density[ip]);
            Pressure[ip] = 0.0f*pressureTait(dMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundTait(dMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (1) : // IDEAL GAS EOS
            Pressure[ip] = pressureGas(dMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundGas(dMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (2) : // MIE-GRUNEISEN POLYNOMIAL EOS
            Density[ip] = densityPoly(dMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressurePoly(dMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundPoly(dMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (3) : // MIE-GRUNEISEN SHOCK EOS
            Density[ip] = densityShock(dMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressureShock(dMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundShock(dMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        case (4) : // TAIT EOS
            Density[ip] = densityTait(dMatProp[iMaterial], Density[ip]);
            Pressure[ip] = pressureTait(dMatProp[iMaterial], Density[ip], Energy[ip]);
            Sound[ip] = soundTait(dMatProp[iMaterial], Density[ip], Energy[ip]);
            break;
        default :
            Pressure[ip] = 0.0f;
        }

    }
}


void updateParticlesHostOld(const float alpha) {

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


void updateLoadsHost(const int pn,
                              const int* Material,
                              float* PosX, float* PosY, float* PosZ,
                              float* VelX, float* VelY, float* VelZ,
                              float* VelDotX, float* VelDotY, float* VelDotZ,
                              float* EnergyDot) {

    int ip, i;

    for (ip = 0; ip < pn; ip++) {
		if (Material[ip] > 0) {
			for (i = 0; i < 10; i++) {
				if ((PosX[ip] > hLoad[i].minX) &&
					(PosX[ip] < hLoad[i].maxX) &&
                    (PosZ[ip] < hLoad[i].maxZ) &&
                    (PosY[ip] > hLoad[i].minY) &&
                    (PosY[ip] < hLoad[i].maxY) &&
                    (PosZ[ip] > hLoad[i].minZ) &&
                    (PosZ[ip] < hLoad[i].maxZ)) {
					VelDotX[ip] += hLoad[i].gx;
					VelDotY[ip] += hLoad[i].gy;
					VelDotZ[ip] += hLoad[i].gz;
					EnergyDot[ip] += hLoad[i].w;
				}
            }
        }

    }

}

__global__ void updateLoadsDevice(const int pn, const int* Material, 
                                  const float* PosX, const float* PosY, const float* PosZ,
                                  float* VelDotX, float* VelDotY, float* VelDotZ,
                                  float* EnergyDot) {

    int ip, i;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if ((ip < pn) && (Material[ip] > 0)) {
        for (i = 0; i < 10; i++) {
            if ((PosX[ip] > dLoad[i].minX) &&
                    (PosX[ip] < dLoad[i].maxX) &&
                    (PosZ[ip] < dLoad[i].maxZ) &&
                    (PosY[ip] > dLoad[i].minY) &&
                    (PosY[ip] < dLoad[i].maxY) &&
                    (PosZ[ip] > dLoad[i].minZ) &&
                    (PosZ[ip] < dLoad[i].maxZ)) {
                VelDotX[ip] += dLoad[i].gx;
                VelDotY[ip] += dLoad[i].gy;
                VelDotZ[ip] += dLoad[i].gz;
                EnergyDot[ip] += dLoad[i].w;
            }
        }
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



// Host code

int initHost(struct model *hm) {

    hm->Material = (int *) malloc(MAXP * sizeof(int));
    hm->Mass = (float *) malloc(MAXP * sizeof(float));
    hm->Smooth = (float *) malloc(MAXP * sizeof(float));
    hm->PosX = (float *) malloc(MAXP * sizeof(float));
    hm->PosY = (float *) malloc(MAXP * sizeof(float));
    hm->PosZ = (float *) malloc(MAXP * sizeof(float));
    hm->VelX = (float *) malloc(MAXP * sizeof(float));
    hm->VelY = (float *) malloc(MAXP * sizeof(float));
    hm->VelZ = (float *) malloc(MAXP * sizeof(float));
    hm->Density = (float *) malloc(MAXP * sizeof(float));
    hm->Energy = (float *) malloc(MAXP * sizeof(float));
    hm->Pressure = (float *) malloc(MAXP * sizeof(float));
    hm->Sound = (float *) malloc(MAXP * sizeof(float));
    hm->VelDotX = (float *) malloc(MAXP * sizeof(float));
    hm->VelDotY = (float *) malloc(MAXP * sizeof(float));
    hm->VelDotZ = (float *) malloc(MAXP * sizeof(float));
    hm->DensityDot = (float *) malloc(MAXP * sizeof(float));
    hm->EnergyDot = (float *) malloc(MAXP * sizeof(float));
    hm->PosX0 = (float *) malloc(MAXP * sizeof(float));
    hm->PosY0 = (float *) malloc(MAXP * sizeof(float));
    hm->PosZ0 = (float *) malloc(MAXP * sizeof(float));
    hm->VelX0 = (float *) malloc(MAXP * sizeof(float));
    hm->VelY0 = (float *) malloc(MAXP * sizeof(float));
    hm->VelZ0 = (float *) malloc(MAXP * sizeof(float));
    hm->Density0 = (float *) malloc(MAXP * sizeof(float));
    hm->Energy0 = (float *) malloc(MAXP * sizeof(float));
    
    hm->Hash = (int *) malloc(MAXP * sizeof(int));
    hm->Index = (int *) malloc(MAXP * sizeof(int));
    hm->List = (int *) malloc(MAXP * MAXN * sizeof(int));
    hm->IntDummy = (int *) malloc(MAXP * sizeof(int));
    hm->FloatDummy = (float *) malloc(MAXP * sizeof(float));
	
    hm->SetStart = (int *) malloc(MAXG * sizeof(int));
    hm->SetStop = (int *) malloc(MAXG * sizeof(int));
    
    
	hMaterial = hm->Material;
	hPosX = hm->PosX;
	hPosY = hm->PosY;
    hVelX = hm->VelX;
    hVelY = hm->VelY;
    hDensity = hm->Density;
    hEnergy = hm->Energy;
    hPressure = hm->Pressure;
    hVelDotX = hm->VelDotX;
    hVelDotY = hm->VelDotY;
    hDensityDot = hm->DensityDot;
    hEnergyDot = hm->EnergyDot;
    hPosX0 = hm->PosX0;
    hPosY0 = hm->PosY0;
    hVelX0 = hm->VelX0;
    hVelY0 = hm->VelY0;
    hDensity0 = hm->Density0;
    hEnergy0 = hm->Energy0;
	
    hList = hm->List;
    hHash = hm->Hash;
    hIndex = hm->Index;
    hSetStart = hm->SetStart;
    hSetStop = hm->SetStop;
	
    return 0;
}

int initDevice(struct model *dm) {
	size_t available, total;
	
    cudaMalloc((void**) &(dm->Material), (MAXP * sizeof(int)));
    cudaMalloc((void**) &(dm->Mass), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Smooth), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosX), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosY), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosZ), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelX), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelY), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelZ), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Density), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Energy), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Pressure), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Sound), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelDotX), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelDotY), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelDotZ), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->DensityDot), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->EnergyDot), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosX0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosY0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->PosZ0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelX0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelY0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->VelZ0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Density0), (MAXP * sizeof(float)));
    cudaMalloc((void**) &(dm->Energy0), (MAXP * sizeof(float)));
    
    cudaMalloc((void**) &(dm->Hash), (MAXP * sizeof(int)));
    cudaMalloc((void**) &(dm->Index), (MAXP * sizeof(int)));
    cudaMalloc((void**) &(dm->List), (MAXP * MAXN * sizeof(int)));
    cudaMalloc((void**) &(dm->IntDummy), (MAXP * sizeof(int)));
    cudaMalloc((void**) &(dm->FloatDummy), (MAXP * sizeof(float)));
	
    cudaMalloc((void**) &(dm->SetStart), (MAXG * sizeof(int)));
    cudaMalloc((void**) &(dm->SetStop), (MAXG * sizeof(int)));
    
	cudaMemGetInfo(&available, &total);
	printf("Available memory %d of %d MB\n", available/1024/1024, total/1024/1024);
	
    return 0;
}

int copyHostToDevice(struct model *hm, struct model *dm) {

    dm->pn = hm->pn;
    
	cudaMemcpy(dm->Material, hm->Material, (MAXP * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Mass, hm->Mass, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Smooth, hm->Smooth, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosX, hm->PosX, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosY, hm->PosY, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosZ, hm->PosZ, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelX, hm->VelX, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelY, hm->VelY, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelZ, hm->VelZ, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Density, hm->Density, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Energy, hm->Energy, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Pressure, hm->Pressure, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Sound, hm->Sound, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelDotX, hm->VelDotX, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelDotY, hm->VelDotY, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelDotZ, hm->VelDotZ, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->DensityDot, hm->DensityDot, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->EnergyDot, hm->EnergyDot, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosX0, hm->PosX0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosY0, hm->PosY0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->PosZ0, hm->PosZ0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelX0, hm->VelX0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelY0, hm->VelY0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->VelZ0, hm->VelZ0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Density0, hm->Density0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Energy0, hm->Energy0, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);

    cudaMemcpy(dm->List, hm->List, (MAXP * MAXN * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Hash, hm->Hash, (MAXP * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->Index, hm->Index, (MAXP * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->IntDummy, hm->IntDummy, (MAXP * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->FloatDummy, hm->FloatDummy, (MAXP * sizeof(float)), cudaMemcpyHostToDevice);

    cudaMemcpy(dm->SetStart, hm->SetStart, (MAXG * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->SetStop, hm->SetStop, (MAXG * sizeof(int)), cudaMemcpyHostToDevice);

    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.oZ = hGrid.oZ;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.nZ = hGrid.nZ;
    dGrid.size = hGrid.size;
    
    for (int i = 0; i < 10; i++) dIn[i] = hIn[i];

    cudaMemcpyToSymbol("dMatType", hMatType, 10 * sizeof(int));
    cudaMemcpyToSymbol("dMatProp", hMatProp, 100 * sizeof(float));
    cudaMemcpyToSymbol("dRun", &hRun, sizeof(struct simulation));
    cudaMemcpyToSymbol("dLoad", &hLoad, 10 * sizeof(struct load));
    cudaMemcpyToSymbol("dFix", &hFix, 10 * sizeof(struct fix));
    cudaMemcpyToSymbol("dOut", &hOut, 10 * sizeof(struct outlet));
    
    return 0;
}

int copyDeviceToHost(struct model *dm, struct model *hm) {

    hm->pn = dm->pn;
    
	cudaMemcpy(hm->Material, dm->Material, (MAXP * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Mass, dm->Mass, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Smooth, dm->Smooth, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosX, dm->PosX, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosY, dm->PosY, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosZ, dm->PosZ, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelX, dm->VelX, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelY, dm->VelY, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelZ, dm->VelZ, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Density, dm->Density, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Energy, dm->Energy, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Pressure, dm->Pressure, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Sound, dm->Sound, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelDotX, dm->VelDotX, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelDotY, dm->VelDotY, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelDotZ, dm->VelDotZ, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->DensityDot, dm->DensityDot, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->EnergyDot, dm->EnergyDot, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosX0, dm->PosX0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosY0, dm->PosY0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->PosZ0, dm->PosZ0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelX0, dm->VelX0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelY0, dm->VelY0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->VelZ0, dm->VelZ0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Density0, dm->Density0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Energy0, dm->Energy0, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);

    cudaMemcpy(hm->List, dm->List, (MAXP * MAXN * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Hash, dm->Hash, (MAXP * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->Index, dm->Index, (MAXP * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->IntDummy, dm->IntDummy, (MAXP * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->FloatDummy, dm->FloatDummy, (MAXP * sizeof(float)), cudaMemcpyDeviceToHost);

    cudaMemcpy(hm->SetStart, dm->SetStart, (MAXG * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->SetStop, dm->SetStop, (MAXG * sizeof(int)), cudaMemcpyDeviceToHost);

    hGrid.oX = dGrid.oX;
    hGrid.oY = dGrid.oY;
    hGrid.oZ = dGrid.oZ;
    hGrid.nX = dGrid.nX;
    hGrid.nY = dGrid.nY;
    hGrid.nZ = dGrid.nZ;
    hGrid.size = dGrid.size;
    
    for (int i = 0; i < 10; i++) hIn[i] = dIn[i];

    return 0;
}


int backupDataHost(struct model *hm) {

    memcpy(hm->PosX0, hm->PosX, MAXP * sizeof(float));
    memcpy(hm->PosY0, hm->PosY, MAXP * sizeof(float));
    memcpy(hm->PosZ0, hm->PosZ, MAXP * sizeof(float));
    memcpy(hm->VelX0, hm->VelX, MAXP * sizeof(float));
    memcpy(hm->VelY0, hm->VelY, MAXP * sizeof(float));
    memcpy(hm->VelZ0, hm->VelZ, MAXP * sizeof(float));
    memcpy(hm->Density0, hm->Density, MAXP * sizeof(float));
    memcpy(hm->Energy0, hm->Energy, MAXP * sizeof(float));

    return 0;
}

int backupDataDevice(struct model *dm) {

    cudaMemcpy(dm->PosX0, dm->PosX, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->PosY0, dm->PosY, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->PosZ0, dm->PosZ, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->VelX0, dm->VelX, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->VelY0, dm->VelY, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->VelZ0, dm->VelZ, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->Density0, dm->Density, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    cudaMemcpy(dm->Energy0, dm->Energy, (MAXP * sizeof(float)), cudaMemcpyDeviceToDevice);
    
    return 0;
}


int backupDataHostOld() {

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

int printDataOld() {
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

int printData(struct model *hm) {
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
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e \n", hm->PosX[i], hm->PosY[i], hm->PosZ[i]);
    fclose(stream);

    // Stream file velocity
    stream = fopen("new_vel.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e \n", hm->VelX[i], hm->VelY[i], hm->VelZ[i]);
    fclose(stream);

    // Stream file info
    stream = fopen("new_info.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%i %+14.8e %+14.8e \n", hm->Material[i], hm->Mass[i], hm->Smooth[i]);
    fclose(stream);

    // Stream file field
    stream = fopen("new_field.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+14.8e %+14.8e %+14.8e \n", hm->Density[i], hm->Pressure[i], hm->Energy[i]);
    fclose(stream);
	/*
    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%d %d %d %f %f %f\n", i, hm->Index[i], hm->Hash[i], hm->PosX[i], hm->PosY[i], hm->PosZ[i]);
    fclose(stream);
	*/
	/*
    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%d %f %f %f %f %f %f\n", i, hm->VelX[i], hm->VelY[i], hm->VelZ[i], hm->Density[i], hm->Energy[i], hm->Pressure[i]);
    fclose(stream);
    */
    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%d %d %d %d %d\n", i, hm->Index[i], hm->Hash[i], hm->SetStart[hm->Hash[i]], hm->SetStop[hm->Hash[i]]);
    fclose(stream);
    /*
    for (i = 0; i < hm->pn; i++) {
    	printf("%d - ", i);
    	for (int j = 0; j < MAXN; j++)
    		printf("%d ", hm->List[i * MAXN +j]);
    	printf("\n");
    }
    */

    return 0;
}


int outputVTK(struct model *hm, int ss) {
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

    fprintf(stream, "POINTS %i float\n", hm->pn);
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+e %+e %+e \n", hm->PosX[i], hm->PosY[i], hm->PosZ[i]);

    fprintf(stream, "CELLS %i %i \n", hm->pn, 2*hm->pn);
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%i %i \n", 1, i);

    fprintf(stream, "CELL_TYPES %i \n", hm->pn);
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%i \n", 1);

    fprintf(stream, "POINT_DATA %i \n", hm->pn);

    fprintf(stream, "SCALARS material int 1 \n", hm->pn);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+d \n", hm->Material[i]);

    fprintf(stream, "SCALARS density float 1 \n", hm->pn);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+e \n", hm->Density[i]);

    fprintf(stream, "SCALARS pressure float 1 \n", hm->pn);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+e \n", hm->Pressure[i]);

    fprintf(stream, "SCALARS energy float 1 \n", hm->pn);
    fprintf(stream, "LOOKUP_TABLE default\n");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+e \n", hm->Energy[i]);

    fprintf(stream, "VECTORS velocity float\n");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%+e %+e %+e \n", hm->VelX[i], hm->VelY[i], hm->VelZ[i]);

    fclose(stream);
    /*
    for (i = 0; i < hm->pn; i++) printf("%d %d \n", i, hm->Hash[i]);
    printf("\n\n\n");
    for (i = 0; i < hm->SetStart.size(); i++) printf("%d %d %d \n", i, hm->SetStart[i], hm->SetStop[i]);

    for (i = 0; i < hm->pn; i++) {
    	printf("%d  -  ", i);
    	for (j = 0; j < MAXN; j++) printf("%d ", hm->List[i*MAXN +j]);
    	printf("\n");
    }
    */
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


int outputVTKOld(int ss) {
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

    dr = 0.05; // x4
    pi = 0;

    for (j = 0; j < 20; j++) {
        for (i = 0; i < 20; i++) {
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
    
    for (j = -3; j < -1; j++) {
        for (i = 0; i < 20; i++) {
            hPosX[pi] = i * dr + 0.0 * dr;
            hPosY[pi] = j * dr + 0.0 * dr;

            hVelX[pi] = 0.0;
            hVelY[pi] = 0.0;
            hMaterial[pi] = -m;
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

    hRun.minX = -1.5;
    hRun.maxX =  2.5;
    hRun.minY = -1.5;
    hRun.maxY =  2.5;
    hRun.minZ = -1.5;
    hRun.maxZ =  2.5;

    hRun.dt = 1e-3; //1.0e-3;
    hRun.tsn = 600; //1000;
    hRun.ssi = 200;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0 * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    

    printf("Freefall\n");
    printf("Particles: %i \n", hPN);
}


void initBox(struct model *hm) {

    int i, j, k, m, b, pi;
    double rho, c0, pmin;
    double dr;

    m = 1;
    b = 2;
    rho = 1000.;
    c0 = 50.;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 0;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = c0;
    hMatProp[b][2] = pmin;

    dr = 0.1; // x4
    pi = 0;

    for (k = 0; k < 10; k++) {
		for (j = 0; j < 10; j++) {
			for (i = 0; i < 10; i++) {
				hm->PosX[pi] = i * dr + 0.0 * dr;
				hm->PosY[pi] = j * dr + 0.0 * dr;
				hm->PosZ[pi] = k * dr + 0.0 * dr;
	
				hm->VelX[pi] = 0.0;
				hm->VelY[pi] = 0.0;
				hm->VelZ[pi] = 0.0;
				hm->Material[pi] = m;
				hm->Density[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
				hm->Energy[pi] = 0.0;
				hm->Pressure[pi] = 1.0;
				pi++;
			}
        }
    }
    
    for (k = 0; k < 10; k++) {
		for (j = -2; j < -1; j++) {
			for (i = 0; i < 10; i++) {
				hm->PosX[pi] = i * dr + 0.0 * dr;
				hm->PosY[pi] = j * dr + 0.0 * dr;
				hm->PosZ[pi] = k * dr + 0.0 * dr;
				
				hm->VelX[pi] = 0.0;
				hm->VelY[pi] = 0.0;
				hm->VelZ[pi] = 0.0;
				hm->Material[pi] = b;
				hm->Density[pi] = rho; //+ (9.81 * rho / c0 / c0 * (50 - j) * dr);
				hm->Energy[pi] = 0.0;
				hm->Pressure[pi] = 1.0;
				pi++;
			}
        }
    }
    
    hPN = pi;
    hSmooth = 1.2 * dr;
    hMass = rho * dr * dr * dr;
    hSound = c0;

    hRun.minX = -2.5;
    hRun.maxX =  2.5;
    hRun.minY = -2.5;
    hRun.maxY =  2.5;
    hRun.minZ = -2.5;
    hRun.maxZ =  2.5;

    hRun.dt = 2.0e-3; //1.0e-3;
    hRun.tsn = 600; //1000;
    hRun.ssi = 200;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0 * hSmooth;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    
	hm->pn = hPN;
	for (i = 0; i < hPN; i++) {
		hm->Mass[i] = hMass;
		hm->Smooth[i] = hSmooth;
	}

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

    // Particles are re ordered
    
    //iSort(hHash, hIndex, hPN);
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

int mapCompare(const void *a, const void *b)
{
	int c;
	struct pair m1, m2;
	
	c = 0;
	m1 = *(struct pair*)a;
	m2 = *(struct pair*)b;
	if (m1.key < m2.key) c = -1;
	if (m1.key > m2.key) c = 1;
  return c;
}



int updateHashHostOld() {
	
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

void updateHashHost(const int pn, const struct grid Grid,
                             const float* PosX, const float* PosY, const float* PosZ,
                             int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip, ix, iy, iz, ic;

    for (ip = 0; ip < pn; ip++) {
        ix = (int) truncf((PosX[ip] - Grid.oX) / Grid.size);
        iy = (int) truncf((PosY[ip] - Grid.oY) / Grid.size);
        iz = (int) truncf((PosZ[ip] - Grid.oZ) / Grid.size);
        ic = ix + iy * Grid.nX + iz * Grid.nX * Grid.nY;
		
		if (ic < 0) ic = 0;
		if (ic > Grid.nX * Grid.nY * Grid.nZ) ic = ic = Grid.nX * Grid.nY * Grid.nZ;
        Hash[ip] = ic;
    }
}


__global__ void updateHashDevice(const int pn, const struct grid Grid,
                                 const float* PosX, const float* PosY, const float* PosZ,
                                 int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip, ix, iy, iz, ic;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < pn) {
        ix = (int) truncf((PosX[ip] - Grid.oX) / Grid.size);
        iy = (int) truncf((PosY[ip] - Grid.oY) / Grid.size);
        iz = (int) truncf((PosZ[ip] - Grid.oZ) / Grid.size);
        ic = ix + iy * Grid.nX + iz * Grid.nX * Grid.nY;

		if (ic < 0) ic = 0;
		if (ic > Grid.nX * Grid.nY * Grid.nZ) ic = ic = Grid.nX * Grid.nY * Grid.nZ;
        Hash[ip] = ic;
    }
}

int updateSetsHostOld() {
	
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


void updateSetsHost(const int pn, int *SetStart, int *SetStop,
                               const int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip;
    int hash, nextHash, prevHash;
    
    for (ip = 0; ip < pn; ip++) {
		hash = Hash[ip];
		if (ip == 0) prevHash = -1;
		else prevHash = Hash[ip -1];
		if (ip == pn -1) nextHash = -1;
		else nextHash = Hash[ip +1];
		
		if (hash != prevHash) SetStart[hash] = ip;
		if (hash != nextHash) SetStop[hash] = ip +1;
	}

}


__global__ void updateSetsDevice(const int pn, int *SetStart, int *SetStop,
                                 const int* Hash) {

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
    if (ip >= pn) return;

    hash = Hash[ip];

    if (threadIdx.x < THREADS -1) prevHash[threadIdx.x +1] = hash;
    if (threadIdx.x > 0) nextHash[threadIdx.x -1] = hash;

    if (threadIdx.x == 0) {
        if (ip == 0) prevHash[threadIdx.x] = -1;
        else prevHash[threadIdx.x] = Hash[ip -1];
    }

    if (threadIdx.x == THREADS -1) {
        if (ip == pn -1) nextHash[threadIdx.x] = -1;
        else nextHash[threadIdx.x] = Hash[ip +1];
    }

    __syncthreads();

    if (hash != prevHash[threadIdx.x]) SetStart[hash] = ip;

    if (hash != nextHash[threadIdx.x]) SetStop[hash] = ip +1;

}


int updateListHostOld(void) {
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


void updateListHost(const int pn, int *List,
                             const int* SetStart, const int* SetStop,
                             const struct grid Grid, const float* Smooth,
                             const float* PosX, const float* PosY, const float* PosZ) {

    int ip, ic, ix, iy, iz, i, j, k, jp, jc, np;
    float dx, dy, dz, dr;

    // Particles list is filled
    for (ip = 0; ip < pn; ip++) {
		ix = (int) ((PosX[ip] - Grid.oX) / Grid.size);
		iy = (int) ((PosY[ip] - Grid.oY) / Grid.size);
		iz = (int) ((PosZ[ip] - Grid.oZ) / Grid.size);
		ic = ix + iy * Grid.nX + iz * Grid.nX * Grid.nY;
		np = 0;
		
		for (k = -1; k <= 1; k++) {
			for (j = -1; j <= 1; j++) {
				for (i = -1; i <= 1; i++) {
					jc = ic + i + j * Grid.nX + k * Grid.nX * Grid.nY;
					
					if (jc >= 0 && jc <= Grid.nX * Grid.nY * Grid.nZ) {
						for (jp = SetStart[jc]; jp < SetStop[jc]; jp++) {
							dx = PosX[ip] - PosX[jp];
							dy = PosY[ip] - PosY[jp];
							dz = PosZ[ip] - PosZ[jp];
							dr = sqrtf(dx * dx + dy * dy + dz * dz);
							
							if ((dr < 2.0f * Smooth[ip]) && (np < MAXN)) {
								List[ip * MAXN + np] = jp;
								np++;
							}
						}
					}
				}
			}
		}
		
		while (np < MAXN) {
			List[ip * MAXN + np] = ip;
			np++;
		}
	}
}


__global__ void updateListDevice(const int pn, int *List,
                                 const int* SetStart, const int* SetStop,
                                 const struct grid Grid, const float* Smooth,
                                 const float* PosX, const float* PosY, const float* PosZ) {

    int ip, ic, ix, iy, iz, i, j, k, jp, jc, np;
    float dx, dy, dz, dr;

    // Particles list is filled
    ip = threadIdx.x + blockDim.x * blockIdx.x;
    if (ip >= pn) return;

    ix = (int) ((PosX[ip] - Grid.oX) / Grid.size);
    iy = (int) ((PosY[ip] - Grid.oY) / Grid.size);
    iz = (int) ((PosZ[ip] - Grid.oZ) / Grid.size);
    ic = ix + iy * Grid.nX + iz * Grid.nX * Grid.nY;
    np = 0;

    for (k = -1; k <= 1; k++) {
        for (j = -1; j <= 1; j++) {
            for (i = -1; i <= 1; i++) {
                jc = ic + i + j * Grid.nX + k * Grid.nX * Grid.nY;
				
				if (jc >= 0 && jc <= Grid.nX * Grid.nY * Grid.nZ) {
					for (jp = SetStart[jc]; jp < SetStop[jc]; jp++) {
						dx = PosX[ip] - PosX[jp];
						dy = PosY[ip] - PosY[jp];
						dz = PosZ[ip] - PosZ[jp];
						dr = sqrtf(dx * dx + dy * dy + dz * dz);
						
						if ((dr < 2.0f * Smooth[ip]) && (np < MAXN)) {
							List[ip * MAXN + np] = jp;
							np++;
						}
					}
                }
            }
        }
    }

    while (np < MAXN) {
        List[ip * MAXN + np] = ip;
        np++;
    }

}

int neighbourListHost(struct model *hm) {
	struct pair map[MAXP];
	int i, ip;
	
	//updateHashHostOld();
    updateHashHost(hm->pn, hGrid, hm->PosX, hm->PosY, hm->PosZ, hm->Hash);
    for (ip = 0; ip < hm->pn; ip++) hm->Index[ip] = ip;
	
	//qsort(hIndex, hPN, sizeof(int), indexCompare);
    for (ip = 0; ip < hm->pn; ip++) {
		map[ip].key = hm->Hash[ip];
		map[ip].value = hm->Index[ip];
	}
	qsort(map, hm->pn, sizeof(struct pair), mapCompare);
    for (ip = 0; ip < hm->pn; ip++) {
		hm->Hash[ip] = map[ip].key;
		hm->Index[ip] = map[ip].value;
	}
	
	//sortArraysHost();
    iSort(hMaterial, hIndex, hPN);
    fSort(hPosX, hIndex, hPN);
    fSort(hPosY, hIndex, hPN);
    fSort(hVelX, hIndex, hPN);
    fSort(hVelY, hIndex, hPN);
    fSort(hDensity, hIndex, hPN);
    fSort(hEnergy, hIndex, hPN);
    fSort(hPressure, hIndex, hPN);
    
    fSort(hm->Mass, hm->Index, hm->pn);
    fSort(hm->Smooth, hm->Index, hm->pn);
    fSort(hm->PosZ, hm->Index, hm->pn);
    fSort(hm->VelZ, hm->Index, hm->pn);
	
	for (i = 0; i < hGrid.nX * hGrid.nY * hGrid.nZ; i++) hm->SetStart[i] = 0;
	for (i = 0; i < hGrid.nX * hGrid.nY * hGrid.nZ; i++) hm->SetStop[i] = 0;
	
	updateSetsHost(hm->pn, hm->SetStart, hm->SetStop, hm->Hash);
	
	//updateSetsHostOld();
	
	updateListHost(hm->pn, hm->List, hm->SetStart, hm->SetStop, hGrid, hm->Smooth,
				   hm->PosX, hm->PosY, hm->PosZ);
	//updateListHostOld();
	
	return 0;
}

int neighbourListDevice(struct model *dm) {
    int blocks, threads;

    blocks = (dm->pn + THREADS - 1) / THREADS;
    threads = THREADS;
    
	thrust::device_ptr<int> tIndex(dm->Index);
	thrust::device_ptr<int> tHash(dm->Hash);
	thrust::device_ptr<int> tMaterial(dm->Material);
	thrust::device_ptr<float> tMass(dm->Mass);
	thrust::device_ptr<float> tSmooth(dm->Smooth);
	thrust::device_ptr<float> tPosX(dm->PosX);
	thrust::device_ptr<float> tPosY(dm->PosY);
	thrust::device_ptr<float> tPosZ(dm->PosZ);
	thrust::device_ptr<float> tVelX(dm->VelX);
	thrust::device_ptr<float> tVelY(dm->VelY);
	thrust::device_ptr<float> tVelZ(dm->VelZ);
	thrust::device_ptr<float> tDensity(dm->Density);
	thrust::device_ptr<float> tEnergy(dm->Energy);
	thrust::device_ptr<int> tIntDummy(dm->IntDummy);
	thrust::device_ptr<float> tFloatDummy(dm->FloatDummy);
	thrust::device_ptr<int> tSetStart(dm->SetStart);
	thrust::device_ptr<int> tSetStop(dm->SetStop);
	
	updateHashDevice <<< blocks, threads >>>
	(dm->pn, dGrid, dm->PosX, dm->PosY, dm->PosZ, dm->Hash);
	
	thrust::sequence(tIndex, tIndex + dm->pn, 1);
	
	thrust::sort_by_key(tHash, tHash + dm->pn, tIndex);
	
	thrust::copy(tMaterial, tMaterial + dm->pn, tIntDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tIntDummy, tMaterial);
	thrust::copy(tMass, tMass + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tMass);
	thrust::copy(tSmooth, tSmooth + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tSmooth);
	thrust::copy(tPosX, tPosX + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tPosX);
	thrust::copy(tPosY, tPosY + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tPosY);
	thrust::copy(tPosZ, tPosZ + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tPosZ);
	thrust::copy(tVelX, tVelX + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tVelX);
	thrust::copy(tVelY, tVelY + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tVelY);
	thrust::copy(tVelZ, tVelZ + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tVelZ);
	thrust::copy(tDensity, tDensity + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tDensity);
	thrust::copy(tEnergy, tEnergy + dm->pn, tFloatDummy);
	thrust::gather(tIndex, tIndex + dm->pn, tFloatDummy, tEnergy);
	
	thrust::fill(tSetStart, tSetStart + dm->pn, 0);
	thrust::fill(tSetStop, tSetStop + dm->pn, 0);
	
	updateSetsDevice <<< blocks, threads >>>
	(dm->pn, dm->SetStart, dm->SetStop, dm->Hash);
	
	updateListDevice <<< blocks, threads >>>
	(dm->pn, dm->List, dm->SetStart, dm->SetStop, dGrid, dm->Smooth,
	dm->PosX, dm->PosY, dm->PosZ);
	
    return 0;
}

int RKstepHost(struct model *hm, float alpha) {
	int ip;
	
	for (ip = 0; ip < hm->pn; ip++) {
		hm->VelDotX[ip] = 0.0f;
		hm->VelDotY[ip] = 0.0f;
		hm->VelDotZ[ip] = 0.0f;
		hm->DensityDot[ip] = 0.0f;
		hm->EnergyDot[ip] = 0.0f;
	}
	
    // External loads
    updateLoadsHost(hm->pn, hm->Material, 
                    hm->PosX, hm->PosY, hm->PosZ,
                    hm->VelX, hm->VelY, hm->VelZ, 
                    hm->VelDotX, hm->VelDotY, hm->VelDotZ, hm->EnergyDot);
	
	// External forces
	//updateForcesHost();
	
	// Calculate particle interactions
    balanceMassMomentumHost(hm->pn, hm->List, hm->Material, hm->Mass, hm->Smooth, 
                            hm->PosX, hm->PosY, hm->PosZ, 
                            hm->VelX, hm->VelY, hm->VelZ, 
                            hm->Density, hm->Pressure, hm->Sound, 
                            hm->DensityDot, hm->VelDotX, hm->VelDotY, hm->VelDotZ);
	//balanceMassMomentumHostOld();
	
    balanceEnergyHost(hm->pn, hm->Pressure, hm->Density, 
                      hm->DensityDot, hm->EnergyDot);
	//balanceEnergyHostOld();
	
    // Update particles
    updateParticlesHost(hm->pn, alpha, hm->Material, 
                        hm->VelDotX, hm->VelDotY, hm->VelDotZ, hm->DensityDot, hm->EnergyDot,
                        hm->PosX0, hm->PosY0, hm->PosZ0, 
                        hm->VelX0, hm->VelY0, hm->VelZ0, hm->Density0, hm->Energy0, 
                        hm->PosX, hm->PosY, hm->PosZ, hm->VelX, hm->VelY, hm->VelZ, 
                        hm->Density, hm->Energy, hm->Pressure, hm->Sound);
	//updateParticlesHostOld(alpha);
	
	return 0;
}

int RKstepDevice(struct model *dm, float alpha) {
    int blocks, threads;

    blocks = (dm->pn + THREADS - 1) / THREADS;
    threads = THREADS;
	
	thrust::device_ptr<float> tVelDotX(dm->VelDotX);
	thrust::device_ptr<float> tVelDotY(dm->VelDotY);
	thrust::device_ptr<float> tVelDotZ(dm->VelDotZ);
	thrust::device_ptr<float> tDensityDot(dm->DensityDot);
	thrust::device_ptr<float> tEnergyDot(dm->EnergyDot);
	
	thrust::fill(tVelDotX, tVelDotX + dm->pn, 0.0f);
	thrust::fill(tVelDotY, tVelDotY + dm->pn, 0.0f);
	thrust::fill(tVelDotZ, tVelDotZ + dm->pn, 0.0f);
	thrust::fill(tDensityDot, tDensityDot + dm->pn, 0.0f);
	thrust::fill(tEnergyDot, tEnergyDot + dm->pn, 0.0f);
	
	// External loads
	updateLoadsDevice <<< blocks, threads >>>
	(dm->pn, dm->Material, dm->PosX, dm->PosY, dm->PosZ, 
	dm->VelDotX, dm->VelDotY, dm->VelDotZ, dm->EnergyDot);
	
	// Calculate particle interactions
	balanceMassMomentumDevice <<< blocks, threads >>>
	(dm->pn, dm->List, dm->Material, dm->Mass, dm->Smooth, dm->PosX, dm->PosY, dm->PosZ,
	dm->VelX, dm->VelY, dm->VelZ, dm->Density, dm->Pressure, dm->Sound,
	dm->DensityDot, dm->VelDotX, dm->VelDotY, dm->VelDotZ);
	
	balanceEnergyDevice <<< blocks, threads >>>
	(dm->pn, dm->Pressure, dm->Density, dm->DensityDot, dm->EnergyDot);
	
	// Update particles
	updateParticlesDevice  <<< blocks, threads >>>
	(dm->pn, alpha, dm->Material, dm->VelDotX, dm->VelDotY, dm->VelDotZ, dm->DensityDot, dm->EnergyDot,
	dm->PosX0, dm->PosY0, dm->PosZ0, dm->VelX0, dm->VelY0, dm->VelZ0, dm->Density0, dm->Energy0,
	dm->PosX, dm->PosY, dm->PosZ, dm->VelX, dm->VelY, dm->VelZ, dm->Density, dm->Energy, dm->Pressure, dm->Sound);
	
    return 0;
}


int RKintegrateHost(struct model *hm) {

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
            printData(hm);
            outputVTK(hm, ts / hRun.ssi);
        }
        
		// Calculate neighbour list
		neighbourListHost(hm);
		
		// Save initial condition
		backupDataHost(hm);
		
        // Step 1
		RKstepHost(hm, 1.0);
		
        // Step 2
		RKstepHost(hm, 1.0 / 4.0);
		
        // Step 3
		RKstepHost(hm, 2.0 / 3.0);
		
	}
	
	return 0;
}


int RKintegrateDevice(struct model *hm, struct model *dm) {

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
            copyDeviceToHost(dm, hm);
            printf("Saving time: %g \n", ts * hRun.dt);
            printData(hm);
            outputVTK(hm, ts / hRun.ssi);
        }
        
		// Calculate neighbour list
		neighbourListDevice(dm);
		
		// Save initial condition
		backupDataDevice(dm);
		
        // Step 1
		RKstepDevice(dm, 1.0);
		
        // Step 2
		RKstepDevice(dm, 1.0 / 4.0);
		
        // Step 3
		RKstepDevice(dm, 2.0 / 3.0);
		
	}
	
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

    struct model hModel, dModel;
    int i;
    
    initHost(&hModel);
    for (i = 0; i < 10; i++) {
        hLoad[i].gx = 0.0f;
        hLoad[i].gy = 0.0f;
        hLoad[i].gz = 0.0f;
        hLoad[i].w = 0.0f;

        hOut[i].nX = 0.0f;
        hOut[i].nY = 0.0f;
        hOut[i].nZ = 0.0f;
    }

    initBox(&hModel);
	
    initDevice(&dModel);
    copyHostToDevice(&hModel, &dModel);
	//RKintegrateDevice(&hModel, &dModel);
	
	RKintegrateHost(&hModel);
	
    return 0;
}
