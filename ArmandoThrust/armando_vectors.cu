#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/count.h>

//#define MAXN 24
#define MAXN 96
#define PI 3.14159f

#define THREADS 256

struct grid {
    float oX, oY, oZ;
    float size;
    int nX, nY, nZ;
};

struct simulation {
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

struct host_model {
    int pn;
    thrust::host_vector<int> Material;
    thrust::host_vector<float> Mass;
    thrust::host_vector<float> Smooth;
    thrust::host_vector<float> PosX;
    thrust::host_vector<float> PosY;
    thrust::host_vector<float> PosZ;
    thrust::host_vector<float> VelX;
    thrust::host_vector<float> VelY;
    thrust::host_vector<float> VelZ;
    thrust::host_vector<float> Density;
    thrust::host_vector<float> Energy;
    thrust::host_vector<float> Pressure;
    thrust::host_vector<float> Sound;
    thrust::host_vector<float> VelDotX;
    thrust::host_vector<float> VelDotY;
    thrust::host_vector<float> VelDotZ;
    thrust::host_vector<float> DensityDot;
    thrust::host_vector<float> EnergyDot;
    thrust::host_vector<float> PosX0;
    thrust::host_vector<float> PosY0;
    thrust::host_vector<float> PosZ0;
    thrust::host_vector<float> VelX0;
    thrust::host_vector<float> VelY0;
    thrust::host_vector<float> VelZ0;
    thrust::host_vector<float> Density0;
    thrust::host_vector<float> Energy0;
    thrust::host_vector<int> List;
    thrust::host_vector<int> Hash;
    thrust::host_vector<int> Index;
    thrust::host_vector<int> SetStart;
    thrust::host_vector<int> SetStop;
    thrust::host_vector<int> IntDummy;
    thrust::host_vector<float> FloatDummy;
};

struct device_model {
    int pn;
    thrust::device_vector<int> Material;
    thrust::device_vector<float> Mass;
    thrust::device_vector<float> Smooth;
    thrust::device_vector<float> PosX;
    thrust::device_vector<float> PosY;
    thrust::device_vector<float> PosZ;
    thrust::device_vector<float> VelX;
    thrust::device_vector<float> VelY;
    thrust::device_vector<float> VelZ;
    thrust::device_vector<float> Density;
    thrust::device_vector<float> Energy;
    thrust::device_vector<float> Pressure;
    thrust::device_vector<float> Sound;
    thrust::device_vector<float> VelDotX;
    thrust::device_vector<float> VelDotY;
    thrust::device_vector<float> VelDotZ;
    thrust::device_vector<float> DensityDot;
    thrust::device_vector<float> EnergyDot;
    thrust::device_vector<float> PosX0;
    thrust::device_vector<float> PosY0;
    thrust::device_vector<float> PosZ0;
    thrust::device_vector<float> VelX0;
    thrust::device_vector<float> VelY0;
    thrust::device_vector<float> VelZ0;
    thrust::device_vector<float> Density0;
    thrust::device_vector<float> Energy0;
    thrust::device_vector<int> List;
    thrust::device_vector<int> Hash;
    thrust::device_vector<int> Index;
    thrust::device_vector<int> SetStart;
    thrust::device_vector<int> SetStop;
    thrust::device_vector<int> IntDummy;
    thrust::device_vector<float> FloatDummy;
};

struct pointer_model {
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

int hMatType[10];
float hMatProp[10][10];
struct simulation hRun;
struct grid hGrid;
struct load hLoad[10];
struct fix hFix[10];
struct outlet hOut[10];
struct inlet hIn[10];

// Device Variables

__device__ __constant__ int dMatType[10];
__device__ __constant__ float dMatProp[10][10];
__device__ __constant__ struct simulation dRun;
__device__ struct grid dGrid;
__device__ __constant__ struct load dLoad[10];
__device__ __constant__ struct fix dFix[10];
__device__ __constant__ struct outlet dOut[10];
__device__ struct inlet dIn[10];


void initPump(struct host_model *hm) {

    int i, j, k, m, p;
    double rho, c0, pmin;
    double dr;

    p = 1;

    m = 1;
    rho = 1000.f;
    c0 = 50.f;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    dr = 0.02f / p;

    for (k = 0; k < 20 * p ; k++) {
        for (j = 0; j < 20 * p ; j++) {
            for (i = 0; i < 40 * p; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = 20 * p; j < 50 * p ; j++) {
            for (i = 0; i < 20 * p -1; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = 20 * p; j < 50 * p ; j++) {
            for (i = 20 * p +1; i < 40 * p; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = -3; j < 0; j++) {
            for (i = -3; i < 40 * p + 3; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = 0; j < 80 * p; j++) {
            for (i = -3; i < 0; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = 0; j < 80 * p; j++) {
            for (i = 40 * p; i < 40 * p +3; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    for (k = 0; k < 20 * p ; k++) {
        for (j = 20 * p; j < 60 * p; j++) {
            for (i = 20 * p -1; i < 20 * p +1; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    for (k = -3; k < 0 ; k++) {
        for (j = -3; j < 80 * p ; j++) {
            for (i = -3; i < 40 * p +3; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    for (k = 20 * p; k < 20 * p +3; k++) {
        for (j = -3; j < 80 * p ; j++) {
            for (i = -3; i < 40 * p +3; i++) {
                hm->PosX.push_back(i * dr);
                hm->PosY.push_back(j * dr);
                hm->PosZ.push_back(k * dr);
                hm->Material.push_back(-m);
            }
        }
    }

    hm->pn = hm->Material.size();
    
    hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hm->VelDotX.resize(hm->pn);
    hm->VelDotY.resize(hm->pn);
    hm->VelDotZ.resize(hm->pn);
    hm->DensityDot.resize(hm->pn);
    hm->EnergyDot.resize(hm->pn);

    hm->PosX0.resize(hm->pn);
    hm->PosY0.resize(hm->pn);
    hm->PosZ0.resize(hm->pn);
    hm->VelX0.resize(hm->pn);
    hm->VelY0.resize(hm->pn);
    hm->VelZ0.resize(hm->pn);
    hm->Density0.resize(hm->pn);
    hm->Energy0.resize(hm->pn);

    hm->List.resize(hm->pn*MAXN);
    hm->Hash.resize(hm->pn);
    hm->Index.resize(hm->pn);

    hRun.minX = -0.5f;
    hRun.maxX =  2.0f;
    hRun.minY = -0.5f;
    hRun.maxY =  2.0f;
    hRun.minZ = -0.5f;
    hRun.maxZ =  2.0f;

    hRun.dt = 2.0e-4 / p; //1.0e-3;
    //hRun.dt = 0.5f * hSmooth / c0;
    hRun.tsn = 4000 * p; //1000;
    hRun.ssi = 100 * p;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hm->IntDummy.resize(hm->pn);
    hm->FloatDummy.resize(hm->pn);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gx = 0.0f;
    hLoad[0].gy = -9.81f;
    hLoad[0].gz = 0.0f;

    hLoad[1].minX = dr * (20 * p);
    hLoad[1].maxX = dr * (40 * p);
    hLoad[1].minY = dr * (20 * p);
    hLoad[1].maxY = dr * (30 * p);
    hLoad[1].minZ = hRun.minZ;
    hLoad[1].maxZ = hRun.maxZ;
    hLoad[1].gx = 0.0f;
    hLoad[1].gy = 4.0f * 9.81f;
    hLoad[1].gz = 0.0f;

    printf("Pump in a box \n");
    printf("Particles: %i \n", hm->pn);
    printf("Neighbours: %i \n", hm->pn*MAXN);
    printf("Grid cells: %i \n", hGrid.nX * hGrid.nY * hGrid.nZ);
}


void initDamBreak(struct host_model *hm) {

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
            hm->PosX.push_back(i * dr + 0.5f * dr);
            hm->PosY.push_back(j * dr + 0.5f * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(m);
            pi++;
        }
    }

    for (j = -3; j <= -1; j++) {
        for (i = -3; i <= 269 * k + 2; i++) {
            hm->PosX.push_back(i * dr);
            hm->PosY.push_back(j * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(-m);
            pi++;
        }
    }

    for (j = -0; j <= 80 * k; j++) {
        for (i = -3; i <= -1; i++) {
            hm->PosX.push_back(i * dr);
            hm->PosY.push_back(j * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(-m);
            pi++;
        }
    }

    for (j = -0; j <= 80 * k; j++) {
        for (i = 269 * k; i <= 269 * k +2; i++) {
            hm->PosX.push_back(i * dr);
            hm->PosY.push_back(j * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(-m);
            pi++;
        }
    }

    hm->pn = hm->Material.size();
    if ((hm->PosX.size() != hm->pn) ||
            (hm->PosY.size() != hm->pn) ||
            (hm->PosZ.size() != hm->pn)) {
        printf("Wrong vector size\n");
        exit(-1);
    }

    hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hm->VelDotX.resize(hm->pn);
    hm->VelDotY.resize(hm->pn);
    hm->VelDotZ.resize(hm->pn);
    hm->DensityDot.resize(hm->pn);
    hm->EnergyDot.resize(hm->pn);
    hm->List.resize(hm->pn*MAXN);
    hm->Hash.resize(hm->pn);
    hm->Index.resize(hm->pn);

    hm->PosX0.resize(hm->pn);
    hm->PosY0.resize(hm->pn);
    hm->PosZ0.resize(hm->pn);
    hm->VelX0.resize(hm->pn);
    hm->VelY0.resize(hm->pn);
    hm->VelZ0.resize(hm->pn);
    hm->Density0.resize(hm->pn);
    hm->Energy0.resize(hm->pn);

    hRun.minX = -1.0f;
    hRun.maxX =  6.0f;
    hRun.minY = -1.0f;
    hRun.maxY =  4.0f;
    hRun.minZ = -1.0f;
    hRun.maxZ =  1.0f;

    hRun.dt = 4.0e-4 / k; //1.0e-3;
    hRun.tsn = 10000 * k; //1000;
    hRun.ssi = 200 * k;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hm->IntDummy.resize(hm->pn);
    hm->FloatDummy.resize(hm->pn);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;

    printf("Dam break in a box \n");
    printf("Particles: %i \n", hm->pn);
}


void initChannel(struct host_model *hm) {

    int i, j, k, m, b, q;
    double rho, c0, pmin;
    double dr;

    q = 2;

    m = 1;
    b = 2;
    rho = 1000.f;
    c0 = 20.0f;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 0;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = 4.0f;
    hMatProp[b][2] = pmin;

    dr = 0.02f / q; // x4
/*
    for (k = 0; k < 10 * q ; k++) {
		for (j = 0; j < 10 * q ; j++) {
			for (i = 0; i < 100 * q; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(m);
				hm->VelX.push_back(1.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
*/
	
    for (k = -1; k < 10 * q +1; k++) {
		for (j = -1; j < -0; j++) {
			for (i = -10; i < 100 * q + 10; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = -1; k < 0 ; k++) {
		for (j = 0; j < 15 * q; j++) {
			for (i = -10; i < 100 * q + 10; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = 10 * q; k < 10 * q +1; k++) {
		for (j = 0; j < 15 * q; j++) {
			for (i = -10; i < 100 * q + 10; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    hm->pn = hm->Material.size();
    
    //hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    //hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    //thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    //thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hRun.minX = -1.0f;
    hRun.maxX =  4.0f;
    hRun.minY = -0.5f;
    hRun.maxY =  1.0f;
    hRun.minZ = -0.5f;
    hRun.maxZ =  1.0f;

    hRun.dt = 1.0f * dr / hMatProp[m][1];
    hRun.tsn = 2000 * q; //1000;
    hRun.ssi = 40 * q;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    
    hOut[0].oX = 105.0f * q * dr;
    hOut[0].oY = 0.0f;
    hOut[0].oZ = 0.0f;
    hOut[0].nX = 1.0f;
    hOut[0].nY = 0.0f;
    hOut[0].nZ = 0.0f;
	
	hIn[0].Material = m;
	hIn[0].Mass = rho * dr * dr * dr;
	hIn[0].Smooth = 1.2f * dr;
	hIn[0].oX = 0.0f;
	hIn[0].oY = 0.0f;
	hIn[0].oZ = 0.0f;
	hIn[0].uX = 0.0f;
	hIn[0].uY = 0.0f;
	hIn[0].uZ = dr;
	hIn[0].nu = 10 * q;
	hIn[0].vX = 0.0f;
	hIn[0].vY = dr;
	hIn[0].vZ = 0.0f;
	hIn[0].nv = 10 * q;
	hIn[0].nX = 1.0f;
	hIn[0].nY = 0.0f;
	hIn[0].nZ = 0.0f;
	hIn[0].Velocity = 1.0f;
	hIn[0].Density = rho;
	hIn[0].Energy = 0.0f;
	hIn[0].Distance = 0.0f;
	
	hFix[0].minX = -0.2f;
	hFix[0].minY = -1.0f;
	hFix[0].minZ = -1.0f;
	hFix[0].maxX =  2.0 * dr;
	hFix[0].maxY =  1.0f;
	hFix[0].maxZ =  1.0f;
	hFix[0].velX =  1.0f;
	hFix[0].velY =  0.0f;
	hFix[0].velZ =  0.0f;
	
    printf("Channel \n");
    printf("Particles: %i \n", hm->pn);
    printf("Grid: %i \n", hm->SetStart.size());
}


void initBath(struct host_model *hm) {

    int i, j, k, m, b, q;
    double rho, c0, pmin;
    double dr;

    q = 3;

    m = 1;
    b = 2;
    rho = 1000.f;
    c0 = 20.0f;
    pmin = -1.e15;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 0;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = c0;
    hMatProp[b][2] = pmin;

    dr = 0.01f / q; // x4
	
    for (k = 0; k < 60 * q +1; k++) {
		for (j = 0; j < 1; j++) {
			for (i = 0; i < 60 * q + 1; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = 0; k < 60 * q +1; k++) {
		for (j = 0; j < 40 * q +1; j++) {
			for (i = 0; i < 1; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = 0; k < 60 * q +1; k++) {
		for (j = 0; j < 40 * q + 1; j++) {
			for (i = 60 * q; i < 60 * q + 1; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = 0; k < 1; k++) {
		for (j = 0; j < 40 * q + 1; j++) {
			for (i = 0; i < 60 * q + 1; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    for (k = 60 * q; k < 60 * q +1; k++) {
		for (j = 0; j < 40 * q + 1; j++) {
			for (i = 0; i < 60 * q + 1; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
				hm->Mass.push_back(rho*dr*dr*dr);
			}
		}
	}
	
    hm->pn = hm->Material.size();
    
    //hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    //hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    //thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    //thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hRun.minX = -0.1f;
    hRun.maxX =  0.7f;
    hRun.minY = -0.1f;
    hRun.maxY =  0.7f;
    hRun.minZ = -0.1f;
    hRun.maxZ =  0.7f;

    hRun.dt = 0.8f * dr / hMatProp[m][1];
    hRun.tsn = 10000 * q; //1000;
    hRun.ssi = 200 * q;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    /*
    hOut[0].oX = 105.0f * q * dr;
    hOut[0].oY = 0.0f;
    hOut[0].oZ = 0.0f;
    hOut[0].nX = 1.0f;
    hOut[0].nY = 0.0f;
    hOut[0].nZ = 0.0f;
	*/
	hIn[0].Material = m;
	hIn[0].Mass = rho * dr * dr * dr;
	hIn[0].Smooth = 1.2f * dr;
	hIn[0].oX = 0.0f;
	hIn[0].oY = 0.45f;
	hIn[0].oZ = 0.45f;
	hIn[0].uX = 0.0f;
	hIn[0].uY = 0.0f;
	hIn[0].uZ = dr;
	hIn[0].nu = 10 * q;
	hIn[0].vX = 0.0f;
	hIn[0].vY = dr;
	hIn[0].vZ = 0.0f;
	hIn[0].nv = 10 * q;
	hIn[0].nX = 1.0f;
	hIn[0].nY = 0.0f;
	hIn[0].nZ = 0.0f;
	hIn[0].Velocity = 1.0f;
	hIn[0].Density = rho;
	hIn[0].Energy = 0.0f;
	hIn[0].Distance = 0.0f;
	
	hFix[0].minX = -dr;
	hFix[0].minY = 0.4f;
	hFix[0].minZ = 0.4f;
	hFix[0].maxX =  2.0 * dr;
	hFix[0].maxY =  0.7f;
	hFix[0].maxZ =  0.7f;
	hFix[0].velX =  1.0f;
	hFix[0].velY =  0.0f;
	hFix[0].velZ =  0.0f;
	
    printf("Bath \n");
    printf("Particles: %i \n", hm->pn);
    printf("Grid: %i \n", hm->SetStart.size());
}


void initBox(struct host_model *hm) {

    int i, j, k, m, b, q;
    double rho, c0, pmin;
    double dr;

    q = 1;

    m = 1;
    b = 2;
    rho = 1000.f;
    c0 = 20.0f;
    pmin = -1.e4;

    hMatType[m] = 3;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 0;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = 1.2f*c0;
    hMatProp[b][2] = pmin;

    dr = 0.02f / q; // x4

    for (k = 0; k < 2 * q ; k++) {
		for (j = 0; j < 2 * q ; j++) {
			for (i = 0; i < 2 * q; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(m);
				hm->VelX.push_back(0.0f);
			}
		}
	}

    for (k = 0; k < 2 * q ; k++) {
		for (j = -1; j < 0; j++) {
			for (i = -0; i < 2 * q; i++) {
				hm->PosX.push_back(i * dr);
				hm->PosY.push_back(j * dr);
				hm->PosZ.push_back(k * dr);
				hm->Material.push_back(b);
				hm->VelX.push_back(0.0f);
			}
		}
	}
	
    hm->pn = hm->Material.size();
    
    hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    //hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    //thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hRun.minX = -1.0f;
    hRun.maxX =  2.0f;
    hRun.minY = -1.0f;
    hRun.maxY =  2.0f;
    hRun.minZ = -1.0f;
    hRun.maxZ =  2.0f;

    hRun.dt = dr / hMatProp[b][1];
    hRun.tsn = 50 * q; //1000;
    hRun.ssi = 1 * q;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    
    printf("Box \n");
    printf("Particles: %i \n", hm->pn);
    printf("Grid: %i \n", hm->SetStart.size());
}


void initChannel2D(struct host_model *hm) {

    int i, j, m, b, k;
    double rho, c0, pmin;
    double dr;

    k = 3;

    m = 1;
    b = 2;
    rho = 1000.f;
    c0 = 20.0f;
    pmin = -1.e12;

    hMatType[m] = 4;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 0;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = 1.2f*c0;
    hMatProp[b][2] = pmin;

    dr = 0.02f / k; // x4

    for (j = 0; j < 10 * k ; j++) {
        for (i = 0; i < 100 * k; i++) {
            hm->PosX.push_back(i * dr);
            hm->PosY.push_back(j * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(m);
            hm->VelX.push_back(1.0f);
        }
    }

    for (j = -3; j < 0; j++) {
        for (i = -10; i < 100 * k + 10; i++) {
            hm->PosX.push_back(i * dr);
            hm->PosY.push_back(j * dr);
            hm->PosZ.push_back(0.0f);
            hm->Material.push_back(b);
            hm->VelX.push_back(0.0f);
        }
    }
	
    hm->pn = hm->Material.size();
    
    hm->Mass.resize(hm->pn);
    hm->Smooth.resize(hm->pn);
    //hm->VelX.resize(hm->pn);
    hm->VelY.resize(hm->pn);
    hm->VelZ.resize(hm->pn);
    hm->Density.resize(hm->pn);
    hm->Energy.resize(hm->pn);
    hm->Pressure.resize(hm->pn);
    hm->Sound.resize(hm->pn);

    thrust::fill(hm->Mass.begin(), hm->Mass.end(), rho * dr * dr);
    thrust::fill(hm->Smooth.begin(), hm->Smooth.end(), 1.2f * dr);
    //thrust::fill(hm->VelX.begin(), hm->VelX.end(), 0.0f);
    thrust::fill(hm->VelY.begin(), hm->VelY.end(), 0.0f);
    thrust::fill(hm->VelZ.begin(), hm->VelZ.end(), 0.0f);
    thrust::fill(hm->Density.begin(), hm->Density.end(), rho);
    thrust::fill(hm->Energy.begin(), hm->Energy.end(), 0.0f);
    thrust::fill(hm->Pressure.begin(), hm->Pressure.end(), 0.0f);
    thrust::fill(hm->Sound.begin(), hm->Sound.end(), c0);

    hRun.minX = -1.0f;
    hRun.maxX =  4.0f;
    hRun.minY = -0.5f;
    hRun.maxY =  1.0f;
    hRun.minZ = -0.1f;
    hRun.maxZ =  0.1f;

    hRun.dt = dr / hMatProp[b][1];
    hRun.tsn = 50000 * k; //1000;
    hRun.ssi = 1000 * k;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
    hGrid.nX = (int) ((hRun.maxX - hRun.minX) / hGrid.size) +1;
    hGrid.nY = (int) ((hRun.maxY - hRun.minY) / hGrid.size) +1;
    hGrid.nZ = (int) ((hRun.maxZ - hRun.minZ) / hGrid.size) +1;

    hm->SetStart.resize(hGrid.nX * hGrid.nY * hGrid.nZ);
    hm->SetStop.resize(hGrid.nX * hGrid.nY * hGrid.nZ);

    hLoad[0].minX = hRun.minX;
    hLoad[0].maxX = hRun.maxX;
    hLoad[0].minY = hRun.minY;
    hLoad[0].maxY = hRun.maxY;
    hLoad[0].minZ = hRun.minZ;
    hLoad[0].maxZ = hRun.maxZ;
    hLoad[0].gy = -9.81f;
    
    hOut[0].oX = 105.0f * k * dr;
    hOut[0].oY = 0.0f;
    hOut[0].oZ = 0.0f;
    hOut[0].nX = 1.0f;
    hOut[0].nY = 0.0f;
    hOut[0].nZ = 0.0f;
	
	hIn[0].Material = m;
	hIn[0].Mass = rho * dr * dr;
	hIn[0].Smooth = 1.2f * dr;
	hIn[0].oX = 0.0f;
	hIn[0].oY = 0.0f;
	hIn[0].oZ = 0.0f;
	hIn[0].uX = 0.0f;
	hIn[0].uY = 0.0f;
	hIn[0].uZ = 1.0f;
	hIn[0].nu = 1;
	hIn[0].vX = 0.0f;
	hIn[0].vY = dr;
	hIn[0].vZ = 0.0f;
	hIn[0].nv = 10 * k;
	hIn[0].nX = 1.0f;
	hIn[0].nY = 0.0f;
	hIn[0].nZ = 0.0f;
	hIn[0].Velocity = 1.0f;
	hIn[0].Density = rho;
	hIn[0].Energy = 0.0f;
	hIn[0].Distance = 0.0f;
	
	hFix[0].minX = -0.2f;
	hFix[0].minY = -1.0f;
	hFix[0].minZ = -1.0f;
	hFix[0].maxX =  2.0 * dr;
	hFix[0].maxY =  1.0f;
	hFix[0].maxZ =  1.0f;
	hFix[0].velX =  1.0f;
	hFix[0].velY =  0.0f;
	hFix[0].velZ =  0.0f;
	
    printf("Channel \n");
    printf("Particles: %i \n", hm->pn);
    printf("Grid: %i \n", hm->SetStart.size());
}


int copyHostToDevice(struct host_model *hm, struct device_model *dm) {

    dm->pn = hm->pn;
    dm->Material = hm->Material;
    dm->Mass = hm->Mass;
    dm->Smooth = hm->Smooth;
    dm->PosX = hm->PosX;
    dm->PosY = hm->PosY;
    dm->PosZ = hm->PosZ;
    dm->VelX = hm->VelX;
    dm->VelY = hm->VelY;
    dm->VelZ = hm->VelZ;
    dm->Density = hm->Density;
    dm->Energy = hm->Energy;
    dm->Pressure = hm->Pressure;
	
    dm->Sound = hm->Sound;
    dm->VelDotX = hm->VelDotX;
    dm->VelDotY = hm->VelDotY;
    dm->VelDotZ = hm->VelDotZ;
    dm->DensityDot = hm->DensityDot;
    dm->EnergyDot = hm->EnergyDot;

    dm->PosX0 = hm->PosX0;
    dm->PosY0 = hm->PosY0;
    dm->PosZ0 = hm->PosZ0;
    dm->VelX0 = hm->VelX0;
    dm->VelY0 = hm->VelY0;
    dm->VelZ0 = hm->VelZ0;
    dm->Density0 = hm->Density0;
    dm->Energy0 = hm->Energy0;

    dm->List = hm->List;
    dm->Hash = hm->Hash;
    dm->Index = hm->Index;

    dm->SetStart = hm->SetStart;
    dm->SetStop = hm->SetStop;

    dm->IntDummy = hm->IntDummy;
    dm->FloatDummy = hm->FloatDummy;

    dGrid.oX = hGrid.oX;
    dGrid.oY = hGrid.oY;
    dGrid.oZ = hGrid.oZ;
    dGrid.nX = hGrid.nX;
    dGrid.nY = hGrid.nY;
    dGrid.nZ = hGrid.nZ;
    dGrid.size = hGrid.size;
    
    for (int i = 0; i < 10; i++) dIn[i] = hIn[i];

    return 0;
}


int copyDeviceToHost(struct device_model *dm, struct host_model *hm) {

    hm->pn = dm->pn;

    thrust::copy(dm->Material.begin(), dm->Material.end(), hm->Material.begin());
    thrust::copy(dm->Mass.begin(), dm->Mass.end(), hm->Mass.begin());
    thrust::copy(dm->Smooth.begin(), dm->Smooth.end(), hm->Smooth.begin());
    thrust::copy(dm->PosX.begin(), dm->PosX.end(), hm->PosX.begin());
    thrust::copy(dm->PosY.begin(), dm->PosY.end(), hm->PosY.begin());
    thrust::copy(dm->PosZ.begin(), dm->PosZ.end(), hm->PosZ.begin());
    thrust::copy(dm->VelX.begin(), dm->VelX.end(), hm->VelX.begin());
    thrust::copy(dm->VelY.begin(), dm->VelY.end(), hm->VelY.begin());
    thrust::copy(dm->VelZ.begin(), dm->VelZ.end(), hm->VelZ.begin());
    thrust::copy(dm->Density.begin(), dm->Density.end(), hm->Density.begin());
    thrust::copy(dm->Energy.begin(), dm->Energy.end(), hm->Energy.begin());
    thrust::copy(dm->Pressure.begin(), dm->Pressure.end(), hm->Pressure.begin());
    thrust::copy(dm->Sound.begin(), dm->Sound.end(), hm->Sound.begin());
    
    thrust::copy(dm->VelDotX.begin(), dm->VelDotX.end(), hm->VelDotX.begin());
    thrust::copy(dm->VelDotY.begin(), dm->VelDotY.end(), hm->VelDotY.begin());
    thrust::copy(dm->VelDotZ.begin(), dm->VelDotZ.end(), hm->VelDotZ.begin());
    thrust::copy(dm->DensityDot.begin(), dm->DensityDot.end(), hm->DensityDot.begin());
    thrust::copy(dm->EnergyDot.begin(), dm->EnergyDot.end(), hm->EnergyDot.begin());

    thrust::copy(dm->List.begin(), dm->List.end(), hm->List.begin());
    thrust::copy(dm->Hash.begin(), dm->Hash.end(), hm->Hash.begin());
    thrust::copy(dm->Index.begin(), dm->Index.end(), hm->Index.begin());

    thrust::copy(dm->SetStart.begin(), dm->SetStart.end(), hm->SetStart.begin());
    thrust::copy(dm->SetStop.begin(), dm->SetStop.end(), hm->SetStop.begin());

    thrust::copy(dm->IntDummy.begin(), dm->IntDummy.end(), hm->IntDummy.begin());
    thrust::copy(dm->FloatDummy.begin(), dm->FloatDummy.end(), hm->FloatDummy.begin());

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


int copyHostToPointer(struct host_model *hm, struct pointer_model *pm) {

    pm->pn = hm->pn;

    pm->Material = thrust::raw_pointer_cast(&hm->Material[0]);
    pm->Mass = thrust::raw_pointer_cast(&hm->Mass[0]);
    pm->Smooth = thrust::raw_pointer_cast(&hm->Smooth[0]);
    pm->PosX = thrust::raw_pointer_cast(&hm->PosX[0]);
    pm->PosY = thrust::raw_pointer_cast(&hm->PosY[0]);
    pm->PosZ = thrust::raw_pointer_cast(&hm->PosZ[0]);
    pm->VelX = thrust::raw_pointer_cast(&hm->VelX[0]);
    pm->VelY = thrust::raw_pointer_cast(&hm->VelY[0]);
    pm->VelZ = thrust::raw_pointer_cast(&hm->VelZ[0]);
    pm->Density = thrust::raw_pointer_cast(&hm->Density[0]);
    pm->Energy = thrust::raw_pointer_cast(&hm->Energy[0]);
    pm->Pressure = thrust::raw_pointer_cast(&hm->Pressure[0]);
    pm->Sound = thrust::raw_pointer_cast(&hm->Sound[0]);
    
    pm->VelDotX = thrust::raw_pointer_cast(&hm->VelDotX[0]);
    pm->VelDotY = thrust::raw_pointer_cast(&hm->VelDotY[0]);
    pm->VelDotZ = thrust::raw_pointer_cast(&hm->VelDotZ[0]);
    pm->DensityDot = thrust::raw_pointer_cast(&hm->DensityDot[0]);
    pm->EnergyDot = thrust::raw_pointer_cast(&hm->EnergyDot[0]);
    pm->PosX0 = thrust::raw_pointer_cast(&hm->PosX0[0]);
    pm->PosY0 = thrust::raw_pointer_cast(&hm->PosY0[0]);
    pm->PosZ0 = thrust::raw_pointer_cast(&hm->PosZ0[0]);
    pm->VelX0 = thrust::raw_pointer_cast(&hm->VelX0[0]);
    pm->VelY0 = thrust::raw_pointer_cast(&hm->VelY0[0]);
    pm->VelZ0 = thrust::raw_pointer_cast(&hm->VelZ0[0]);
    pm->Density0 = thrust::raw_pointer_cast(&hm->Density0[0]);
    pm->Energy0 = thrust::raw_pointer_cast(&hm->Energy0[0]);
    pm->Hash = thrust::raw_pointer_cast(&hm->Hash[0]);
    pm->Index = thrust::raw_pointer_cast(&hm->Index[0]);
    pm->SetStart = thrust::raw_pointer_cast(&hm->SetStart[0]);
    pm->SetStop = thrust::raw_pointer_cast(&hm->SetStop[0]);
    pm->List = thrust::raw_pointer_cast(&hm->List[0]);
    pm->IntDummy = thrust::raw_pointer_cast(&hm->IntDummy[0]);
    pm->FloatDummy = thrust::raw_pointer_cast(&hm->FloatDummy[0]);

    return 0;
}

int copyDeviceToPointer(struct device_model *dm, struct pointer_model *pm) {

    pm->pn = dm->pn;

    pm->Material = thrust::raw_pointer_cast(&dm->Material[0]);
    pm->Mass = thrust::raw_pointer_cast(&dm->Mass[0]);
    pm->Smooth = thrust::raw_pointer_cast(&dm->Smooth[0]);
    pm->PosX = thrust::raw_pointer_cast(&dm->PosX[0]);
    pm->PosY = thrust::raw_pointer_cast(&dm->PosY[0]);
    pm->PosZ = thrust::raw_pointer_cast(&dm->PosZ[0]);
    pm->VelX = thrust::raw_pointer_cast(&dm->VelX[0]);
    pm->VelY = thrust::raw_pointer_cast(&dm->VelY[0]);
    pm->VelZ = thrust::raw_pointer_cast(&dm->VelZ[0]);
    pm->Density = thrust::raw_pointer_cast(&dm->Density[0]);
    pm->Energy = thrust::raw_pointer_cast(&dm->Energy[0]);
    pm->Pressure = thrust::raw_pointer_cast(&dm->Pressure[0]);
    pm->Sound = thrust::raw_pointer_cast(&dm->Sound[0]);
    
    pm->VelDotX = thrust::raw_pointer_cast(&dm->VelDotX[0]);
    pm->VelDotY = thrust::raw_pointer_cast(&dm->VelDotY[0]);
    pm->VelDotZ = thrust::raw_pointer_cast(&dm->VelDotZ[0]);
    pm->DensityDot = thrust::raw_pointer_cast(&dm->DensityDot[0]);
    pm->EnergyDot = thrust::raw_pointer_cast(&dm->EnergyDot[0]);
    pm->PosX0 = thrust::raw_pointer_cast(&dm->PosX0[0]);
    pm->PosY0 = thrust::raw_pointer_cast(&dm->PosY0[0]);
    pm->PosZ0 = thrust::raw_pointer_cast(&dm->PosZ0[0]);
    pm->VelX0 = thrust::raw_pointer_cast(&dm->VelX0[0]);
    pm->VelY0 = thrust::raw_pointer_cast(&dm->VelY0[0]);
    pm->VelZ0 = thrust::raw_pointer_cast(&dm->VelZ0[0]);
    pm->Density0 = thrust::raw_pointer_cast(&dm->Density0[0]);
    pm->Energy0 = thrust::raw_pointer_cast(&dm->Energy0[0]);
    pm->Hash = thrust::raw_pointer_cast(&dm->Hash[0]);
    pm->Index = thrust::raw_pointer_cast(&dm->Index[0]);
    pm->SetStart = thrust::raw_pointer_cast(&dm->SetStart[0]);
    pm->SetStop = thrust::raw_pointer_cast(&dm->SetStop[0]);
    pm->List = thrust::raw_pointer_cast(&dm->List[0]);
    pm->IntDummy = thrust::raw_pointer_cast(&dm->IntDummy[0]);
    pm->FloatDummy = thrust::raw_pointer_cast(&dm->FloatDummy[0]);

    return 0;
}

int resizeHost(struct host_model *hm) {

    hm->Material.resize(hm->pn, 0);
    hm->Mass.resize(hm->pn), 0.0f;
    hm->Smooth.resize(hm->pn, 0.0f);
    hm->PosX.resize(hm->pn, 0.0f);
    hm->PosY.resize(hm->pn, 0.0f);
    hm->PosZ.resize(hm->pn, 0.0f);
    hm->VelX.resize(hm->pn, 0.0f);
    hm->VelY.resize(hm->pn, 0.0f);
    hm->VelZ.resize(hm->pn, 0.0f);
    hm->Density.resize(hm->pn, 0.0f);
    hm->Energy.resize(hm->pn, 0.0f);
    hm->Pressure.resize(hm->pn, 0.0f);
    hm->Sound.resize(hm->pn, 0.0f);
    hm->VelDotX.resize(hm->pn, 0.0f);
    hm->VelDotY.resize(hm->pn, 0.0f);
    hm->VelDotZ.resize(hm->pn, 0.0f);
    hm->DensityDot.resize(hm->pn, 0.0f);
    hm->EnergyDot.resize(hm->pn, 0.0f);
    hm->PosX0.resize(hm->pn, 0.0f);
    hm->PosY0.resize(hm->pn, 0.0f);
    hm->PosZ0.resize(hm->pn, 0.0f);
    hm->VelX0.resize(hm->pn, 0.0f);
    hm->VelY0.resize(hm->pn, 0.0f);
    hm->VelZ0.resize(hm->pn, 0.0f);
    hm->Density0.resize(hm->pn, 0.0f);
    hm->Energy0.resize(hm->pn, 0.0f);
    hm->Hash.resize(hm->pn, 0);
    hm->Index.resize(hm->pn, 0);
    hm->List.resize(hm->pn * MAXN, 0);
    hm->IntDummy.resize(hm->pn, 0);
    hm->FloatDummy.resize(hm->pn, 0.0f);

    return 0;
}


int resizeDevice(struct device_model *dm) {
	
    dm->Material.resize(dm->pn, 0);
    dm->Mass.resize(dm->pn, 0.0f);
    dm->Smooth.resize(dm->pn, 0.0f);
    dm->PosX.resize(dm->pn, 0.0f);
    dm->PosY.resize(dm->pn, 0.0f);
    dm->PosZ.resize(dm->pn, 0.0f);
    dm->VelX.resize(dm->pn, 0.0f);
    dm->VelY.resize(dm->pn, 0.0f);
    dm->VelZ.resize(dm->pn, 0.0f);
    dm->Density.resize(dm->pn, 0.0f);
    dm->Energy.resize(dm->pn, 0.0f);
    dm->Pressure.resize(dm->pn, 0.0f);
    dm->Sound.resize(dm->pn, 0.0f);
    
    dm->VelDotX.resize(dm->pn, 0.0f);
    dm->VelDotY.resize(dm->pn, 0.0f);
    dm->VelDotZ.resize(dm->pn, 0.0f);
    dm->DensityDot.resize(dm->pn, 0.0f);
    dm->EnergyDot.resize(dm->pn, 0.0f);
    dm->PosX0.resize(dm->pn, 0.0f);
    dm->PosY0.resize(dm->pn, 0.0f);
    dm->PosZ0.resize(dm->pn, 0.0f);
    dm->VelX0.resize(dm->pn, 0.0f);
    dm->VelY0.resize(dm->pn, 0.0f);
    dm->VelZ0.resize(dm->pn, 0.0f);
    dm->Density0.resize(dm->pn, 0.0f);
    dm->Energy0.resize(dm->pn, 0.0f);
    dm->Hash.resize(dm->pn, 0);
    dm->Index.resize(dm->pn, 0);
    dm->List.resize(dm->pn * MAXN, 0);
    dm->IntDummy.resize(dm->pn, 0);
    dm->FloatDummy.resize(dm->pn, 0.0f);
	
    return 0;
}

int printData(struct host_model *hm) {
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

    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%d %d %d\n", i, hm->Index[i], hm->Hash[i]);
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


int outputVTK(struct host_model *hm, int ss) {
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


__host__ void addInletHost(const int pn, const int* Index, const float* Smooth, 
                           float* PosX, float* PosY, float* PosZ,
                           float* PosX0, float* PosY0, float* PosZ0) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

	int ip;
	float dx, dy, dz, d;
	
    for (ip = 0; ip < pn; ip++) 
		if (Index[ip] == 1) {
			dx = PosX[ip] - PosX0[ip];
			dy = PosY[ip] - PosY0[ip];
			dz = PosZ[ip] - PosZ0[ip];
			d = sqrtf(dx*dx + dy*dy + dz*dz);
			
			dx *= Smooth[ip] / (1.2f * d);
			dy *= Smooth[ip] / (1.2f * d);
			dz *= Smooth[ip] / (1.2f * d);
			
			PosX[ip] -= dx;
			PosY[ip] -= dy;
			PosZ[ip] -= dz;
			PosX0[ip] = PosX[ip];
			PosY0[ip] = PosY[ip];
			PosZ0[ip] = PosZ[ip];
		}
}


__host__ void updateHashHost(const int pn, const struct grid Grid,
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

        Hash[ip] = ic;
    }
}

__host__ void checkBoundariesHost(const int pn, const struct grid Grid,
                                  const float* PosX, const float* PosY, const float* PosZ,
                                  int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip;
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
	
	minX = Grid.oX; maxX = Grid.oX + Grid.size * Grid.nX;
	minY = Grid.oY; maxY = Grid.oY + Grid.size * Grid.nY;
	minZ = Grid.oZ; maxZ = Grid.oZ + Grid.size * Grid.nZ;
	
    for (ip = 0; ip < pn; ip++) {
        if ((PosX[ip] < minX) || (PosX[ip] > maxX) ||
            (PosY[ip] < minY) || (PosY[ip] > maxY) ||
            (PosZ[ip] < minZ) || (PosZ[ip] > maxZ)) {
            Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
        }
    }
}

__global__ void checkBoundariesDevice(const int pn, const struct grid Grid,
                                      float* PosX, float* PosY, float* PosZ,
                                      int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    int ip;
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
	
    ip = threadIdx.x + blockDim.x * blockIdx.x;
    
	minX = Grid.oX; maxX = Grid.oX + Grid.size * Grid.nX;
	minY = Grid.oY; maxY = Grid.oY + Grid.size * Grid.nY;
	minZ = Grid.oZ; maxZ = Grid.oZ + Grid.size * Grid.nZ;
	
    if (ip < pn) {
		
        if ((PosX[ip] < minX) || (PosX[ip] > maxX) ||
            (PosY[ip] < minY) || (PosY[ip] > maxY) ||
            (PosZ[ip] < minZ) || (PosZ[ip] > maxZ)) {
			//PosX[ip] = 0.5f*(minX + maxX);
			//PosY[ip] = 0.5f*(minY + maxY);
			//PosZ[ip] = 0.5f*(minZ + maxZ);
			Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
        }
        /*
        if (PosX[ip] < minX) PosX[ip] = minX;
        if (PosY[ip] < minY) PosY[ip] = minY;
        if (PosZ[ip] < minZ) PosZ[ip] = minZ;
        if (PosX[ip] > maxX) PosX[ip] = maxX;
        if (PosY[ip] > maxY) PosY[ip] = maxY;
        if (PosZ[ip] > maxZ) PosZ[ip] = maxZ;
        if (Hash[ip] > Grid.nX * Grid.nY * Grid.nZ) Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
        */
    }
}

__host__ void checkOutletHost(const int pn, const struct grid Grid, 
                              const float* PosX, const float* PosY, const float* PosZ,
                              int* Hash) {
    
    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    float d;
    int ip, i;
    
    for (ip = 0; ip < pn; ip++) {
		for (i = 0; i < 10; i++) {
	        d = 0.0f;
			d += (PosX[ip] - hOut[i].oX) * hOut[i].nX;
			d += (PosY[ip] - hOut[i].oY) * hOut[i].nY;
			d += (PosZ[ip] - hOut[i].oZ) * hOut[i].nZ;
			if (d > 0.0f) {
				Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
			}
		}
	}
}

__global__ void checkOutletDevice(const int pn, const struct grid Grid, 
	const float* PosX, const float* PosY, const float* PosZ,
    int* Hash) {

    /**
     * \brief Update particles
     *
     * \date Jan 6, 2010
     * \author Luca Massidda
     */

    float d;
    int ip, i;

    ip = threadIdx.x + blockDim.x * blockIdx.x;

    if (ip < pn) {
		for (i = 0; i < 10; i++) {
	        d = 0.0f;
			d += (PosX[ip] - dOut[i].oX) * dOut[i].nX;
			d += (PosY[ip] - dOut[i].oY) * dOut[i].nY;
			d += (PosZ[ip] - dOut[i].oZ) * dOut[i].nZ;
			if (d > 0.0f) {
				Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
			}
		}
	}
}


__host__ void makeInletHost(struct host_model *hm, struct pointer_model *pm) {
	
	int i, j, iu, iv;
	
    for (i = 0; i < 10; i++) if (hIn[i].Material != 0) {
		hIn[i].Distance += hIn[i].Velocity * hRun.dt;
		
		if ((1.2f * hIn[i].Distance) > hIn[i].Smooth) {
			hIn[i].Distance = 0.0f;
			
			for (iv = 0; iv < hIn[i].nv; iv++) {
				for (iu = 0; iu < hIn[i].nu; iu++) {
					hm->pn++;
					pm->pn++;
					
					hm->Material.push_back(hIn[i].Material);
					hm->Mass.push_back(hIn[i].Mass);
					hm->Smooth.push_back(hIn[i].Smooth);
					
					hm->PosX.push_back(hIn[i].oX + iu * hIn[i].uX + iv * hIn[i].vX);
					hm->PosY.push_back(hIn[i].oY + iu * hIn[i].uY + iv * hIn[i].vY);
					hm->PosZ.push_back(hIn[i].oZ + iu * hIn[i].uZ + iv * hIn[i].vZ);
					
					hm->VelX.push_back(hIn[i].Velocity * hIn[i].nX);
					hm->VelY.push_back(hIn[i].Velocity * hIn[i].nY);
					hm->VelZ.push_back(hIn[i].Velocity * hIn[i].nZ);
					
					hm->Density.push_back(hIn[i].Density);
					hm->Energy.push_back(hIn[i].Energy);
					
					hm->PosX0.push_back(hIn[i].oX + iu * hIn[i].uX + iv * hIn[i].vX);
					hm->PosY0.push_back(hIn[i].oY + iu * hIn[i].uY + iv * hIn[i].vY);
					hm->PosZ0.push_back(hIn[i].oZ + iu * hIn[i].uZ + iv * hIn[i].vZ);
					
					hm->VelX0.push_back(hIn[i].Velocity * hIn[i].nX);
					hm->VelY0.push_back(hIn[i].Velocity * hIn[i].nY);
					hm->VelZ0.push_back(hIn[i].Velocity * hIn[i].nZ);
					
					hm->Density0.push_back(hIn[i].Density);
					hm->Energy0.push_back(hIn[i].Energy);
					
					hm->VelDotX.push_back(0.0f);
					hm->VelDotY.push_back(0.0f);
					hm->VelDotZ.push_back(0.0f);
		
					hm->DensityDot.push_back(0.0f);
					hm->EnergyDot.push_back(0.0f);
					
					hm->Pressure.push_back(0.0f);
					hm->Sound.push_back(0.0f);
					
					for (j = 0; j < MAXN; j++) hm->List.push_back(0);
					hm->Hash.push_back(0);
					hm->Index.push_back(0);
					hm->IntDummy.push_back(0);
					hm->FloatDummy.push_back(0.0f);
				}
			}
		}
	}
	
	copyHostToPointer(hm, pm);
}


void makeInletDevice(struct device_model *dm, struct pointer_model *pm) {
	
	int i, j, iu, iv;
	
    for (i = 0; i < 10; i++) if (dIn[i].Material != 0) {
		dIn[i].Distance += dIn[i].Velocity * hRun.dt;
		
		if ((1.2f * dIn[i].Distance) > dIn[i].Smooth) {
			dIn[i].Distance = 0.0f;
			
			for (iv = 0; iv < dIn[i].nv; iv++) {
				for (iu = 0; iu < dIn[i].nu; iu++) {
					dm->pn++;
					pm->pn++;
					
					dm->Material.push_back(dIn[i].Material);
					dm->Mass.push_back(dIn[i].Mass);
					dm->Smooth.push_back(dIn[i].Smooth);
					
					dm->PosX.push_back(dIn[i].oX + iu * dIn[i].uX + iv * dIn[i].vX);
					dm->PosY.push_back(dIn[i].oY + iu * dIn[i].uY + iv * dIn[i].vY);
					dm->PosZ.push_back(dIn[i].oZ + iu * dIn[i].uZ + iv * dIn[i].vZ);
					
					dm->VelX.push_back(dIn[i].Velocity * dIn[i].nX);
					dm->VelY.push_back(dIn[i].Velocity * dIn[i].nY);
					dm->VelZ.push_back(dIn[i].Velocity * dIn[i].nZ);
					
					dm->Density.push_back(dIn[i].Density);
					dm->Energy.push_back(dIn[i].Energy);
					
					dm->PosX0.push_back(dIn[i].oX + iu * dIn[i].uX + iv * dIn[i].vX);
					dm->PosY0.push_back(dIn[i].oY + iu * dIn[i].uY + iv * dIn[i].vY);
					dm->PosZ0.push_back(dIn[i].oZ + iu * dIn[i].uZ + iv * dIn[i].vZ);
					
					dm->VelX0.push_back(dIn[i].Velocity * dIn[i].nX);
					dm->VelY0.push_back(dIn[i].Velocity * dIn[i].nY);
					dm->VelZ0.push_back(dIn[i].Velocity * dIn[i].nZ);
					
					dm->Density0.push_back(dIn[i].Density);
					dm->Energy0.push_back(dIn[i].Energy);
					
					dm->VelDotX.push_back(0.0f);
					dm->VelDotY.push_back(0.0f);
					dm->VelDotZ.push_back(0.0f);
		
					dm->DensityDot.push_back(0.0f);
					dm->EnergyDot.push_back(0.0f);
					
					dm->Pressure.push_back(0.0f);
					dm->Sound.push_back(0.0f);
					
					for (j = 0; j < MAXN; j++) dm->List.push_back(0);
					dm->Hash.push_back(0);
					dm->Index.push_back(0);
					dm->IntDummy.push_back(0);
					dm->FloatDummy.push_back(0.0f);
				}
			}
		}
	}
	
	copyDeviceToPointer(dm, pm);
}


__host__ void updateSetsHost(const int pn, int *SetStart, int *SetStop,
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


__host__ void updateListHost(const int pn, int *List,
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

    while (np < MAXN) {
        List[ip * MAXN + np] = ip;
        np++;
    }

}

struct is_in
{
    __host__ __device__
    bool operator()(int x)
    {
        return x == 1;
    }
};

struct is_out
{
    __host__ __device__
    bool operator()(int x)
    {
        return x == -1;
    }
};

int neighbourListDevice(struct device_model *dm, struct pointer_model *pm) {
	int pout;

    int blocks, threads;

    blocks = (pm->pn + THREADS - 1) / THREADS;
    threads = THREADS;
    
	try
	{
		updateHashDevice <<< blocks, threads >>>
		(pm->pn, dGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in updateHash: %s", e.what());
	}
	
    try
    {
		checkOutletDevice <<< blocks, threads >>>
		(pm->pn, dGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in checkOutlet: %s", e.what());
	}
	
    try
    {
		checkBoundariesDevice <<< blocks, threads >>>
		(pm->pn, dGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in checkBoundaries: %s", e.what());
	}
	
	
    try
    {
		thrust::sequence(dm->Index.begin(), dm->Index.end(), 1);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in sequence: %s", e.what());
	}
	
    try
    {
		thrust::sort_by_key(dm->Hash.begin(), dm->Hash.end(), dm->Index.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in sort: %s", e.what());
	}
    
    try
    {
		thrust::copy(dm->Material.begin(), dm->Material.end(), dm->IntDummy.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in copy: %s\n", e.what());
		printf("material: %d %d\n", dm->Material.size(), dm->Material.end() - dm->Material.begin());
		printf("dummy: %d %d\n", dm->IntDummy.size(), dm->IntDummy.end() - dm->IntDummy.begin());
	}
    try
    {
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->IntDummy.begin(), dm->Material.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in gather: %s", e.what());
	}
		
    try
    {
		thrust::copy(dm->Mass.begin(), dm->Mass.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->Mass.begin());
		
		thrust::copy(dm->Smooth.begin(), dm->Smooth.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->Smooth.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in copy-gather 1: %s", e.what());
	}
		
    try
    {
		thrust::copy(dm->PosX.begin(), dm->PosX.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->PosX.begin());
		
		thrust::copy(dm->PosY.begin(), dm->PosY.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->PosY.begin());
		
		thrust::copy(dm->PosZ.begin(), dm->PosZ.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->PosZ.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in copy-gather 2: %s", e.what());
	}
		
    try
    {
		thrust::copy(dm->VelX.begin(), dm->VelX.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->VelX.begin());
		
		thrust::copy(dm->VelY.begin(), dm->VelY.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->VelY.begin());
		
		thrust::copy(dm->VelZ.begin(), dm->VelZ.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->VelZ.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in copy-gather 3: %s", e.what());
	}
		
    try
    {
		thrust::copy(dm->Density.begin(), dm->Density.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->Density.begin());
		
		thrust::copy(dm->Energy.begin(), dm->Energy.end(), dm->FloatDummy.begin());
		thrust::gather(dm->Index.begin(), dm->Index.end(), dm->FloatDummy.begin(), dm->Energy.begin());
	}
	catch(thrust::system_error &e)
	{
		printf("Error in copy-gather 4: %s", e.what());
	}
	
    try
    {
		pout = thrust::count(dm->Hash.begin(), dm->Hash.end(), (dGrid.nX * dGrid.nY * dGrid.nZ));
		dm->pn -= pout;
		pm->pn -= pout;
		resizeDevice(dm);
		copyDeviceToPointer(dm, pm);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in resize: %s", e.what());
	}
	
    try
    {
		thrust::fill(dm->SetStart.begin(), dm->SetStart.end(), 0);
		thrust::fill(dm->SetStop.begin(), dm->SetStop.end(), 0);
		
		updateSetsDevice <<< blocks, threads >>>
		(pm->pn, pm->SetStart, pm->SetStop, pm->Hash);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in updateSets: %s", e.what());
	}
	
    try
    {
		updateListDevice <<< blocks, threads >>>
		(pm->pn, pm->List, pm->SetStart, pm->SetStop, dGrid, pm->Smooth,
		pm->PosX, pm->PosY, pm->PosZ);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in updateList: %s", e.what());
	}
	
    return 0;
}

int neighbourListHost(struct host_model *hm, struct pointer_model *pm) {
	int pout;
	
    updateHashHost(pm->pn, hGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	
	checkOutletHost(pm->pn, hGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	
	checkBoundariesHost(pm->pn, hGrid, pm->PosX, pm->PosY, pm->PosZ, pm->Hash);
	
    try {
        thrust::sequence(hm->Index.begin(), hm->Index.end());
        thrust::sort_by_key(hm->Hash.begin(), hm->Hash.end(), hm->Index.begin());
    } catch(std::bad_alloc &e) {
        printf("Ran out of memory while sorting: %s\n", e.what());
        exit(-1);
    } catch(thrust::system_error &e) {
        printf("Error accessing vector element: %s \n", e.what());
        exit(-1);
    }
	
	thrust::copy(hm->Material.begin(), hm->Material.end(), hm->IntDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->IntDummy.begin(), hm->Material.begin());
	
	thrust::copy(hm->Mass.begin(), hm->Mass.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->Mass.begin());
	
	thrust::copy(hm->Smooth.begin(), hm->Smooth.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->Smooth.begin());
	
	thrust::copy(hm->PosX.begin(), hm->PosX.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->PosX.begin());
	
	thrust::copy(hm->PosY.begin(), hm->PosY.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->PosY.begin());
	
	thrust::copy(hm->PosZ.begin(), hm->PosZ.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->PosZ.begin());
	
	thrust::copy(hm->VelX.begin(), hm->VelX.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->VelX.begin());
	
	thrust::copy(hm->VelY.begin(), hm->VelY.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->VelY.begin());
	
	thrust::copy(hm->VelZ.begin(), hm->VelZ.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->VelZ.begin());
	
	thrust::copy(hm->Density.begin(), hm->Density.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->Density.begin());
	
	thrust::copy(hm->Energy.begin(), hm->Energy.end(), hm->FloatDummy.begin());
	thrust::gather(hm->Index.begin(), hm->Index.end(), hm->FloatDummy.begin(), hm->Energy.begin());
	
	pout = thrust::count(hm->Hash.begin(), hm->Hash.end(), (hGrid.nX * hGrid.nY * hGrid.nZ));
    hm->pn -= pout;
    pm->pn -= pout;
	resizeHost(hm);
	copyHostToPointer(hm, pm);
	
    thrust::fill(hm->SetStart.begin(), hm->SetStart.end(), 0);
    thrust::fill(hm->SetStop.begin(), hm->SetStop.end(), 0);
	
    updateSetsHost(pm->pn, pm->SetStart, pm->SetStop, pm->Hash);

    updateListHost(pm->pn, pm->List, pm->SetStart, pm->SetStop, hGrid, pm->Smooth,
                   pm->PosX, pm->PosY, pm->PosZ);

    return 0;
}

int backupDataHost(struct host_model *hm) {

    thrust::copy(hm->PosX.begin(), hm->PosX.end(), hm->PosX0.begin());
    thrust::copy(hm->PosY.begin(), hm->PosY.end(), hm->PosY0.begin());
    thrust::copy(hm->PosZ.begin(), hm->PosZ.end(), hm->PosZ0.begin());
    thrust::copy(hm->VelX.begin(), hm->VelX.end(), hm->VelX0.begin());
    thrust::copy(hm->VelY.begin(), hm->VelY.end(), hm->VelY0.begin());
    thrust::copy(hm->VelZ.begin(), hm->VelZ.end(), hm->VelZ0.begin());
    thrust::copy(hm->Density.begin(), hm->Density.end(), hm->Density0.begin());
    thrust::copy(hm->Energy.begin(), hm->Energy.end(), hm->Energy0.begin());

    return 0;
}

int backupDataDevice(struct device_model *dm) {

    thrust::copy(dm->PosX.begin(), dm->PosX.end(), dm->PosX0.begin());
    thrust::copy(dm->PosY.begin(), dm->PosY.end(), dm->PosY0.begin());
    thrust::copy(dm->PosZ.begin(), dm->PosZ.end(), dm->PosZ0.begin());
    thrust::copy(dm->VelX.begin(), dm->VelX.end(), dm->VelX0.begin());
    thrust::copy(dm->VelY.begin(), dm->VelY.end(), dm->VelY0.begin());
    thrust::copy(dm->VelZ.begin(), dm->VelZ.end(), dm->VelZ0.begin());
    thrust::copy(dm->Density.begin(), dm->Density.end(), dm->Density0.begin());
    thrust::copy(dm->Energy.begin(), dm->Energy.end(), dm->Energy0.begin());

    return 0;
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


__host__ void updateParticlesHost(const int pn, const float alpha,
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
            //Pressure[ip] = pressureShock(dMatProp[iMaterial], iDensity, iEnergy);
            //Sound[ip] = soundShock(dMatProp[iMaterial], iDensity, iEnergy);
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


__host__ __device__ float kernelGauss(float r, float h) {

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


    alpha = 1.0 / (pow(h, 1) * pow(3.14, 0.5));
    //alpha = 1.0f / (3.14f * h2);

    w = 0.0f;
    //dwdr = 0.0;

    if (q2 < 4.0f) {
        w = alpha * expf(-q2);
        //dwdr = w * (-2.0 * r / h2);
    }

    return w;
}


__host__ __device__ float kernelDerivGauss(float r, float h) {

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


    alpha = 1.0f / (h * powf(3.14f, 0.5f));
    //alpha = 1.0f / (3.14f * h2);

    w = 0.0f;
    dwdr = 0.0f;

    if (q2 < 4.0f) {
        w = alpha * expf(-q2);
        dwdr = w * (-2.0f * r / h2);
    }

    return dwdr;
}


__host__ __device__ float kernelSpiky(float r, float h) {

    float q, alpha, w;
    /**
     * \brief Spiky kernel
     *
     * \date Dec 21, 2010
     * \author Luca Massidda
     */

    q = r / h;

    alpha = 15.0f / (64.0f * 3.14f * pow(h, 3));

    w = 0.0f;

    if (q < 2.0f) {
        w = alpha * powf(2.0f - q, 3);
    }

    return w;
}


__host__ __device__ float kernelDerivSpiky(float r, float h) {

    float q, alpha, dwdr;
    /**
     * \brief Gauss kernel
     *
     * \date Dec 21, 2010
     * \author Luca Massidda
     */

    q = r / h;


    alpha = -45.0f / (64.0f * 3.14f * pow(h, 4));

    dwdr = 0.0;

    if (q < 2.0f) {
        dwdr = alpha * powf(2.0f - q, 2);
    }

    return dwdr;
}


__host__ void updateLoadsHost(const int pn,
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

__host__ void balanceMassMomentumHost(const int pn, const int* List,
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
//        for (jp = 0; jp < pn; jp++) {
            
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
				f = 0.02f * w / w0 * Sound[ip] * Sound[jp] / dr;
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

__host__ void balanceEnergyHost(const int pn,
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

int RKstepHost(struct host_model *hm, struct pointer_model *pm, float alpha) {
	
    thrust::fill(hm->VelDotX.begin(), hm->VelDotX.end(), 0.0f);
    thrust::fill(hm->VelDotY.begin(), hm->VelDotY.end(), 0.0f);
    thrust::fill(hm->VelDotZ.begin(), hm->VelDotZ.end(), 0.0f);
    thrust::fill(hm->DensityDot.begin(), hm->DensityDot.end(), 0.0f);
    thrust::fill(hm->EnergyDot.begin(), hm->EnergyDot.end(), 0.0f);

    // External loads
    updateLoadsHost(pm->pn, pm->Material, 
                    pm->PosX, pm->PosY, pm->PosZ,
                    pm->VelX, pm->VelY, pm->VelZ, 
                    pm->VelDotX, pm->VelDotY, pm->VelDotZ, pm->EnergyDot);

    // Calculate particle interactions
    balanceMassMomentumHost(pm->pn, pm->List, pm->Material, pm->Mass, pm->Smooth, 
                            pm->PosX, pm->PosY, pm->PosZ, 
                            pm->VelX, pm->VelY, pm->VelZ, 
                            pm->Density, pm->Pressure, pm->Sound, 
                            pm->DensityDot, pm->VelDotX, pm->VelDotY, pm->VelDotZ);

    balanceEnergyHost(pm->pn, pm->Pressure, pm->Density, 
                      pm->DensityDot, pm->EnergyDot);

        
    // Update particles
    updateParticlesHost(pm->pn, alpha, pm->Material, 
                        pm->VelDotX, pm->VelDotY, pm->VelDotZ, pm->DensityDot, pm->EnergyDot,
                        pm->PosX0, pm->PosY0, pm->PosZ0, 
                        pm->VelX0, pm->VelY0, pm->VelZ0, pm->Density0, pm->Energy0, 
                        pm->PosX, pm->PosY, pm->PosZ, pm->VelX, pm->VelY, pm->VelZ, 
                        pm->Density, pm->Energy, pm->Pressure, pm->Sound);
    /*
    printf("pos - %f %f\n", pm->PosX[101], pm->PosY[101]);
    printf("vel - %f %f\n", pm->VelX[101], pm->VelY[101]);
    printf("st  - %f %f\n", pm->Density[101], pm->Pressure[101]);
    printf("dot - %f %f %f\n", pm->DensityDot[101], pm->VelDotX[101], pm->VelDotY[101]);
    */
    return 0;
}


int RKstepDevice(struct device_model *dm, struct pointer_model *pm, float alpha) {
    int blocks, threads;

    blocks = (pm->pn + THREADS - 1) / THREADS;
    threads = THREADS;

	try
	{
		thrust::fill(dm->VelDotX.begin(), dm->VelDotX.end(), 0.0f);
		thrust::fill(dm->VelDotY.begin(), dm->VelDotY.end(), 0.0f);
		thrust::fill(dm->VelDotZ.begin(), dm->VelDotZ.end(), 0.0f);
		thrust::fill(dm->DensityDot.begin(), dm->DensityDot.end(), 0.0f);
		thrust::fill(dm->EnergyDot.begin(), dm->EnergyDot.end(), 0.0f);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in fill: %s\n", e.what());
		printf("pn: %d %d\n", dm->pn, pm->pn);
		printf("vdx: %d %d\n", dm->VelDotX.size(), dm->VelDotX.end() - dm->VelDotX.begin());
		printf("vdy: %d %d\n", dm->VelDotY.size(), dm->VelDotY.end() - dm->VelDotY.begin());
		printf("vdz: %d %d\n", dm->VelDotZ.size(), dm->VelDotZ.end() - dm->VelDotZ.begin());
		printf("dd: %d %d\n", dm->DensityDot.size(), dm->DensityDot.end() - dm->DensityDot.begin());
		printf("ed: %d %d\n", dm->EnergyDot.size(), dm->EnergyDot.end() - dm->EnergyDot.begin());
	}
	
	try
	{
		// External loads
		updateLoadsDevice <<< blocks, threads >>>
		(pm->pn, pm->Material, pm->PosX, pm->PosY, pm->PosZ, 
		pm->VelDotX, pm->VelDotY, pm->VelDotZ, pm->EnergyDot);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in loads: %s\n", e.what());
	}
	
	try
	{
		// Calculate particle interactions
		balanceMassMomentumDevice <<< blocks, threads >>>
		(pm->pn, pm->List, pm->Material, pm->Mass, pm->Smooth, pm->PosX, pm->PosY, pm->PosZ,
		pm->VelX, pm->VelY, pm->VelZ, pm->Density, pm->Pressure, pm->Sound,
		pm->DensityDot, pm->VelDotX, pm->VelDotY, pm->VelDotZ);
	
		balanceEnergyDevice <<< blocks, threads >>>
		(pm->pn, pm->Pressure, pm->Density, pm->DensityDot, pm->EnergyDot);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in balance: %s\n", e.what());
	}
	
	try
	{
		// Update particles
		updateParticlesDevice  <<< blocks, threads >>>
		(pm->pn, alpha, pm->Material, pm->VelDotX, pm->VelDotY, pm->VelDotZ, pm->DensityDot, pm->EnergyDot,
		pm->PosX0, pm->PosY0, pm->PosZ0, pm->VelX0, pm->VelY0, pm->VelZ0, pm->Density0, pm->Energy0,
		pm->PosX, pm->PosY, pm->PosZ, pm->VelX, pm->VelY, pm->VelZ, pm->Density, pm->Energy, pm->Pressure, pm->Sound);
	}
	catch(thrust::system_error &e)
	{
		printf("Error in update: %s\n", e.what());
	}

    return 0;
}


int RKintegrateHost(struct host_model *hm, struct pointer_model *pm) {

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

    // Save initial condition
    backupDataHost(hm);

    // TIME CYCLE
    for (ts = 0; ts <= hRun.tsn; ts++) {

        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
            resizeHost(hm);
            printData(hm);
            outputVTK(hm, ts / hRun.ssi);
        }
        
        // Calculate neighbour list
        neighbourListHost(hm, pm);
		
        // Save initial condition
        backupDataHost(hm);

        // Step 1
        RKstepHost(hm, pm, 1.0f);

        // Step 2
        RKstepHost(hm, pm, 1.0f / 4.0f);

        // Step 3
        RKstepHost(hm, pm, 2.0f / 3.0f);
        
        makeInletHost(hm, pm);
    }

    return 0;
}


int RKintegrateDevice(struct host_model *hm, struct device_model *dm, struct pointer_model *pm) {

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
	size_t available, total;


    // TIME CYCLE
    for (ts = 0; ts <= hRun.tsn; ts++) {

        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
			hm->pn = dm->pn;
            resizeHost(hm);
            copyDeviceToHost(dm, hm);
            printf("Particles: %d \n", hm->pn);
			cudaMemGetInfo(&available, &total);
			printf("Available memory %d MB\n", available/1024/1024);
            printData(hm);
            outputVTK(hm, ts / hRun.ssi);
        }
        
        try
        {
			// Calculate neighbour list
			neighbourListDevice(dm, pm);
		}
		catch(thrust::system_error &e)
		{
			printf("Error in neighbourList: %s", e.what());
			exit(-1);
		}
		
        try
        {
			// Save initial condition
			backupDataDevice(dm);
		}
		catch(thrust::system_error &e)
		{
			printf("Error in backupData: %s", e.what());
			exit(-1);
		}
		
        try
        {
			// Step 1
			RKstepDevice(dm, pm, 1.0f);
			
			// Step 2
			RKstepDevice(dm, pm, 1.0f / 4.0f);
			
			// Step 3
			RKstepDevice(dm, pm, 2.0f / 3.0f);
		}
		catch(thrust::system_error &e)
		{
			printf("Error in RKstep: %s", e.what());
			exit(-1);
		}
        
        try
        {
			makeInletDevice(dm, pm);
		}
		catch(thrust::system_error &e)
		{
			printf("Error in makeInlet: %s", e.what());
			exit(-1);
		}
        
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
     * \date May 2, 2012
     * \author Luca Massidda
     */

    struct host_model hModel;
    struct device_model dModel;
    struct pointer_model pModel;
	thrust :: host_vector <int > v;
	
	size_t available, total;
	
	cudaMemGetInfo(&available, &total);
	printf("Total memory %d MB\n", total/1024/1024);
	printf("Available memory %d MB\n", available/1024/1024);
	 
    for (int i = 0; i < 10; i++) {
        hLoad[i].gx = 0.0f;
        hLoad[i].gy = 0.0f;
        hLoad[i].gz = 0.0f;
        hLoad[i].w = 0.0f;

        hOut[i].nX = 0.0f;
        hOut[i].nY = 0.0f;
        hOut[i].nZ = 0.0f;
    }

    //initSingle(&hModel);
    //initPump(&hModel);
    //initBlock(&hModel);
    //initDamBreak(&hModel);
    //initChannel(&hModel);
    initBath(&hModel);
    
    resizeHost(&hModel);
	/*
    copyHostToDevice(&hModel, &dModel);
    copyDeviceToPointer(&dModel, &pModel);
    
    cudaMemcpyToSymbol("dMatType", hMatType, 10 * sizeof(int));
    cudaMemcpyToSymbol("dMatProp", hMatProp, 100 * sizeof(float));
    cudaMemcpyToSymbol("dRun", &hRun, sizeof(struct simulation));
    cudaMemcpyToSymbol("dLoad", &hLoad, 10 * sizeof(struct load));
    cudaMemcpyToSymbol("dFix", &hFix, 10 * sizeof(struct fix));
    cudaMemcpyToSymbol("dOut", &hOut, 10 * sizeof(struct outlet));
	
	cudaMemGetInfo(&available, &total);
	printf("Available memory %d MB\n", available/1024/1024);
	
	printf("%d\n", (&dModel)->Index.end()-(&dModel)->Index.begin());
	v.resize(hModel.pn);
    thrust::copy((&dModel)->Index.begin(), (&dModel)->Index.end(), v.begin());
	printf("%d\n", v[1]);
	thrust::sequence((&dModel)->Index.begin(), (&dModel)->Index.end(), 1);
	
    RKintegrateDevice(&hModel, &dModel, &pModel);
    */
    copyHostToPointer(&hModel, &pModel);
    RKintegrateHost(&hModel, &pModel);
	
    return 0;
}
