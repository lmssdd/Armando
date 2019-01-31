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
#define MAXP 300000
#define MAXS 3000000
#define PI 3.14159f

#define THREADS 256

struct pair {
	int key;
	int value;
};

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


void initBox(struct model *hm) {

    int i, j, k, m, b, q, ip;
    double rho, c0, pmin;
    double dr;
	
    q = 2;

    m = 1;
    b = 2;
    rho = 1000.f;
    c0 = 20.0f;
    pmin = -1.e4;

    hMatType[m] = 3;
    hMatProp[m][0] = rho;
    hMatProp[m][1] = c0;
    hMatProp[m][2] = pmin;

    hMatType[b] = 3;
    hMatProp[b][0] = rho;
    hMatProp[b][1] = c0;
    hMatProp[b][2] = pmin;

    dr = 0.02f / q; // x4
	
	ip = 0;
    for (k = 0; k < 2 * q ; k++) {
		for (j = 0; j < 2 * q ; j++) {
			for (i = 0; i < 2 * q; i++) {
				hm->PosX[ip] = i * dr;
				hm->PosY[ip] = j * dr;
				hm->PosZ[ip] = k * dr;
				hm->Material[ip] = m;
				hm->VelX[ip] = 0.0f;
				ip++;
			}
		}
	}

    for (k = 0; k < 2 * q ; k++) {
		for (j = -2; j < -1; j++) {
			for (i = -0; i < 2 * q; i++) {
				hm->PosX[ip] = i * dr;
				hm->PosY[ip] = j * dr;
				hm->PosZ[ip] = k * dr;
				hm->Material[ip] = b;
				hm->VelX[ip] = 0.0f;
				ip++;
			}
		}
	}
	
    hm->pn = ip;
    
    thrust::fill(hm->Mass, hm->Mass + hm->pn, rho * dr * dr * dr);
    thrust::fill(hm->Smooth, hm->Smooth + hm->pn, 1.2f * dr);
    thrust::fill(hm->VelY, hm->VelY + hm->pn, 0.0f);
    thrust::fill(hm->VelZ, hm->VelZ + hm->pn, 0.0f);
    thrust::fill(hm->Density, hm->Density + hm->pn, rho);
    thrust::fill(hm->Energy, hm->Energy + hm->pn, 0.0f);
    thrust::fill(hm->Pressure, hm->Pressure + hm->pn, 0.0f);
    thrust::fill(hm->Sound, hm->Sound + hm->pn, c0);
	
    hRun.minX = -0.2f;
    hRun.maxX =  1.0f;
    hRun.minY = -0.2f;
    hRun.maxY =  1.0f;
    hRun.minZ = -0.2f;
    hRun.maxZ =  1.0f;

    hRun.dt = dr / hMatProp[m][1];
    hRun.tsn = 1000 * q; //1000;
    hRun.ssi = 100 * q;

    hGrid.oX = hRun.minX;
    hGrid.oY = hRun.minY;
    hGrid.oZ = hRun.minZ;
    hGrid.size = 2.0f * 1.2f * dr;
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
    
    printf("Box \n");
    printf("Particles: %i \n", hm->pn);
    printf("Grid: %i \n", hGrid.nX * hGrid.nY * hGrid.nZ);
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

    cudaMemcpy(dm->SetStart, hm->SetStart, (MAXS * sizeof(int)), cudaMemcpyHostToDevice);
    cudaMemcpy(dm->SetStop, hm->SetStop, (MAXS * sizeof(int)), cudaMemcpyHostToDevice);

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

    cudaMemcpy(hm->SetStart, dm->SetStart, (MAXS * sizeof(int)), cudaMemcpyDeviceToHost);
    cudaMemcpy(hm->SetStop, dm->SetStop, (MAXS * sizeof(int)), cudaMemcpyDeviceToHost);

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
	
    hm->SetStart = (int *) malloc(MAXS * sizeof(int));
    hm->SetStop = (int *) malloc(MAXS * sizeof(int));
    
    return 0;
}


int initDevice(struct model *dm) {
	
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
	
    cudaMalloc((void**) &(dm->SetStart), (MAXS * sizeof(int)));
    cudaMalloc((void**) &(dm->SetStop), (MAXS * sizeof(int)));
    
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
    // Stream file add1
    stream = fopen("new_debug.txt", "w");
    for (i = 0; i < hm->pn; i++)
        fprintf(stream, "%d %f %f %f %f %f %f\n", i, hm->VelX[i], hm->VelY[i], hm->VelZ[i], hm->Density[i], hm->Energy[i], hm->Pressure[i]);
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
			PosX[ip] = 0.5f*(minX + maxX);
			PosY[ip] = 0.5f*(minY + maxY);
			PosZ[ip] = 0.5f*(minZ + maxZ);
			//Hash[ip] = Grid.nX * Grid.nY * Grid.nZ;
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

/*
__host__ void makeInletHost(struct model *hm, struct pointer_model *pm) {
	
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


void makeInletDevice(struct model *dm) {
	
	int i, j, iu, iv;
	
    for (i = 0; i < 10; i++) if (dIn[i].Material != 0) {
		dIn[i].Distance += dIn[i].Velocity * hRun.dt;
		
		if ((1.2f * dIn[i].Distance) > dIn[i].Smooth) {
			dIn[i].Distance = 0.0f;
			printf("Inlet!\n");
			for (iv = 0; iv < dIn[i].nv; iv++) {
				for (iu = 0; iu < dIn[i].nu; iu++) {
					dm->pn++;
					
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

}
*/

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

int neighbourListDevice(struct model *dm) {

    int blocks, threads;

    blocks = (dm->pn + THREADS - 1) / THREADS;
    threads = THREADS;
    
	updateHashDevice <<< blocks, threads >>>
	(dm->pn, dGrid, dm->PosX, dm->PosY, dm->PosZ, dm->Hash);
	
	checkOutletDevice <<< blocks, threads >>>
	(dm->pn, dGrid, dm->PosX, dm->PosY, dm->PosZ, dm->Hash);
	
	checkBoundariesDevice <<< blocks, threads >>>
	(dm->pn, dGrid, dm->PosX, dm->PosY, dm->PosZ, dm->Hash);
	
	// wrap raw pointer with a device_ptr 
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
	
	// use device_ptr in thrust algorithms
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
	
	thrust::device_ptr<int> tSetStart(dm->SetStart);
	thrust::device_ptr<int> tSetStop(dm->SetStop);
	thrust::fill(tSetStart, tSetStart + dm->pn, 0);
	thrust::fill(tSetStop, tSetStop + dm->pn, 0);
	
	updateSetsDevice <<< blocks, threads >>>
	(dm->pn, dm->SetStart, dm->SetStop, dm->Hash);
	
	updateListDevice <<< blocks, threads >>>
	(dm->pn, dm->List, dm->SetStart, dm->SetStop, dGrid, dm->Smooth,
	dm->PosX, dm->PosY, dm->PosZ);
	
    return 0;
}

int iSort(int *array, int *perm, int *dummy, int n) {
    int i;

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


int neighbourListHost(struct model *hm) {
	struct pair map[MAXP];
	
    updateHashHost(hm->pn, hGrid, hm->PosX, hm->PosY, hm->PosZ, hm->Hash);
	/*
	checkOutletHost(hm->pn, hGrid, hm->PosX, hm->PosY, hm->PosZ, hm->Hash);
	
	checkBoundariesHost(hm->pn, hGrid, hm->PosX, hm->PosY, hm->PosZ, hm->Hash);
	*/
	
	//thrust::sequence(hm->Index, hm->Index + hm->pn);
    for (int ip = 0; ip < hm->pn; ip++) hm->Index[ip] = ip;
	//for (int ip = 0; ip < hm->pn; ip++) printf("%d\n", hm->Hash[ip]),
    for (int ip = 0; ip < hm->pn; ip++) {
		map[ip].key = hm->Hash[ip];
		map[ip].value = hm->Index[ip];
	}
	qsort(map, hm->pn, sizeof(struct pair), mapCompare);
    for (int ip = 0; ip < hm->pn; ip++) {
		hm->Hash[ip] = map[ip].key;
		hm->Index[ip] = map[ip].value;
	}
    //for (int ip = 0; ip < hm->pn; ip++) hm->Index[ip] = map[ip].value;
	
	//thrust::sort_by_key(hm->Hash, hm->Hash + hm->pn, hm->Index);
	/*
	thrust::copy(hm->Material, hm->Material + hm->pn, hm->IntDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->IntDummy, hm->Material);
	
	thrust::copy(hm->Mass, hm->Mass + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->Mass);
	
	thrust::copy(hm->Smooth, hm->Smooth + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->Smooth);
	
	thrust::copy(hm->PosX, hm->PosX + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->PosX);
	
	thrust::copy(hm->PosY, hm->PosY + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->PosY);
	
	thrust::copy(hm->PosZ, hm->PosZ + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->PosZ);
	
	thrust::copy(hm->VelX, hm->VelX + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->VelX);
	
	thrust::copy(hm->VelY, hm->VelY + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->VelY);
	
	thrust::copy(hm->VelZ, hm->VelZ + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->VelZ);
	
	thrust::copy(hm->Density, hm->Density + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->Density);
	
	thrust::copy(hm->Energy, hm->Energy + hm->pn, hm->FloatDummy);
	thrust::gather(hm->Index, hm->Index + hm->pn, hm->FloatDummy, hm->Energy);
	*/
    //iSort(hm->Hash, hm->Index, hm->pn);
    iSort(hm->Material, hm->Index, hm->IntDummy, hm->pn);
    fSort(hm->Mass, hm->Index, hm->pn);
    fSort(hm->Smooth, hm->Index, hm->pn);
    fSort(hm->PosX, hm->Index, hm->pn);
    fSort(hm->PosY, hm->Index, hm->pn);
    fSort(hm->PosZ, hm->Index, hm->pn);
    fSort(hm->VelX, hm->Index, hm->pn);
    fSort(hm->VelY, hm->Index, hm->pn);
    fSort(hm->VelZ, hm->Index, hm->pn);
    fSort(hm->Density, hm->Index, hm->pn);
    fSort(hm->Energy, hm->Index, hm->pn);
	
	//thrust::fill(hm->SetStart, hm->SetStart + hm->pn, 0);
    for (int i = 0; i < hGrid.nX * hGrid.nY * hGrid.nZ; i++) hm->SetStart[i] = 0;
	//thrust::fill(hm->SetStop, hm->SetStop + hm->pn, 0);
    for (int i = 0; i < hGrid.nX * hGrid.nY * hGrid.nZ; i++) hm->SetStop[i] = 0;
	
	updateSetsHost(hm->pn, hm->SetStart, hm->SetStop, hm->Hash);
	
	updateListHost(hm->pn, hm->List, hm->SetStart, hm->SetStop, hGrid, hm->Smooth,
				   hm->PosX, hm->PosY, hm->PosZ);
	
    return 0;
}

int backupDataHost(struct model *hm) {

    thrust::copy(hm->PosX, hm->PosX + hm->pn, hm->PosX0);
    thrust::copy(hm->PosY, hm->PosY + hm->pn, hm->PosY0);
    thrust::copy(hm->PosZ, hm->PosZ + hm->pn, hm->PosZ0);
    thrust::copy(hm->VelX, hm->VelX + hm->pn, hm->VelX0);
    thrust::copy(hm->VelY, hm->VelY + hm->pn, hm->VelY0);
    thrust::copy(hm->VelZ, hm->VelZ + hm->pn, hm->VelZ0);
    thrust::copy(hm->Density, hm->Density + hm->pn, hm->Density0);
    thrust::copy(hm->Energy, hm->Energy + hm->pn, hm->Energy0);

    return 0;
}

int backupDataHostOld(struct model *hm) {

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

	// wrap raw pointer with a device_ptr 
	thrust::device_ptr<float> tPosX(dm->PosX);
	thrust::device_ptr<float> tPosY(dm->PosY);
	thrust::device_ptr<float> tPosZ(dm->PosZ);
	thrust::device_ptr<float> tVelX(dm->VelX);
	thrust::device_ptr<float> tVelY(dm->VelY);
	thrust::device_ptr<float> tVelZ(dm->VelZ);
	thrust::device_ptr<float> tDensity(dm->Density);
	thrust::device_ptr<float> tEnergy(dm->Energy);
	thrust::device_ptr<float> tPosX0(dm->PosX0);
	thrust::device_ptr<float> tPosY0(dm->PosY0);
	thrust::device_ptr<float> tPosZ0(dm->PosZ0);
	thrust::device_ptr<float> tVelX0(dm->VelX0);
	thrust::device_ptr<float> tVelY0(dm->VelY0);
	thrust::device_ptr<float> tVelZ0(dm->VelZ0);
	thrust::device_ptr<float> tDensity0(dm->Density0);
	thrust::device_ptr<float> tEnergy0(dm->Energy0);
	
	// use device_ptr in thrust algorithms
    thrust::copy(tPosX, tPosX + dm->pn, tPosX0);
    thrust::copy(tPosY, tPosY + dm->pn, tPosY0);
    thrust::copy(tPosZ, tPosZ + dm->pn, tPosZ0);
    thrust::copy(tVelX, tVelX + dm->pn, tVelX0);
    thrust::copy(tVelY, tVelY + dm->pn, tVelY0);
    thrust::copy(tVelZ, tVelZ + dm->pn, tVelZ0);
    thrust::copy(tDensity, tDensity + dm->pn, tDensity0);
    thrust::copy(tEnergy, tEnergy + dm->pn, tEnergy0);
    
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

int RKstepHost(struct model *hm, float alpha) {
	/*
    thrust::fill(hm->VelDotX, hm->VelDotX + hm->pn, 0.0f);
    thrust::fill(hm->VelDotY, hm->VelDotY + hm->pn, 0.0f);
    thrust::fill(hm->VelDotZ, hm->VelDotZ + hm->pn, 0.0f);
    thrust::fill(hm->DensityDot, hm->DensityDot + hm->pn, 0.0f);
    thrust::fill(hm->EnergyDot, hm->EnergyDot + hm->pn, 0.0f);
	*/
	for (int ip = 0; ip < hm->pn; ip++) {
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
	
    // Calculate particle interactions
    balanceMassMomentumHost(hm->pn, hm->List, hm->Material, hm->Mass, hm->Smooth, 
                            hm->PosX, hm->PosY, hm->PosZ, 
                            hm->VelX, hm->VelY, hm->VelZ, 
                            hm->Density, hm->Pressure, hm->Sound, 
                            hm->DensityDot, hm->VelDotX, hm->VelDotY, hm->VelDotZ);

    balanceEnergyHost(hm->pn, hm->Pressure, hm->Density, 
                      hm->DensityDot, hm->EnergyDot);

	
    // Update particles
    updateParticlesHost(hm->pn, alpha, hm->Material, 
                        hm->VelDotX, hm->VelDotY, hm->VelDotZ, hm->DensityDot, hm->EnergyDot,
                        hm->PosX0, hm->PosY0, hm->PosZ0, 
                        hm->VelX0, hm->VelY0, hm->VelZ0, hm->Density0, hm->Energy0, 
                        hm->PosX, hm->PosY, hm->PosZ, hm->VelX, hm->VelY, hm->VelZ, 
                        hm->Density, hm->Energy, hm->Pressure, hm->Sound);
	
    for (int ip = 0; ip < hm->pn; ip++) {
		printf("%d %d %f %f %f %f %f \n", hm->Index[ip], hm->Hash[ip], hm->VelX[ip], hm->VelY[ip], hm->VelZ[ip], hm->Density[ip], hm->Pressure[ip]);
	}
    
    return 0;
}


int RKstepDevice(struct model *dm, float alpha) {
    int blocks, threads;

    blocks = (dm->pn + THREADS - 1) / THREADS;
    threads = THREADS;

	// wrap raw pointer with a device_ptr 
	thrust::device_ptr<float> tVelDotX(dm->VelDotX);
	thrust::device_ptr<float> tVelDotY(dm->VelDotY);
	thrust::device_ptr<float> tVelDotZ(dm->VelDotZ);
	thrust::device_ptr<float> tDensityDot(dm->DensityDot);
	thrust::device_ptr<float> tEnergyDot(dm->EnergyDot);
	
	// use device_ptr in thrust algorithms
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
        RKstepHost(hm, 1.0f);
		/*
        // Step 2
        RKstepHost(hm, 1.0f / 4.0f);

        // Step 3
        RKstepHost(hm, 2.0f / 3.0f);
        */
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
	size_t available, total;
	
	
    // TIME CYCLE
//    for (ts = 0; ts <= hRun.tsn; ts++) {
    for (ts = 0; ts < 1; ts++) {
		
		// Calculate neighbour list
		neighbourListDevice(dm);
		
		// Save initial condition
		backupDataDevice(dm);
		
		// Step 1
		RKstepDevice(dm, 1.0f);
		/*
		// Step 2
		RKstepDevice(dm, 1.0f / 4.0f);
		
		// Step 3
		RKstepDevice(dm, 2.0f / 3.0f);
		*/
        // Output data
        if ((ts % hRun.ssi) == 0) {
            printf("Saving time: %g \n", ts * hRun.dt);
            copyDeviceToHost(dm, hm);
            printf("Particles: %d \n", hm->pn);
			cudaMemGetInfo(&available, &total);
			printf("Available memory %d MB\n", available/1024/1024);
            printData(hm);
            outputVTK(hm, ts / hRun.ssi);
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

    struct model hModel, dModel;
	
	size_t available, total;
	
	cudaMemGetInfo(&available, &total);
	printf("Occupied memory %d of %dMB\n", available/1024/1024, total/1024/1024);
	 
    for (int i = 0; i < 10; i++) {
        hLoad[i].gx = 0.0f;
        hLoad[i].gy = 0.0f;
        hLoad[i].gz = 0.0f;
        hLoad[i].w = 0.0f;

        hOut[i].nX = 0.0f;
        hOut[i].nY = 0.0f;
        hOut[i].nZ = 0.0f;
    }
    
    initHost(&hModel);
    
    //initSingle(&hModel);
    //initPump(&hModel);
    //initBlock(&hModel);
    //initDamBreak(&hModel);
    //initChannel(&hModel);
    initBox(&hModel);
    /*
    initDevice(&dModel);
	cudaMemGetInfo(&available, &total);
	printf("Available memory %d MB\n", available/1024/1024);
	
    //copyHostToDevice(&hModel, &dModel);
    
	cudaMemGetInfo(&available, &total);
	printf("Available memory %d MB\n", available/1024/1024);
	
	thrust::device_ptr<int> tIndex((&dModel)->Index);
	thrust::sequence(tIndex, tIndex + (&dModel)->pn, 1);
	
	
	neighbourListDevice(&dModel);
	copyDeviceToHost(&dModel, &hModel);
	
	printData(&hModel);
	*/
	//RKintegrateDevice(&hModel, &dModel);
    
    RKintegrateHost(&hModel);
	
    return 0;
}
