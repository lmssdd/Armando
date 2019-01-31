//      untitled.c
//      
//      Copyright 2011 Luca Massidda <lucam@pinco>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

typedef struct tpoint {
	float x;
	float y;
	float z;
} point;

typedef struct ttriangle {
	point v[3];
} triangle;

float DotProduct(point a, point b) {
	float r;
	
	r = a.x*b.x + a.y*b.y + a.z*b.z;
	
	return r;
}

point CrossProduct(point a, point b) {
	point c;
	
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	
	return c;
}

point Normalize(point a) {
	point c;
	float d;
	
	d = sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
	c.x = a.x / d; 
	c.y = a.y / d; 
	c.z = a.z / d;
	
	return c;
}

int SameSide(point a, triangle t) {
	
	return 1;
}

int LineFacet(point p0, point p1, triangle t, point *p) {
	point u, v, w, n;
	float td, tn;
	float dotuu, dotvv, dotuv, dotwu, dotwv;
	float s1, s2, invDenom;
	
	u.x = t.v[1].x - t.v[0].x;
	u.y = t.v[1].y - t.v[0].y;
	u.z = t.v[1].z - t.v[0].z;
	
	v.x = t.v[2].x - t.v[0].x;
	v.y = t.v[2].y - t.v[0].y;
	v.z = t.v[2].z - t.v[0].z;
	
	n = CrossProduct(u, v);
	n = Normalize(n);
	
	w.x = t.v[0].x - p0.x;
	w.y = t.v[0].y - p0.y;
	w.z = t.v[0].z - p0.z;
	tn = DotProduct(w, n);
	
	w.x = p1.x - p0.x;
	w.y = p1.y - p0.y;
	w.z = p1.z - p0.z;
	td = DotProduct(w, n);
	if (fabs(td) < FLT_EPSILON) td = FLT_EPSILON;
	
	p->x = p0.x + (tn * w.x / td);
	p->y = p0.y + (tn * w.y / td);
	p->z = p0.z + (tn * w.z / td);
	
	
	w.x = p->x - t.v[0].x;
	w.y = p->y - t.v[0].y;
	w.z = p->z - t.v[0].z;
	
	// Compute dot products
	dotuu = DotProduct(u, u);
	dotvv = DotProduct(v, v);
	dotuv = DotProduct(u, v);
	dotwu = DotProduct(w, u);
	dotwv = DotProduct(w, v);
	
	// Compute barycentric coordinates
	invDenom = 1.0 / (dotuv * dotuv - dotuu * dotvv);
	s1 = (dotuv * dotwv - dotvv * dotwu) * invDenom;
	s2 = (dotuv * dotwu - dotuu * dotwv) * invDenom;
	
	// Check if point is in triangle
	return ((s1 >= 0.0) && (s2 >= 0.0) && (s1 + s2 <= 1.0));
}

int FacetRasterization(triangle t, point O, int n, float h, int *voxels) {
	point u, v, p;
	float a, b, da, db;
	int i, j, k;
	
	u.x = t.v[1].x - t.v[0].x;
	u.y = t.v[1].y - t.v[0].y;
	u.z = t.v[1].z - t.v[0].z;
	
	v.x = t.v[2].x - t.v[0].x;
	v.y = t.v[2].y - t.v[0].y;
	v.z = t.v[2].z - t.v[0].z;
	
	da = h / DotProduct(u, u);
	db = h / DotProduct(v, v);
	
	for (b = 0.0; b < 1.0; b += db) {
		for (a = 0.0; (a < 1.0) && (a + b <= 1.0); a += da) {
			p.x = t.v[0].x + a*u.x + b*v.x;
			p.y = t.v[0].y + a*u.y + b*v.y;
			p.z = t.v[0].z + a*u.z + b*v.z;
			
			i = (int) roundf((p.x - O.x) / h);
			j = (int) roundf((p.y - O.y) / h);
			k = (int) roundf((p.z - O.z) / h);
			
			if ((i >= 0) && (i < n) && 
				(j >= 0) && (j < n) && 
				(k >= 0) && (k < n)) voxels[i*n*n + j*n + k] = 1;
		}
	}
	
	return 0;
}


int main(int argc, char **argv) {
	FILE *fp;
	char infilename[80];
	char outfilename[80];
	char buffer[80];
	char *bufferp;
	int issolid = 0;
	int isfacet = 0;
	int nfacets = 0;
	int i, j, k, f;
	
	triangle *facet;
	point O, p;
	int n = 200;
	int *voxels;
	const float h = 1.0;
	
	O.x = -100.0; O.y = -100.0; O.z = -100.0;
	
	facet = (triangle*) malloc(1000000 * sizeof(triangle));
	voxels = (int*) malloc(n*n*n * sizeof(int));
	
	strcpy(infilename, argv[1]);
	strcpy(outfilename, argv[2]);
	
	if (!(fp = fopen(infilename, "r"))) {
		printf("Error:  Unable to open file %s\n", infilename);
		return 1;
	}
	printf("Opening %s\n", infilename);
	
	i = 0; j = 0;
	while (fgets(buffer, sizeof(buffer), fp) != NULL) {
		
		if (strstr(buffer, "solid") != NULL) issolid = 1;
		if (strstr(buffer, "endsolid") != NULL) issolid = 0;
		if (strstr(buffer, "facet") != NULL) isfacet = 1;
		if (strstr(buffer, "endfacet") != NULL) {
			i++;
			j = 0;
			isfacet = 0;
		}
		if (strstr(buffer, "vertex") != NULL) {
			if (!issolid) {
				printf("Error : not a solid!\n");
				return 1;
			}
			if (!isfacet) {
				printf("Error : not a facet!\n");
				return 1;
			}
			if (j > 2) {
				printf("Error : too many verts!\n");
				printf("Facet : %d\n",i);
				return 1;
			}
			bufferp = strstr(buffer, "vertex");
			
			sscanf(bufferp, "vertex %f %f %f", &p.x, &p.y, &p.z);
			facet[i].v[j] = p;
			j++;
		}
	}
	fclose(fp);
	nfacets = i;
	
	printf("Voxelising...\n");
	
	for (i = 0; i < n*n*n; i++) voxels[i] = 0;
	
	for (f = 0; f < nfacets; f++) {
		if (! (f % 1000)) printf("Facet %d / %d\n", f, nfacets);
		
		FacetRasterization(facet[f], O, n, h, voxels);
	}
	
	
	if (!(fp = fopen(outfilename, "w"))) {
		printf("Error:  Unable to open file %s\n", outfilename);
		return 1;
	}
	printf("Opening %s\n", outfilename);
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				if (voxels[i*n*n + j*n + k] > 0) {
					fprintf(fp, "%12f, %12f, %12f\n", 
						O.x + h * i, O.y + h * j, O.z + h * k);
				}
			}
		}
	}
	fclose(fp);
	
	
	return 0;
}

