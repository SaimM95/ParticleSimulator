#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <time.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222


void drawParticle(unsigned char* img, int w,int h,int x, int y, int r, int cr,int cg, int cb){
  int minX = x - r;
  int maxX = x + r;
  printf("going to draw particle\n");
  for(int x = minX; x < maxX; x++){
    if(x < 0) continue;
    if(x >= w) break;
    //int minY = y - sqrt(r*r - x*x);
    //int maxY = y + sqrt(z);
    int minY = y - r;
    int maxY = y + r;
    for(int y = minY; y < maxY; y++){
      if(y < 0) continue;
      if(y >= h) break;
      img[(x + y*w)*3] = cr;
      img[(x + y*w)*3+1] = cg;
      img[(x + y*w)*3+2] = cb;
    }
  }
}
unsigned char* createImage(double* pos, int w,int h, int nl, int nm, int nh){
  int size = w*h*3;
  unsigned char* img = (unsigned char *)malloc(sizeof(char)*size);
  for(int i =0; i < size; i++){
    img[i] = 0;
  }

  int radius = 2;
  for(int i =0; i < nl; i++){
    int x = (int)pos[i*2];
    int y = (int)pos[i*2+1];
    drawParticle(img,w,h,x,y,radius,255,255,255);
    printf("setting value at (%i,%i) to 255\n", x,y);
  }

  return img;
}




double* genTestArr(int size, int start) {
  double* positions = (double*) malloc(sizeof(double) * size);
  int val = start;
  for (int i = 0; i < size; ++i) {
    if (i > 0 && i%2 == 0) {
      val++;
    }
    positions[i] = (double) val;
  }
  return positions;
}

double* genRandomArr(int size, double min, double max) {
  double* positions = (double*) malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i) {
    positions[i] = (drand48() * (max - min)) + min;
  }
  return positions;
}

void printVecArr(double *arr, int size) {
  for (int i = 0; i < size; i += 2) {
    printf("(%.4f,%.4f)\n", arr[i], arr[i+1]);
  }
}

int main(int argc, char* argv[]){
  if( argc != 10){
    printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
  }

  MPI_Init(&argc,&argv);

  int p, my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int nLight = atoi(argv[1]);
  int nMedium = atoi(argv[2]);
  int nHeavy = atoi(argv[3]);
  int nSteps = atoi(argv[4]);
  int subSteps = atoi(argv[5]);
  double timeSubStep = atof(argv[6]);
  int width = atoi(argv[7]);
  int height = atoi(argv[8]);

  double * pos = (double*)malloc(sizeof(double) * 2 * 10);
  for(int i =0; i < 10; i++){
    pos[i*2] = 5;
    pos[i*2+1] = 5;
  }
  printf("width: %i, height:%i\n", width, height);
  unsigned char* img;

  //root node stuff goes here
  if(my_rank == 0){
    // set seed for random number generation
    srand48(time(NULL));

    int numParticles = nLight + nMedium + nHeavy;

    int size = numParticles * 2;
    double *positions = genTestArr(size, 0);
    double *velocities = genTestArr(size, 1);

    printf("Positions:\n");
    printVecArr(positions, size);

    printf("Velocities:\n");
    printVecArr(velocities, size);

    unsigned char* img = createImage(pos, width, height, nLight,nMedium,nHeavy);

    //almost done, just save the img
    //saveBMP(argv[9], img, width, height);

    //almost done, just save the img
    // saveBMP(argv[9], img, width, height);
  }
  //all other nodes do this
  else {

  }

  //free(img);

  MPI_Finalize();
  return 0;
}
