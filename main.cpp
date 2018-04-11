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

  int numParticlesLight = 2;//atoi(argv[1]);
  int numParticleMedium = 2;//atoi(argv[1]);
  int numParticleHeavy = 2;//atoi(argv[1]);
  // int numSteps = argv[4];
  // int subSteps = argv[5];
  // double timeSubStep = argv[6];
  // int imageWidth = argv[7];
  // int imageHeight = argv[8];
  // int imageFilenamePrefix = argv[9];

  int width, height;

  unsigned char* image;

  //root node stuff goes here
  if(my_rank == 0){
    // set seed for random number generation
    srand48(time(NULL));

    int numParticles = 4;//numParticlesLight + numParticleMedium + numParticleHeavy;

    printf("Positions:\n");
    int size = numParticles * 2;
    double *positions = genTestArr(size, 0);// genRandomArr(size, 0, 10);
    printVecArr(positions, size);    

    printf("Velocities:\n");
    double *velocities = genTestArr(size, 1);// genRandomArr(size, 0, 10);
    printVecArr(velocities, size);

    //almost done, just save the image
    // saveBMP(argv[9], image, width, height);
  }
  //all other nodes do this
  else {

  }

  // free(image);

  MPI_Finalize();
  return 0;
}
