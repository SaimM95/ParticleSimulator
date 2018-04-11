#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){

  if( argc != 10){
    printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
  }

  MPI_Init(&argc,&argv);

  int p, my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int numParticlesLight = atoi(argv[1]);
  // int numParticleMedium = argv[2];
  // int numParticleHeavy = argv[3];
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



    //almost done, just save the image
    saveBMP(argv[9], image, width, height);
  }
  //all other nodes do this
  else{

  }

  free(image);

  MPI_Finalize();
  return 0;
}
