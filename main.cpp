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

void printVecArr(double *arr, int size, int P) {
  for (int i = 0; i < size; i += 2) {
    printf("%d: (%.4f,%.4f)\n", P, arr[i], arr[i+1]);
  }
}

void printArr(int *arr, int size) {
  for (int i = 0; i < size; ++i) {
    printf("%d,", arr[i]);
  }
  printf("\n");
}

void calcRowDistributions(int rows, int proc, int *rowDistributions) {
    int rowsLeft = rows;
    int p;

    for (p = proc; p > 0; p--) {
        int avgRowBlock = (int) ceil((double)rowsLeft / p); // roundup
        int procIndex = proc - p;
        rowDistributions[procIndex] = avgRowBlock;
        rowsLeft -= avgRowBlock;

        if (rowsLeft <= 0) {
            break;
        }
    }
}

int* scatter(double *matrix, double *local_matrix, int rows, int cols, int p, int rank) {
  int *rowDistributions, *counts, *displacements;
  int local_rows, local_count, i;

  // how many rows per processor
  rowDistributions = (int*) malloc(sizeof(int) * p);
  calcRowDistributions(rows, p, rowDistributions);

  if (rank == p-1) {
    printf("Row Distributions:\n");
    printArr(rowDistributions, p);
  }

  // 1D mapping of rowDistributions
  // e.g. For a 5x4 matrix, {1,1,1,2} => {4,4,4,8}
  counts = (int*) malloc(sizeof(int) * p);
  for (i = 0; i < p; ++i) {
      counts[i] = rowDistributions[i] * cols;
  }

  // starting indices of each row in the 1D mapping
  // e.g. For counts = {4,4,4,8} => {0,4,8,12}
  displacements = (int*) malloc(sizeof(int) * p);
  displacements[0] = 0;
  for (i = 1; i < p; ++i) {
      displacements[i] = displacements[i-1] + counts[i-1];
  }

  local_rows = rowDistributions[rank];
  local_count = counts[rank];
  
  MPI_Scatterv(matrix, counts, displacements, MPI_DOUBLE, // send local_n rows, which contains m vals
              local_matrix, local_count, MPI_DOUBLE, // receive local_n rows, which contains m vals
              0, MPI_COMM_WORLD); // sent from root node 0
  
  return rowDistributions;
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
  char* filePrefix = argv[9];

  double * pos = (double*)malloc(sizeof(double) * 2 * 10);
  for(int i =0; i < 10; i++){
    pos[i*2] = 5;
    pos[i*2+1] = 5;
  }
  printf("width: %i, height:%i\n", width, height);
  unsigned char* img;

  int numParticles = 7;//numParticlesLight + numParticleMedium + numParticleHeavy;
  int size = numParticles * 2;
  double *positions, *velocities;

  //root node stuff goes here
  if(my_rank == 0){
    // set seed for random number generation
    srand48(time(NULL));

    int size = numParticles * 2;
    positions = genTestArr(size, 0);
    velocities = genTestArr(size, 1);

    printf("Positions:\n");
    printVecArr(positions, size, my_rank);

    printf("Velocities:\n");
    printVecArr(velocities, size, my_rank);
    printf("\n");

    unsigned char* img = createImage(pos, width, height, nLight,nMedium,nHeavy);

    //almost done, just save the img
    saveBMP(filePrefix, img, width, height);
    printf("\n");
  }

  int l_size = (int) ceil((double)numParticles / p) * 2;
  double *l_pos = (double*) malloc(sizeof(double) * l_size);
  double *l_vel = (double*) malloc(sizeof(double) * l_size);

  // FOR DEBUGGING ONLY
  for (int i = 0; i < l_size; ++i) {
    l_pos[i] = -1;
    l_vel[i] = -1;
  }

  int *rowDistributions;
  rowDistributions = scatter(positions, l_pos, numParticles, 2, p, my_rank);
  // following call should have same result as call to scatter above
  scatter(velocities, l_vel, numParticles, 2, p, my_rank);

  int l_size_actual = rowDistributions[my_rank] * 2;

  printf("\n");
  printf("Process %d has size:%d\n", my_rank, l_size_actual);
  printf("\n");

  printf("P - Process %d with:\n", my_rank);
  printVecArr(l_pos, l_size_actual, my_rank);
  printf("\n");

  printf("V - Process %d with:\n", my_rank);
  printVecArr(l_vel, l_size_actual, my_rank);

  // free(img);

  MPI_Finalize();
  return 0;
}
