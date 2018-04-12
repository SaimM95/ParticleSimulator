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


/* calls the f function on the amount of time this process needs
 *
 */
void calculate(double *start_pos, double * local_pos, double *forces){
  ///
}

/* calcualtes the force particle one exerts on particle 2
 */
double f(double m1, double m2, double *p1, double *p2){
  double G = -0.00000000006673;
  double totalMass = m1 * m2;
  double dist = sqrt((p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[0]-p2[0])*(p1[0]-p2[0]));
  return G * totalMass / (dist*dist*dist);
}

// updates the positions of the particles
void updatePos(double *forces, double *pos, double *vel, int w, int h){

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

void scatter(double* arr, double* l_arr, int l_size) {
  MPI_Scatter(arr, l_size, MPI_DOUBLE, // send one row, which contains l_size integers
              l_arr, l_size, MPI_DOUBLE, // receive one row, which contains l_size integers
              0, MPI_COMM_WORLD); // sent from root node 0
}

void sendForces(double *forces, int myRank, int p){

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

  printf("width: %i, height:%i\n", width, height);
  unsigned char* img;

  int numParticles = 4;//nLight + nMedium + nHeavy;
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
    printVecArr(positions, size);

    printf("Velocities:\n");
    printVecArr(velocities, size);
    printf("\n");

  }


  int l_size = size / p;
  double *l_pos = (double*) malloc(sizeof(double) * l_size);
  double *l_vel = (double*) malloc(sizeof(double) * l_size);

  scatter(positions, l_pos, l_size);
  scatter(velocities, l_vel, l_size);

  // init stuff
  int blockSize = ceil(numParticles/p);
  double * forces = (double *)malloc(sizeof(double) * numParticles * blockSize);

  for(int step = 0; step < nSteps; step++){
    for(int substep = 0; substep < subSteps; substep++){
      // recieve

      double * localPos;
      calculate(positions, localPos, forces);

      // send
    }

    sendForces(forces, my_rank, p);

    if(my_rank == 0){
      //unsigned char* img = createImage(pos, width, height, nLight,nMedium,nHeavy);
      //saveBMP(argv[9], img, width, height);
      //free(img);
    }
  }

  MPI_Finalize();
  return 0;
}

