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


void drawParticle(unsigned char* img, int w,int h,int x, int y, int r, vec3 col){
  int minX = x - r;
  int maxX = x + r;
  for(int cx = minX; cx < maxX; cx++){
    if(x < 0) continue;
    if(x >= w) break;
    int diff = r*r - (cx-x)*(cx-x);
    int val = (int)sqrt(abs(diff));
    int minY = y - val;
    int maxY = y + val;
    for(int y = minY; y < maxY; y++){
      if(y < 0) continue;
      if(y >= h) break;
      img[(cx + y*w)*3] = col.x*255;
      img[(cx + y*w)*3+1] = col.y*255;
      img[(cx + y*w)*3+2] = col.z*255;
    }
  }
}

void createImage(double* pos, int w,int h, int nl, int nm, int nh, unsigned char* img){
  int size = w*h*3;
  for(int i =0; i < size; i++){
    img[i] = 0;
  }

  int radius = 5;
  vec3 col = colourLight;
  int numPart = nl+nm+nh;
  for(int i =0; i < numPart; i++){
    if(i >= nm + nl) col = colourHeavy;
    else if(i >= nl) col = colourMedium;
    int x = (int)pos[i*3];
    int y = (int)pos[i*3+1];
    drawParticle(img,w,h,x,y,radius,col);
  }
}

int pixelsInMeter = 1000;
/* calcualtes the force particle one exerts on particle 2
 */
double calcForce(double *p1, double *p2){
  double G = -0.00000000006673;
  //double G = -0.1;
  double totalMass = p1[2] * p2[2];
  double dist = sqrt((p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[0]-p2[0])*(p1[0]-p2[0]));
  //if(dist < 1) return 0; // equal force in all direction
  dist = dist / pixelsInMeter;
  //return G * totalMass / (dist*dist*dist);
  //return G * totalMass / (dist*dist);
  //int pixelsToMeters = 10;
  //int distMeters = dist / (pixelsToMeters*pixelsToMeters);
  //double force = totalMass;
  //double force = G * totalMass / (dist*dist*dist);
  //double force = G * totalMass / (dist);
  double smallConstant = 0.0000033;
  double force = G * totalMass / (dist*dist*dist + smallConstant);
  // printf("the calculated force: %f dist:%f, totalMas: %f\n", force, dist, totalMass);
  return force;
  //return (dist*dist*dist);
}

int inline getNodeForProc(int proc, int p, int n){
  int rem = n%p;
  int avgBlockSize = n/p;
  if(proc < rem) return proc * (avgBlockSize+1);
  return rem* (avgBlockSize+1) + (proc-rem)*(avgBlockSize);
}

/* calls the f function on the amount of time this process needs
 *
 */
void calculate(double *startPos, double * localPos, double *forces, int blockSize, int rank, int n, int rankOther, int p){
  // printf("%d: starting calculate numpar:%d blockSize:%i\n", rank, n, blockSize);
  for(int i =0; i < blockSize; i++){
    int startIndex=0;
    if(startPos == localPos) startIndex = i;
    //if(startPos == localPos) startIndex = i+1;
    for(int j =startIndex; j < blockSize; j++){
      double * p1 = &startPos[i*3];
      double * p2 = &localPos[j*3];
      //if(p1 == p2) continue;
      // printf("%d: updating force value p1: (%f,%f,%f), p2:(%f,%f,%f)\n",rank, p1[0],p1[1],p1[2],  p2[0],p2[1],p2[2] );
      double f = calcForce(p1,p2);
      int newR = i;
      int newC = j + getNodeForProc(rankOther, p, n);
      int index = (newR * n + newC)*2;
      forces[index] = f * (p1[0] - p2[0])/pixelsInMeter;
      forces[index+1] = f * (p1[1] - p2[1])/pixelsInMeter;
      // printf("%d: calculate force f:%f updating value at (%i,%i) rankOther:%i n:%i\n", rank, f, newR, newC, rankOther, n);
    }
  }
}

// updates the positions of the particles
void updatePos(double *forces, double *pos, double *vel, int w, int h, int n, int blockSize){
  double *F = (double*) malloc(sizeof(double) * blockSize * 2);
  int timeStep = 1;

  for (int i = 0; i < blockSize; ++i) {
    F[i*2] = 0;
    F[i*2+1] = 0;
  }

  for (int i = 0; i < blockSize; ++i) {
    double forcesX = 0, forcesY = 0;

    for (int j = 0; j < n; ++j) {
      int ind = i*n*2 + j*2;
      forcesX += forces[ind];
      forcesY += forces[ind + 1];
      // printf("forcesX:%.4f forcesY:%.4f\n", forcesX, forcesY);
    }

    F[i*2] += forcesX;
    F[i*2+1] += forcesY;
  }

  // printf("Printing F\n");
  // for (int i = 0; i < blockSize; ++i) {
  //   printf("[%.4f,%.4f]\n", F[i*2], F[i*2+1]);
  // }

  for (int i = 0; i < blockSize; ++i) {
    pos[i*3] += timeStep * vel[i*2];
    pos[i*3+1] += timeStep * vel[i*2+1];

    double mass = pos[i*3+2];
    vel[i*2] += timeStep * F[i*2] * (1/mass);
    vel[i*2+1] += timeStep * F[i*2+1] * (1/mass);
  }

  // update velocities
  // for(int p = 0; p < blockSize; p++){
  //   double totalForceX = 0;
  //   double totalForceY = 0;
  //   for(int i =0; i < n; i++){
  //     totalForceX += forces[(i + n*p)*2 + 0];
  //     totalForceY += forces[(i + n*p)*2 + 1];
  //   }
  //   double mass = pos[p*3+2];
  //   vel[p*2] += totalForceX/mass;
  //   vel[p*2+1] += totalForceY/mass;
  // }

  // // update position
  // for(int i = 0; i < blockSize; i++){
  //   pos[i*3] += vel[i*2];
  //   pos[i*3+1] += vel[i*2+1];
  // }
}


double* genTestArr(int size, int start, int width) {
  double* positions = (double*) malloc(sizeof(double) * size);
  int val = start;
  for (int i = 0; i < size; ++i) {
    if (i > 0 && i%width == 0) {
      val++;
    }
    positions[i] = (double) val;
  }
  return positions;
}

double* genRandomPos(int width, int height, int nl, int nm, int nh) {
  int n = nl+nm+nh;
  double* positions = (double*) malloc(sizeof(double) * n *3);
  int max = massLightMax;
  int min = massLightMin;
  for (int i = 0; i < n; ++i) {
    if(i ==nl){
      max = massMediumMax;
      min = massMediumMin;
    }else if(i == nm){
      max = massHeavyMax;
      min = massHeavyMin;
    }
    positions[i*3] = (drand48()*width);
    positions[i*3+1] = (drand48()*height);
    positions[i*3+2] = (drand48()*(max-min)) + min;
  }
  return positions;
}

double* genRandomVel(int nl, int nm, int nh) {
  int n = nl+nm+nh;
  double* vel = (double*) malloc(sizeof(double) * n *2);
  int max = velocityLightMin;
  int min = velocityLightMax;
  for (int i = 0; i < n; ++i) {
    if(i ==nl){
      max = velocityMediumMax;
      min = velocityMediumMin;
    }else if(i == nm){
      max = velocityHeavyMax;
      min = velocityHeavyMin;
    }
    vel[i*2] = (drand48()*(max-min)) + min;
    vel[i*2+1] = (drand48()*(max-min)) + min;

    // velocity can be negative
    if(drand48()*2 < 1) vel[i*2] = vel[i*2]*-1;
    if(drand48()*2 < 1) vel[i*2+1] = vel[i*2+1]*-1;
  }
  return vel;
}

double* genTestPos2(){
  // test for testing forces
  int size = 2 * 3;
  double* positions = (double*) malloc(sizeof(double) * size);
  positions[0] = 500;
  positions[1] = 480;
  positions[2] = 1000.0f;

  positions[3] = 500;
  positions[4] = 520;
  positions[5] = 3000.0f;

  /*
  positions[6] = 700;
  positions[7] = 300;
  positions[8] = 3;

  positions[9] = 700;
  positions[10] = 500;
  positions[11] = 5;
  */
  return positions;
}
double * genVelZeros(int n){
  int size = n * 2;
  double* vel = (double*) malloc(sizeof(double) * size);
  for(int i =0; i < n*2; i++){
    vel[i] = 0;
  }
  return vel;
}

void printVecArr(double *arr, int size, int P, int width) {
  for (int i = 0; i < size; i += width) {
    printf("%d: (", P);
    for(int j=0; j < width; j++){
      printf("%.4f,", arr[i+j]);
    }
    printf(")\n");
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
  int local_count, i;

  // how many rows per processor
  rowDistributions = (int*) malloc(sizeof(int) * p);
  calcRowDistributions(rows, p, rowDistributions);

  // if (rank == p-1) {
  //   printf("Row Distributions:\n");
  //   printArr(rowDistributions, p);
  // }

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

  local_count = counts[rank];

  MPI_Scatterv(matrix, counts, displacements, MPI_DOUBLE, // send local_n rows, which contains m vals
              local_matrix, local_count, MPI_DOUBLE, // receive local_n rows, which contains m vals
              0, MPI_COMM_WORLD); // sent from root node 0

  return rowDistributions;
}
void gather(double *pos, double *pos_loc, int rows, int cols, int p, int rank){
  int *rowDistributions, *counts, *displacements;
  int i;

  rowDistributions = (int*) malloc(sizeof(int) * p);
  calcRowDistributions(rows, p, rowDistributions);

  counts = (int*) malloc(sizeof(int) * p);
  for (i = 0; i < p; ++i) {
      counts[i] = rowDistributions[i] * cols;
  }

  displacements = (int*) malloc(sizeof(int) * p);
  displacements[0] = 0;
  for (i = 1; i < p; ++i) {
      displacements[i] = displacements[i-1] + counts[i-1];
  }

  MPI_Gatherv(pos_loc, rowDistributions[rank]*3, MPI_DOUBLE,
    pos, counts, displacements,
    MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int inline getProcessForNode(int sendNode, int p, int n){
  int rem = n%p;
  int avgBlockSize = n/p;
  if(sendNode < rem + avgBlockSize * rem) return sendNode/(avgBlockSize+1);
  return rem + (sendNode - (avgBlockSize+1)*rem)/avgBlockSize;
}
void sendForces(double *forces, int rank, int p, int n, int blockSize){
  /*
  int num =0;
  int size = ceil(n/(float)p);
  double * tempArr = (double*)malloc(sizeof(double) * size * n*2);
  for(int i=0; i < size; i++){
    for(int j =0; j < n; j++){
      tempArr[(i*n+j)*2+0] = num;
      tempArr[(i*n+j)*2+1] = num;
      num++;
    }
  }
  printf("going to print the temp array\n");
  printVecArr(tempArr, blockSize*n*2, rank, n*2);
  */
  // recieve everything
  for(int r =0; r < blockSize; r++){
    int actualRow = r + getNodeForProc(rank, p, n);
    for(int c =0; c < actualRow; c++){
      int actualRow = r + getNodeForProc(rank, p, n);
      int transferIndex = getProcessForNode(c, p, n);
      int index = (r*n+c)*2;
      // printf("%d: - sendForces - going to recieve (%i,%i) transferIndex: %i\n", rank, r,c,transferIndex);
      if(transferIndex == rank) continue;
      //printf("%d: recieved value (%i,%i)\n", rank, r,c);
      MPI_Recv(&forces[index], 2, MPI_DOUBLE, transferIndex, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      forces[index] = forces[index]*-1;
      forces[index+1] = forces[index+1]*-1;
      // printf("%d: - sendForces - going to recieved the value for (%i,%i) from: %i value:(%f,%f)\n", rank, r,c,transferIndex, forces[index], forces[index+1]);
      //printf("%d: recieved value (%i,%i) of (%f,%f)\n", rank, r,c, forces[index], forces[index+1]);
      //printf("%d: updateForce - rec  - (%i,%i) actualRow: [%i] to %i got values:(%f,%f) index:%i\n", rank, r,c, actualRow, transferIndex, forces[index], forces[index+1], index);
    }
  }

  // send everything
  int actualBaseRow = 0 + getNodeForProc(rank, p, n);
  for(int c =0; c < n; c++){
    int rowEnd = c - actualBaseRow;
    if(rowEnd > blockSize) rowEnd = blockSize;
    //printf("%d: updateForce - send - rowEnd: %i, (_,%i)\n", rank, rowEnd, c);
    for(int r =0; r < rowEnd; r++){
      int transferIndex = getProcessForNode(c, p, n);
      int index = (r*n+c)*2;
      //printf("%d: sendForces - send - (%i:%i) transferIndex: %i\n", rank, r,c,transferIndex);
      if(transferIndex == rank){
        // move element forces[r][c] to forces[c][r] where c > r
        int newR = c - actualBaseRow;
        int newC = r + actualBaseRow;
        //printf("%d: sendForces - send - updating self (%i,%i) -> (%i,%i)\n", rank, r,c,newR, newC);
        forces[(newR * n + newC)*2 + 0] = forces[(r*n+c)*2]*-1;
        forces[(newR * n + newC)*2 + 1] = forces[(r*n+c)*2+1]*-1;
        // forces[(newR * n + newC)*2 + 0] = tempArr[(r*n+c)*2]*-1;
        // forces[(newR * n + newC)*2 + 1] = tempArr[(r*n+c)*2+1]*-1;
        continue;
      }
      // printf("%d: - sendForces - going to send the value for (%i,%i) from: %i value:(%f,%f)\n", rank, r,c,transferIndex, forces[index], forces[index+1]);
      MPI_Send(&forces[index], 2, MPI_DOUBLE, transferIndex, 2, MPI_COMM_WORLD);
      //MPI_Send(&tempArr[index], 2, MPI_DOUBLE, transferIndex, 2, MPI_COMM_WORLD);
      //printf("%d: updateForce - send - (%i,%i) actualRow: [%i] to %i got values:(%f,%f) index:%i\n", rank, r,c, actualRow, transferIndex, forces[index], forces[index+1], index);
    }
  }
}

void saveImage(char* filePrefix, int step, double* pos, int width, int height, int nLight, int nMedium, int nHeavy){
  // char * fileName = (char *)malloc(strlen(filePrefix)+30);
  int size = width*height*3;
  unsigned char* img = (unsigned char *)malloc(sizeof(unsigned char)*size);
  createImage(pos, width, height,nLight,nMedium,nHeavy, img);
  
  char fileName[30];
  strcpy(fileName, filePrefix);
  // printf("file prefix is: %s, fileName:%s\n", filePrefix, fileName);

  char snum[5];
  sprintf(snum, "%05d", step);
  strcat(fileName, snum);
  strcat(fileName, ".bmp");

  saveBMP(fileName, img, width, height);
  free(img);
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
  int numParticles = nLight + nMedium + nHeavy;
  // printf("width: %i, height:%i\n", width, height);
  printf("light:%d med:%d heavy:%d\n", nLight, nMedium, nHeavy);


  double *positions, *velocities;
  // positions = (double*)malloc(sizeof(double)*3*numParticles);
  //double positions[6] = {500,300,1, 500,500,pow(10,17)};
  //double velocities[4] = {0,0, 0,0};

  //root node stuff goes here
  if(my_rank == 0){
    // set seed for random number generation
    srand48(time(NULL));

    // printf("%d: my_rank is \n", my_rank);
    positions = genRandomPos(width, height, nLight, nMedium, nHeavy);
    velocities = genRandomVel(nLight, nMedium, nHeavy);
    //positions = genTestPos2();
    //velocities = genVelZeros(2);

    // positions = genTestPos2();
    // velocities = genVelZeros(2);

    // positions = genTestArr(numParticles*3, 0,3);
    // velocities = genTestArr(numParticles*2, 1,2);

    // printf("Printing positions [mass]\n");
    // for (int i = 0; i < numParticles; ++i) {
    //   printf("(p):(%.4f,%.4f)[%f]\n", positions[i*3], positions[i*3+1], positions[i*3+2]);
    // }

    // printf("Printing velocities\n");
    // for (int i = 0; i < numParticles; ++i) {
    //   printf("(v):(%.4f,%.4f)\n", velocities[i*2], velocities[i*2+1]);
    // }

    saveImage(filePrefix, 0, positions, width, height, nLight, nMedium, nHeavy);
  }

  int maxBlockSize = (int) ceil(numParticles / (double)p);
  //int l_size = size / p;
  /* double *l_pos = (double*) malloc(sizeof(double) * l_size); */
  /* double *l_vel = (double*) malloc(sizeof(double) * l_size); */
  double *l_pos = (double*) malloc(sizeof(double) * maxBlockSize * 3);
  double *l_vel = (double*) malloc(sizeof(double) * maxBlockSize * 2);
  double * iterPos = l_pos;
  double * emptArr = (double *)malloc(sizeof(double) * maxBlockSize * 3);
  // printf("&lpos:%p, iterpos:%p blockSize:%i\n", l_pos, iterPos, maxBlockSize);

  int *rowDistributions;
  rowDistributions = scatter(positions, l_pos, numParticles, 3, p, my_rank);
  // following call should have same result as call to scatter above
  scatter(velocities, l_vel, numParticles, 2, p, my_rank);

  // init stuff
  //int maxBlockSize = numParticles/p+1;
  // printf("numParticles:%i, maxBlockSize:%i\n", numParticles, maxBlockSize);
  double * forces = (double *)malloc(sizeof(double) * numParticles * maxBlockSize*2);
  int blockSize = rowDistributions[my_rank];
  int recIndex = (my_rank-1+p) % p;
  int sendIndex = (my_rank+1) % p;

  /*
  //test send forces
  int num =0;
  for(int i=0; i < maxBlockSize; i++){
    for(int j =0; j < numParticles; j++){
      forces[(i*numParticles+j)*2] = -1;
      forces[(i*numParticles+j)*2+1] = -1;
      num++;
    }
  }
  sendForces(forces, my_rank, p, numParticles, blockSize);
  printf("going to print forces array\n");
  printVecArr(forces, blockSize*numParticles*2, my_rank, numParticles*2);
  */

  // printf("%d: l_pos is \n", my_rank);
  // printVecArr(l_pos, blockSize*3, my_rank, 3);

  for(int step = 1; step < nSteps; step++){
     // printf("%d: starting step:%i\n", my_rank, step);
    for(int substep = 0; substep < subSteps; substep++){
      iterPos = l_pos;
      // printf("  %d: starting substep:%i step:%i\n", my_rank, substep, step);
      int rankOther = my_rank;
      for(int iter = 0; iter < p; iter++){
        // printf("    %d: starting substep:%i step:%i iter:%i\n", my_rank, substep, step, iter);
        // printf("%d: going to process\n", my_rank);
        //printVecArr(iterPos, maxBlockSize*3, my_rank, 3);
        //calculate(l_pos, iterPos, forces, rowDistributions[my_rank], my_rank);
        calculate(l_pos, iterPos, forces, rowDistributions[my_rank], my_rank, numParticles, rankOther, p);
        // printf("%d: after calculating print forces\n", my_rank);
        // printVecArr(forces, numParticles*maxBlockSize*2, my_rank,numParticles*2);
        // send
        if(sendIndex != my_rank && iter != p-1){
          rankOther = recIndex;
          if(my_rank %2==0){
            // recieve
            // printf("      %d: waiting to recieve from: %d\n", my_rank, recIndex);
            iterPos = emptArr;
            MPI_Recv(iterPos, maxBlockSize*3, MPI_DOUBLE, recIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("      %d: recieved the array\n", my_rank);

            // printf("      %d: going to send to: %d sending\n", my_rank, sendIndex);
            //printVecArr(iterPos, maxBlockSize*3, my_rank, 3);
            MPI_Send(iterPos, maxBlockSize*3, MPI_DOUBLE, sendIndex, 1, MPI_COMM_WORLD);

            //printVecArr(iterPos, maxBlockSize*3, my_rank, 3);
          }else{
            // printf("      %d: going to send to: %d\n", my_rank, sendIndex);
            //printVecArr(iterPos, maxBlockSize*3, my_rank, 3);
            MPI_Send(iterPos, maxBlockSize*3, MPI_DOUBLE, sendIndex, 1, MPI_COMM_WORLD);

            // recieve
            // printf("      %d: waiting to recieve from: %d\n", my_rank, recIndex);
            iterPos = emptArr;
            MPI_Recv(iterPos, maxBlockSize*3, MPI_DOUBLE, recIndex, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("      %d: recieved the array\n", my_rank);
            //printVecArr(iterPos, maxBlockSize*3, my_rank, 3);
          }
        }
      }
      // printf("%d: after done the iteration print forces\n", my_rank);
      // printVecArr(forces, numParticles*maxBlockSize*2, my_rank,numParticles*2);

      // printf("%d: starting sendForces\n", my_rank);
      sendForces(forces, my_rank, p, numParticles, blockSize);
      // printf("%d: done sending forces\n", my_rank);

      // printf("%d: printing forces after send\n", my_rank);
      // printVecArr(forces, numParticles*maxBlockSize*2, my_rank,numParticles*2);

      // printf("%d: going to updatePosition\n", my_rank);
      updatePos(forces, l_pos, l_vel, width, height, numParticles, blockSize);
      // printf("%d: done updatePos \n", my_rank);
      // printf("%d: done iteration\n", my_rank);
      //printf("%d: printing l_pos\n", my_rank);
      //printVecArr(l_pos, maxBlockSize*3, my_rank, 3);
    }

    // printf("%d: starting the gather\n", my_rank);
    gather(positions, l_pos, numParticles, 3, p, my_rank);

    // printf("%d: ending the gather\n", my_rank);
    if(my_rank == 0){
      // if (step == nSteps-1) {
      //   printf("Printing positions\n");
      //   for (int i = 0; i < blockSize; ++i) {
      //     printf("%d:(%.4f,%.4f)\n", my_rank, l_pos[i*3], l_pos[i*3+1]);
      //   }
      // }
      // send all the positions to process 0, do a gather
      // printf("post gather positions:\n");
      // printVecArr(positions, numParticles*3, my_rank,3);
      saveImage(filePrefix, step, positions, width, height, nLight,nMedium,nHeavy);
    }
  }

  MPI_Finalize();
  return 0;
}

