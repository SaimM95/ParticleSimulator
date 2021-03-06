#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
// double velocityLightMin = 11;
// double velocityLightMax = 15;
double velocityLightMin = 0.8;
double velocityLightMax = 1.2;

double velocityMediumMin = 0.01;
double velocityMediumMax = 0.05;

//heavy particles are the slowest
double velocityHeavyMin = 0.001;
double velocityHeavyMax = 0.005;

// double velocityMediumMin = 0.01;
// double velocityMediumMax = 0.1;
//
// //heavy particles are the slowest
// double velocityHeavyMin = 0.001;
// double velocityHeavyMax = 0.005;

//mass
double massLightMin = 10000;
double massLightMax = 50000;

double massMediumMin = 600000;
double massMediumMax = 1000000;

double massHeavyMin = 11000000;
double massHeavyMax = 15000000;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif
