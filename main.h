//
//  main.h
//  gtap
//
//  Created by Mark Williamsen on 1/26/14.
//  Copyright (c) 2014 Williamsonic. All rights reserved.
//
// http://www.williamsonic.com/DipoleMic/
//

#ifndef gtap_main_h
#define gtap_main_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
#include <cmath>
#include <map>

#include "AIFFchunk.h"
#include "WAVEchunk.h"

using namespace std;

// enumerations
enum burstType
{
  simpleSine, harmonicSine, simpleCosine, harmonicCosine, unknownBurst
};
enum fileType
{
  typeAIFF, typeWAVE, unknownFile
};

// methods
void
showParams (ostream&);
void
showHelp (ostream&);
void
reset (void);
void
writeOutputFile (istream&);
void
readInputFile (istream&);
void
readPilotFile (istream&);

extern const int sampleRate;

#endif
