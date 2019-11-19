//
//  WAVEchunk.h
//  gtap
//
//  Created by Mark Williamsen on 10/25/13.
//  Copyright (c) 2013 Williamsonic. All rights reserved.
//

#ifndef __gtap__WAVEchunk__
#define __gtap__WAVEchunk__

#include <iostream>
using namespace std;

//----------------------------
// WAVE file format chunks
struct chunkWAVEchunk
{
	// data members
	char chunkID[4];	// 4 characters, not null-terminated
	int  chunkSize;		// remaining byte count in this file
    
    // method members
    chunkWAVEchunk(void);
    chunkWAVEchunk(istream &);
    void showDetails(ostream &);
};

// total 12 bytes
struct riffWAVEchunk : chunkWAVEchunk
{
	// data member
	char format[4];		// 4 characters, not null-terminated
	
    // method members
    riffWAVEchunk(void);
    riffWAVEchunk(int);
    riffWAVEchunk(istream &);
    void showDetails(ostream &);
};

// total 24 bytes
struct formatWAVEchunk : chunkWAVEchunk
{
	// data members
	short fmtCode;		// data format code
	short numChan;		// number of audio channels
	int   sampRate;     // sample rate per second
	int   byteRate;     // byte rate per second
	short blockAlign;	// byte count per sample
	short bitsSamp;		// bits count per sample
	
    // method members
    formatWAVEchunk(void);
    formatWAVEchunk(int);
    formatWAVEchunk(istream &);
    void showDetails(ostream &);
};

// total 8 bytes
struct dataWAVEchunk : chunkWAVEchunk
{
	// actual sound data starts here, but we won't try to load it into memory
	
    // method members
    dataWAVEchunk(void);
    dataWAVEchunk(int);
    dataWAVEchunk(istream &);
    void showDetails(ostream &);
};

#endif /* defined(__gtap__WAVEchunk__) */
