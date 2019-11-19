//
//  AIFFchunk.h
//  gtap
//
//  Created by Mark Williamsen on 10/25/13.
//  Copyright (c) 2013 Williamsonic. All rights reserved.
//

#ifndef __gtap__AIFFchunk__
#define __gtap__AIFFchunk__

#include <iostream>
#include <CoreServices/CoreServices.h>

using namespace std;

//---------------------------
// AIFF file format chunks
struct chunkAIFFchunk
{
    // data members
    char chunkID[4];    // 4 bytes
    int chunkSize;      // 4 bytes
    
    // method members
    chunkAIFFchunk(void);
    chunkAIFFchunk(istream &);
    void showDetails(ostream &);
};

// total 12 bytes
struct formatAIFFchunk : chunkAIFFchunk
{
    // data members
    char formType[4];     // 4 bytes
    
    // method members
    formatAIFFchunk(void);
    formatAIFFchunk(int);
    formatAIFFchunk(istream &);
    void showDetails(ostream &);
};

// total 26 bytes
struct commonAIFFchunk : chunkAIFFchunk
{
    // data members
    short channelCount;   // 2 bytes
    int   frameCount;     // 4 bytes
    short sampleSize;     // 2 bytes
    Float80 sampleRate;   // 10 bytes
    
    // method members
    commonAIFFchunk(void);
    commonAIFFchunk(int, double);
    commonAIFFchunk(istream &);
    void showDetails(ostream &);
} __attribute__ ((packed));

// total 16 bytes
struct soundAIFFchunk : chunkAIFFchunk
{
    // data members
    int offset;       // 4 bytes
    int blockSize;    // 4 bytes
    
    // method members
    soundAIFFchunk(void);
    soundAIFFchunk(int);
    soundAIFFchunk(istream &);
    void showDetails(ostream &);
};

// NOTE swap is only needed for AIFF files on Intel machines
template<typename T> void byteSwap(T &u)
{
    T t     = u;
    int   i = sizeof(u), j = 0;
    char *p = (char *)&u, *q = (char *)&t;
    while (i) {p[--i] = q[j++];}
}

#endif /* defined(__gtap__AIFFchunk__) */
