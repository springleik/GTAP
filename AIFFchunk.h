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

using namespace std;

//---------------------------
// AIFF file format chunks
// note that all chunks have even length, conforms to standards
struct chunkAIFFchunk
{
  // data members
  char chunkID[4];    // 4 bytes
  int chunkSize;      // 4 bytes

  // method members
  chunkAIFFchunk (void);
  chunkAIFFchunk (istream&);
  void
  showDetails (ostream&);
};

// total 12 bytes
struct formatAIFFchunk : chunkAIFFchunk
{
  // data members
  char formType[4];     // 4 bytes

  // method members
  formatAIFFchunk (void);
  formatAIFFchunk (int);
  formatAIFFchunk (istream&);
  void
  showDetails (ostream&);
};

// total 26 bytes
#pragma pack(2)
struct commonAIFFchunk : chunkAIFFchunk
{
  // data members
  short channelCount;   // 2 bytes
  int frameCount;       // 4 bytes
  short sampleSize;     // 2 bytes
  char sampRate[10];    // 10 bytes

  // method members
  commonAIFFchunk (void);
  commonAIFFchunk (int, double);
  commonAIFFchunk (istream&);
  void
  showDetails (ostream&);
};

// total 16 bytes
struct soundAIFFchunk : chunkAIFFchunk
{
  // data members
  int offset;       // 4 bytes
  int blockSize;    // 4 bytes

  // method members
  soundAIFFchunk (void);
  soundAIFFchunk (int);
  soundAIFFchunk (istream&);
  void
  showDetails (ostream&);
};

// NOTE swap is only needed for AIFF files on Intel machines
template<typename T>
  void
  byteSwap (T &u)
  {
    T t = u;
    int i = sizeof(u), j = 0;
    char *p = (char*) &u, *q = (char*) &t;
    while (i)
      {
	p[--i] = q[j++];
      }
  }

#endif /* defined(__gtap__AIFFchunk__) */
