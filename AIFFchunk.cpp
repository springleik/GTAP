//
//  AIFFchunk.cpp
//  gtap
//
//  Created by Mark Williamsen on 10/25/13.
//  Copyright (c) 2013 Williamsonic. All rights reserved.
//

#include <cstring>

#include "AIFFchunk.h"

chunkAIFFchunk::chunkAIFFchunk (void)
{
  memset (this, 0, sizeof(*this));
}
formatAIFFchunk::formatAIFFchunk (void)
{
  memset (formType, 0, sizeof(formType));
}
commonAIFFchunk::commonAIFFchunk (void) :
    channelCount (0), frameCount (0), sampleSize (0)
{
  memset (sampRate, 0, sizeof(sampRate));
}
soundAIFFchunk::soundAIFFchunk (void) :
    offset (0), blockSize (0)
{
}

chunkAIFFchunk::chunkAIFFchunk (istream &inp)
{
  memset (this, 0, sizeof(*this));
  inp.read ((char*) this, sizeof(*this));
  byteSwap (chunkSize);
}

formatAIFFchunk::formatAIFFchunk (istream &inp)
{
  inp.read ((char*) this, sizeof(*this));
  byteSwap (chunkSize);
}

commonAIFFchunk::commonAIFFchunk (istream &inp)
{
  inp.read ((char*) this, sizeof(*this));
  byteSwap (chunkSize);
  byteSwap (channelCount);
  byteSwap (frameCount);
  byteSwap (sampleSize);
  cout << "common chunkSize: " << chunkSize << endl;
}

soundAIFFchunk::soundAIFFchunk (istream &inp)
{
  inp.read ((char*) this, sizeof(*this));
  byteSwap (chunkSize);
  byteSwap (offset);
  byteSwap (blockSize);
  cout << "sound chunkSize: " << chunkSize << endl;
}

formatAIFFchunk::formatAIFFchunk (int theSize)
{
  // populate fields
  memcpy (chunkID, "FORM", 4);   // not a null-terminated string
  chunkSize = theSize + 46;
  memcpy (formType, "AIFF", 4);  // not a null-terminated string

  // swap bytes as needed
  byteSwap (chunkSize);
}

commonAIFFchunk::commonAIFFchunk (int theSize, double theRate)
{
  // assume always 16 bit samples, 2 channel stereo
  memcpy (chunkID, "COMM", 4);   // not a null-terminated string
  chunkSize = 18;
  channelCount = 2;
  frameCount = theSize / 4;
  sampleSize = 16;
  union
  {
    long double ldRate;
    char cRate[10];
  } rateUnion =
    { theRate };

  // swap bytes as needed
  int i = 10;
  while (i)
    {
      --i;
      sampRate[i] = rateUnion.cRate[9 - i];
    }
  byteSwap (chunkSize);
  byteSwap (channelCount);
  byteSwap (frameCount);
  byteSwap (sampleSize);
}

soundAIFFchunk::soundAIFFchunk (int theSize)
{
  // populate fields
  memcpy (chunkID, "SSND", 4); // not a null-terminated string
  chunkSize = theSize + 8;
  offset = 0;
  blockSize = 0;

  // swap bytes as needed
  byteSwap (chunkSize);
  byteSwap (offset);
  byteSwap (blockSize);
}

void
chunkAIFFchunk::showDetails (ostream &out)
{
  // copy 4 character ID to an array
  char id[] = "ABCD";
  memcpy (id, chunkID, 4);
  out << "chunkID," << id << endl;
  out << "chunkSize," << chunkSize << endl;
}

void
formatAIFFchunk::showDetails (ostream &out)
{
  chunkAIFFchunk::showDetails (out);
  char id[] = "ABCD";
  memcpy (id, formType, 4);
  out << "formType," << id << endl;
}

void
commonAIFFchunk::showDetails (ostream &out)
{
  chunkAIFFchunk::showDetails (out);
  out << "channelCount," << channelCount << endl;
  out << "frameCount," << frameCount << endl;
  out << "sampleSize," << sampleSize << endl;

  // swap bytes as needed
  union
  {
    long double ldRate;
    char cRate[10];
  } rateUnion =
    { 0 };
  int i = 10;
  while (i)
    {
      --i;
      rateUnion.cRate[i] = sampRate[9 - i];
    }
  out << "sampRate," << rateUnion.ldRate << endl;
}

void
soundAIFFchunk::showDetails (ostream &out)
{
  chunkAIFFchunk::showDetails (out);
  out << "offset," << offset << endl;
  out << "blockSize," << blockSize << endl;
}

