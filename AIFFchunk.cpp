//
//  AIFFchunk.cpp
//  gtap
//
//  Created by Mark Williamsen on 10/25/13.
//  Copyright (c) 2013 Williamsonic. All rights reserved.
//

#include "AIFFchunk.h"

chunkAIFFchunk::chunkAIFFchunk(void)   {memset(this, 0, sizeof(chunkAIFFchunk));}
formatAIFFchunk::formatAIFFchunk(void) {memset(this, 0, sizeof(formatAIFFchunk));}
commonAIFFchunk::commonAIFFchunk(void) {memset(this, 0, sizeof(commonAIFFchunk));}
soundAIFFchunk::soundAIFFchunk(void)   {memset(this, 0, sizeof(soundAIFFchunk));}

chunkAIFFchunk::chunkAIFFchunk(istream &inp)
{
    memset(this, 0, sizeof(chunkAIFFchunk));
    inp.read((char *)this, sizeof(chunkAIFFchunk));
    byteSwap(chunkSize);
}

formatAIFFchunk::formatAIFFchunk(istream &inp)
{
    memset(this, 0, sizeof(formatAIFFchunk));
    inp.read((char *)this, sizeof(formatAIFFchunk));
    byteSwap(chunkSize);
}

commonAIFFchunk::commonAIFFchunk(istream &inp)
{
    memset(this, 0, sizeof(commonAIFFchunk));
    inp.read((char *)this, sizeof(commonAIFFchunk));
    byteSwap(chunkSize);
    byteSwap(channelCount);
    byteSwap(frameCount);
    byteSwap(sampleSize);
}

soundAIFFchunk::soundAIFFchunk(istream &inp)
{
    memset(this, 0, sizeof(soundAIFFchunk));
    inp.read((char *)this, sizeof(soundAIFFchunk));
    byteSwap(chunkSize);
    byteSwap(offset);
    byteSwap(blockSize);
}

formatAIFFchunk::formatAIFFchunk(int theSize)
{
    // populate fields
    memcpy(chunkID, "FORM", 4);   // not a null-terminated string
    chunkSize = theSize + 46;
    memcpy(formType, "AIFF", 4);  // not a null-terminated string
    
    // swap bytes as needed
    byteSwap(chunkSize);
}

commonAIFFchunk::commonAIFFchunk(int theSize, double theRate)
{
    // assume always 16 bit samples, 2 channel stereo
    memcpy(chunkID, "COMM", 4);   // not a null-terminated string
    chunkSize = 18;
    channelCount = 2;
    frameCount = theSize / 4;
    sampleSize = 16;
    dtox80(&theRate, &sampleRate);
    
    // swap bytes as needed
    byteSwap(chunkSize);
    byteSwap(channelCount);
    byteSwap(frameCount);
    byteSwap(sampleSize);
}

soundAIFFchunk::soundAIFFchunk(int theSize)
{
    // populate fields
    memcpy(chunkID, "SSND", 4); // not a null-terminated string
    chunkSize = theSize + 8;
    offset = 0;
    blockSize = 0;
    
    // swap bytes as needed
    byteSwap(chunkSize);
    byteSwap(offset);
    byteSwap(blockSize);
}

void chunkAIFFchunk::showDetails(ostream &out)
{
    // copy 4 character ID to an array
    char id[] = "ABCD";
    memcpy(id, chunkID, 4);
    out << "chunkID," << id << endl;
    out << "chunkSize," << chunkSize << endl;
}

void formatAIFFchunk::showDetails(ostream &out)
{
    chunkAIFFchunk::showDetails(out);
    char id[] = "ABCD";
    memcpy(id, formType, 4);
    out << "formType," << id << endl;
}

void commonAIFFchunk::showDetails(ostream &out)
{
    chunkAIFFchunk::showDetails(out);
    out << "channelCount," << channelCount << endl;
    out << "frameCount," << frameCount << endl;
    out << "sampleSize," << sampleSize << endl;
    out << "sampleRate," << x80tod(&sampleRate) << endl;
}

void soundAIFFchunk::showDetails(ostream &out)
{
    chunkAIFFchunk::showDetails(out);
    out << "offset," << offset << endl;
    out << "blockSize," << blockSize << endl;
}

