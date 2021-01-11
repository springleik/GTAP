//
//  main.cpp
//  gtap
//
//  Created by Mark Williamsen on 1/26/14.
//  Copyright (c) 2014 Williamsonic. All rights reserved.
//

#include <cstring>
#include <cassert>

#include "gtap.h"

// data members
burstType theBurstType = harmonicSine;  // waveform type
fileType  theFileType = typeAIFF;       // file format type
const int sampleRate = 44100;   // samples per second
int averaging = 1;              // number of bursts to average over
int numBursts = 201;            // number of bursts
int numCycles = 1;              // number of cycles per burst
int delay = 22050;              // samples to start of first burst
int burstInterval = 22050;      // total samples per burst
int burstMin = 100;             // minimum number of samples within any burst
int duration = 100;             // actual sample count within burst
int burstCount = 0;             // index to count bursts
int cycleCount = 0;             // index to count cycles
double startFreq = 100.;        // frequency of first tone burst
double stopFreq = 10000.;       // frequency of last tone burst
double leftAmpl = 0.5;          // left channel gain factor relative to unity
double rightAmpl = 0.5;         // right channel gain factor relative to unity
double freqIncrement = 1.0;     // constant factor relating adjacent frequencies
double nominalFreq = 100.;      // nominal tone burst frequency
double actualFreq = 100.;       // actual tone burst frequency
double factor = 100.;           // frequency in sample-based units
bool verbose = true;            // verbose output mode
bool pilot = true;              // pilot tone enabled
bool phase = true;              // output in phase

int main(int argc, const char * argv[])
{
  // check object sizes
  assert(12 == sizeof(formatAIFFchunk));
  assert(26 == sizeof(commonAIFFchunk));
  assert(16 == sizeof(soundAIFFchunk));
  assert(12 == sizeof(riffWAVEchunk));
  assert(24 == sizeof(formatWAVEchunk));
  assert(8  == sizeof(dataWAVEchunk));
  
  // process input commands
  cout <<
  "Gated Toneburst Analysis Program, http://www.williamsonic.com"
  "\nusage: gtap < script.txt > output.csv, built: "  __DATE__ ", " __TIME__ "."
  "\ngtap: ";
  string theLine;
  while (getline(cin, theLine))
  {
    // identify first token
    stringstream str(theLine);
    string cmd, path;
    int delta = 0;
    str >> cmd;
    if (!cmd.length())          {continue;} // ignore blank lines
    else if ("#" == cmd)        {continue;} // ignore comments
    else if ("quit" == cmd)     {break;}    // done with program
    else if ("exit" == cmd)     {break;}    // done with program
    
    // change a setting, no arguments expected
    else if ("aif" == cmd)      {theFileType = typeAIFF;}
    else if ("aiff" == cmd)     {theFileType = typeAIFF;}
    else if ("wav" == cmd)      {theFileType = typeWAVE;}
    else if ("wave" == cmd)     {theFileType = typeWAVE;}
    else if ("sin1" == cmd)     {theBurstType = simpleSine;}
    else if ("sine1" == cmd)    {theBurstType = simpleSine;}
    else if ("cos1" == cmd)     {theBurstType = simpleCosine;}
    else if ("cosine1" == cmd)  {theBurstType = simpleCosine;}
    else if ("sin2" == cmd)     {theBurstType = harmonicSine;}
    else if ("sine2" == cmd)    {theBurstType = harmonicSine;}
    else if ("cos2" == cmd)     {theBurstType = harmonicCosine;}
    else if ("cosine2" == cmd)  {theBurstType = harmonicCosine;}
    
    // change a setting, expect one numeric argument
    else if ("avg" == cmd)      {str >> averaging;}
    else if ("average" == cmd)  {str >> averaging;}
    else if ("count" == cmd)   {str >> numBursts;}
    else if ("cycles" == cmd)   {str >> numCycles;}
    else if ("delay" == cmd)    {str >> delay;}
    else if ("int" == cmd)      {str >> burstInterval;}
    else if ("interval" == cmd) {str >> burstInterval;}
    else if ("min" == cmd)      {str >> burstMin;}
    else if ("minimum" == cmd)  {str >> burstMin;}
    else if ("start" == cmd)    {str >> startFreq;}
    else if ("stop" == cmd)     {str >> stopFreq;}
    else if ("left" == cmd)     {str >> leftAmpl;}
    else if ("right" == cmd)    {str >> rightAmpl;}
    else if ("ampl" == cmd)     {str >> leftAmpl; rightAmpl = leftAmpl;}
    else if ("amplitude" == cmd){str >> leftAmpl; rightAmpl = leftAmpl;}
    else if ("verbose" == cmd)  {str >> verbose;}
    else if ("verb" == cmd)     {str >> verbose;}
    else if ("pilot" == cmd)    {str >> pilot;}
    else if ("phase" == cmd)    {str >> phase;}
    
    // take an action, may consume additional arguments
    else if ("show" == cmd)     {reset(); showParams(cout);}
    else if ("help" == cmd)     {showHelp(cout);}
    else if ("write" == cmd)    {writeOutputFile(str);}
    else if ("read" == cmd)     {readInputFile(str);}
    else if ("find" == cmd)     {readPilotFile(str);}
    else if ("delta" == cmd)    {str >> delta; delay += delta;}
    
    else {cout << "Unknown command." << endl;}
    cout << "gtap: ";
  }
  return 0;
}

string getText(burstType b)
{
  switch(b)
  {
    case simpleSine:     return "simpleSine";
    case simpleCosine:   return "simpleCosine";
    case harmonicSine:   return "harmonicSine";
    case harmonicCosine: return "harmonicCosine";
    default:             return "Unknown burst type";
  }
}

string getText(fileType f)
{
  switch(f)
  {
    case typeAIFF: return "typeAIFF";
    case typeWAVE: return "typeWAVE";
    default:       return "Unknown file type";
  }
}

void showParams(ostream &out)
{
  out << "User settings:"
    << "\n delay,         " << delay
    << "\n burstInterval, " << burstInterval
    << "\n startFreq,     " << startFreq
    << "\n stopFreq,      " << stopFreq
    << "\n averaging,     " << averaging
    << "\n numBursts,     " << numBursts
    << "\n numCycles,     " << numCycles
    << "\n leftAmpl,      " << leftAmpl
    << "\n rightAmpl,     " << rightAmpl
    << "\n theBurstType,  " << getText(theBurstType)
    << "\n theFileType,   " << getText(theFileType)
    << "\n burstMin,      " << burstMin
    << "\n phase,         " << (phase? "true": "false")
    << "\n pilot,         " << (pilot? "true": "false")
    << "\n verbose,       " << (verbose? "true": "false") << endl;
  
  out << "Not settings:"
    << "\n sampleRate,    " << sampleRate
    << "\n duration,      " << duration
    << "\n nominalFreq,   " << nominalFreq
    << "\n actualFreq,    " << actualFreq
    << "\n freqIncrement, " << freqIncrement
    << "\n burstCount,    " << burstCount
    << "\n factor,        " << factor << endl;
}

void showHelp(ostream &out)
{
  out <<
  "Available commands:"
  "\n aiff      (aiff file format)"
  "\n amplitude gain factor"
  "\n average   bursts to average"
  "\n cos1      (simple cosine)"
  "\n cos2      (harmonic cosine)"
  "\n count     number of bursts"
  "\n cycles    minimum cycles per burst"
  "\n delay     samples to start of first burst"
  "\n delta     samples to add to delay"
  "\n find      name of input file, left|right"
  "\n help      (this list)"
  "\n interval  total samples per interval"
  "\n left      left channel gain"
  "\n minimum   minimum samples per burst"
  "\n phase     0|1 (output in phase)"
  "\n pilot     0|1 (pilot tone enabled)"
  "\n read      name of input file, name of log file"
  "\n right     right channel gain"
  "\n show      (measurement setup)"
  "\n sin1      (simple sine"
  "\n sin2      (harmonic sine)"
  "\n start     first frequency (Hz)"
  "\n stop      last frequency (Hz)"
  "\n verbose   0|1 (verbose output mode)"
  "\n wave      (wave file format)"
  "\n write     name of output file, name of log file"
  << endl;
}

// update data members for a particular frequency
void calc(void)
{
  // find least number of full cycles whose duration exceeds burst minimum
  while ((sampleRate / nominalFreq) * cycleCount < burstMin) {cycleCount++;}

  // round off (truncate) duration to exact samples
  duration = int((sampleRate / nominalFreq) * cycleCount);

  // find frequency that exactly fills the duration
  actualFreq = 1.0 * sampleRate * cycleCount / duration;

  // calculate common factor for sine and cosine
  factor = 2.0 * M_PI * actualFreq / sampleRate;
}

// read tone burst from disk, matched filter technique
// looks only for the exact frequency being measured
void readInputBurst(istream &inFile, ostream &logStream)
{
  int i = 0, j = 0;		// local loop indices
  short a = 0, b = 0;
  complex<double> sum1(0,0);          // channel 1 response
  complex<double> sum2(0,0);          // channel 2 response
  complex<double> sum3(0,0);          // channel 1 background
  complex<double> sum4(0,0);          // channel 2 background
  complex<double> cfactor(0, factor); // complex frequency in rad/sec
  complex<double> shift(0, -1.0);     // 90 deg. phase shift

  // iterate over averaging, burst interval
  for (i = 0; i < averaging; i++)
    for (j = 0; j < burstInterval; j++)
    {
      // read data from both channels
      inFile.read((char *)&a, 2);
      inFile.read((char *)&b, 2);
      
      // swap bytes for AIFF files
      if (typeAIFF == theFileType) {byteSwap(&a); byteSwap(&b);}
      
      // analyze burst response
      // this is a single frequency discrete Fourier transform
      // phase angle is referred to start of burst
      if (j < duration)
      {
        complex<double> ccoeff = exp(cfactor * double(j)) * shift;
        sum1 += double(a) * ccoeff;
        sum2 += double(b) * ccoeff;
      }
      
      // analyze background level, near end of burst burstInterval
      if ((j < (burstInterval - duration)) && (j >= (burstInterval -  (2 * duration))))
      {
        complex<double> ccoeff = exp(cfactor * double(j)) * shift;
        sum3 += double(a) * ccoeff;
        sum4 += double(b) * ccoeff;
      }
    }

  // factor out sample count and averaging, normalize to +0 dB
  sum1 /= (duration * averaging * leftAmpl  * 32767.0 / 2.0);
  sum2 /= (duration * averaging * rightAmpl * 32767.0 / 2.0);
  sum3 /= (duration * averaging * leftAmpl  * 32767.0 / 2.0);
  sum4 /= (duration * averaging * rightAmpl * 32767.0 / 2.0);

  // report results to console
  // 0.0 dB reference level when analyzing original generated file
  logStream << ',' << abs(sum1)          // magnitude channel 1
  << ',' << abs(sum2)               // magnitude channel 2
  << ',' << 20.0*log10(abs(sum1))   // dB channel 1
  << ',' << 20.0*log10(abs(sum2))   // dB channel 2
  << ',' << 20.0*log10(abs(sum1)/abs(sum2))	// dB difference
  << ',' << arg(sum1)               // phase channel 1
  << ',' << arg(sum2)               // phase channel 2
  << ',' << arg(sum1) -arg(sum2)		// phase difference
  << ',' << 20.0*log10(abs(sum3))		// dB background 1
  << ',' << 20.0*log10(abs(sum4))		// dB background 2
  << endl;
  return;
}

void reset(void)
{
  burstCount = numBursts;
  cycleCount = numCycles;
  nominalFreq = startFreq;
  freqIncrement = pow(stopFreq/startFreq, 1.0/(numBursts - 1));
  calc();
}

// tone burst iterator
void next(void)
{
  // only needed for freq sweep case
  nominalFreq *= freqIncrement;
  calc();
  burstCount--;
  return;
}

// test to see if we are done
bool good(void) {return (burstCount > 0);}

// true if signs of a and b differ
bool signDiff(int a, int b) {return (0 > int(a ^ b));}

// scan sound data to find 441 Hz pilot tone
void readPilotData(istream &inFile, ostream &logStream)
{
  // skip over nominal delay time
  inFile.seekg(2 * 2 * delay, inFile.cur);
  
  // allocate and clear sample FIFOs and axis crossing arrays
  short a[100], b[100];
  memset(a, 0, sizeof(a));
  memset(b, 0, sizeof(b));
  enum {I1, Q1, I2, Q2, NChan};
  int next[NChan] = {0, 0, 0, 0},   prev[NChan] = {0, 0, 0, 0};
  
  // allocate and clear axis-crossing FIFOs for each channel
  struct xing{int samp; int ampl;} xingI1 = {0, 0},
  xingQ1 = {0, 0}, xingI2 = {0, 0}, xingQ2 = {0, 0};
  
  // look at 10 seconds of data
  ofstream outFile("outPut.csv");     // overwrite temp file each time through
  outFile << "n,I1,Q1,I2,Q2" << endl; // column headings
  bool keep = false;
  int offset = 0;
  for (int i = 0; i < (10 * sampleRate); i++)
  {
    // slide fifo along input stream
    int j = 100 - 1;
    while(j--)
    {
      a[j+1] = a[j];
      b[j+1] = b[j];
    }
    
    // read data from both channels
    inFile.read((char *)a, 2); inFile.read((char *)b, 2);
    
    // swap bytes for AIFF files
    if (typeAIFF == theFileType) {byteSwap(a); byteSwap(b);}
    
    // keep previous sums, update next
    prev[I1] = next[I1]; next[I1] = prev[I1]+a[0]-a[49]-a[50]+a[99];
    prev[I2] = next[I2]; next[I2] = prev[I2]+b[0]-b[49]-b[50]+b[99];
    prev[Q1] = next[Q1]; next[Q1] = prev[Q1]+a[0]-a[24]-a[25]+a[74]+a[75]-a[99];
    prev[Q2] = next[Q2]; next[Q2] = prev[Q2]+b[0]-b[24]-b[25]+b[74]+b[75]-b[99];
    
    // write debug output
    outFile << i << ',' << next[I1] << ',' << next[Q1] << ',' << next[I2] << ',' << next[Q2] << endl;
    
    // capture zero crossings, look for pilot tone signature
    if (signDiff(prev[I1], next[I1]))
    {
      keep = true;
      xingI1.samp = i;
      int diff = xingI1.samp - xingQ1.samp;
      if (diff < 18 || diff > 32) {keep = false;}
      double ratio = log(fabs(double(next[Q1]) / double(xingI1.ampl)));
      if (ratio < 0.5) {keep = false;}
      int val = xingI1.ampl = next[Q1];
      if (abs(val) < 10000) {keep = false;}
      if (keep)
      {
        // exit on first match
        cout << i << ", I1, " << val << ", " << diff << ", " << ratio << endl;
        offset = xingI1.samp;
        break;
      }
    }
    if (signDiff(prev[Q1], next[Q1]))
    {
      keep = true;
      xingQ1.samp = i;
      int diff = xingQ1.samp - xingI1.samp;
      if (diff < 18 || diff > 32) {keep = false;}
      double ratio = log(fabs(double(next[I1]) / double(xingQ1.ampl)));
      if (ratio < 0.5) {keep = false;}
      int val = xingQ1.ampl = next[I1];
      if (abs(val) < 10000) {keep = false;}
      if (keep)
      {
        // exit on first match
        cout << i << ", Q1, " << val << ", " << diff << ", " << ratio << endl;
        offset = xingQ1.samp;
        break;
      }
    }
    if (signDiff(prev[I2], next[I2]))
    {
      keep = true;
      xingI2.samp = i;
      int diff = xingI2.samp - xingQ2.samp;
      if (diff < 18 || diff > 32) {keep = false;}
      double ratio = log(fabs(double(next[Q2]) / double(xingI2.ampl)));
      if (ratio < 0.5) {keep = false;}
      int val = xingI2.ampl = next[Q2];
      if (abs(val) < 10000) {keep = false;}
      if (keep)
      {
        // exit on first match
        cout << i << ", I2, " << val << ", " << diff << ", " << ratio << endl;
        offset = xingI2.samp;
        break;
      }
    }
    if (signDiff(prev[Q2], next[Q2]))
    {
      keep = true;
      xingQ2.samp = i;
      int diff = xingQ2.samp - xingI2.samp;
      if (diff < 18 || diff > 32) {keep = false;}
      double ratio = log(fabs(double(next[I2]) / double(xingQ2.ampl)));
      if (ratio < 0.5) {keep = false;}
      int val = xingQ2.ampl = next[I2];
      if (abs(val) < 10000) {keep = false;}
      if (keep)
      {
        // exit on first match
        cout << i << ", Q2, " << val << ", " << diff << ", " << ratio << endl;
        offset = xingQ2.samp;
        break;
      }
    }
  }
  
  // compute total delay including pilot tone
  if (offset)
  {
    delay += (offset - 99);
    logStream << "Pilot tone found at sample," << delay << endl;
  }
}

// show burst parameters on console
void showBurst(ostream &out)
{
  out << cycleCount
  << ',' << duration
  << ',' << nominalFreq
  << ',' << actualFreq;
}

void readInputData(istream &inFile, ostream &logStream)
{
  // skip over pilot tone if enabled
  inFile.seekg(2 * 2 * delay, inFile.cur);
  if (pilot) {inFile.seekg(2 * 2 * burstInterval, inFile.cur);}
  
  // show column headings here
  logStream << "numCyc,duration,nomFreq,actFreq,"
  "abs1,abs2,dB1,dB2,dBDiff,"
  "phase1,phase2,phaseDiff,bkg1,bkg2" << endl;

  // iterate over tone bursts while reading from disk
  for(reset(); good(); next())
  {
    // check input file before reading
    if (inFile.eof())
    {
      logStream << "Failed to read tone bursts from disk." << endl;
      return;
    }
    showBurst(logStream);
    readInputBurst(inFile, logStream);
  }
  logStream << "Next file position," << inFile.tellg() << endl;
}

bool readAIFFheader(istream &inFile, ostream &logStream)
{
  // read AIFF header from disk
  formatAIFFchunk theFormat(inFile);
  if (strncmp(theFormat.chunkID, "FORM", 4)) {return false;}
  if (strncmp(theFormat.formType, "AIFF", 4)) {return false;}
  logStream << "AIFF file detected." << endl;
  
  // obtain a list of chunks and file positions
  map<string, long int> theMap;
  theMap["FORM"] = 0;
  while (inFile)
  {
    long thePos = inFile.tellg();
    chunkAIFFchunk theChunk(inFile);
    if (!inFile) {break;}
    char s[] = "ABCD";
    memcpy(s, theChunk.chunkID, 4);
    theMap[s] = thePos;
    inFile.seekg(theChunk.chunkSize, inFile.cur);
  };
  
  // show list contents on console
  logStream << "Chunk List:" << endl;
  for (auto i = theMap.begin(); i != theMap.end(); i++)
  {
    logStream << ' ' << i->first << ',' << i->second << endl;
  }
  
  // populate required common chunk
  // TODO find() may crash if chunk is missing (4 pl.)
  long chunkPos = 0;
  chunkPos = theMap.find("COMM")->second;
  if (!chunkPos)
  {
    logStream << "Required COMM chunk not found!" << endl;
    return false;
  }
  inFile.clear();
  inFile.seekg(chunkPos, inFile.beg);
  commonAIFFchunk theCommon(inFile);
  if (!inFile)
  {
    logStream << "Failed to read COMM chunk." << endl;
    return false;
  }
  
  // populate required sound chunk
  chunkPos = 0;
  chunkPos = theMap.find("SSND")->second;
  if (!chunkPos)
  {
    logStream << "Required SSND chunk not found!" << endl;
    return false;
  }
  inFile.seekg(chunkPos, inFile.beg);
  soundAIFFchunk theSound(inFile);
  if (!inFile)
  {
    logStream << "Failed to read SSND chunk." << endl;
    return false;
  }
  
  // show results on console
  if (verbose)
  {
    theFormat.showDetails(logStream);
    theCommon.showDetails(logStream);
    theSound.showDetails(logStream);
  }
  
  // leave file pointer on start of sound data
  return true;
}

bool readWAVEheader(istream &inFile, ostream &logStream)
{
  // read WAVE header from disk
  riffWAVEchunk theRiff(inFile);
  if (strncmp(theRiff.chunkID, "RIFF", 4)) {return false;}
  if (strncmp(theRiff.format, "WAVE", 4)) {return false;}
  logStream << "WAVE file detected." << endl;
  
  // obtain a list of chunks and file positions
  map<string, long int> theMap;
  theMap["RIFF"] = 0;
  while (inFile)
  {
    long thePos = inFile.tellg();
    chunkWAVEchunk theChunk(inFile);
    if (!inFile) {break;}
    char s[] = "ABCD";
    memcpy(s, theChunk.chunkID, 4);
    theMap[s] = thePos;
    inFile.seekg(theChunk.chunkSize, inFile.cur);
  };
  
  // show list contents on console
  logStream << "Chunk List:" << endl;
  for (auto i = theMap.begin(); i != theMap.end(); i++)
  {
    logStream << ' ' << i->first << ',' << i->second << endl;
  }
  
  // populate required common chunk
  // TODO find() may crash if chunk is missing (4 pl.)
  long chunkPos = 0;
  chunkPos = theMap.find("fmt ")->second;
  if (!chunkPos)
  {
    logStream << "Required fmt chunk not found!" << endl;
    return false;
  }
  inFile.clear();
  inFile.seekg(chunkPos, inFile.beg);
  formatWAVEchunk theFormat(inFile);
  if (!inFile)
  {
    logStream << "Failed to read fmt chunk." << endl;
    return false;
  }
  
  // populate required sound chunk
  chunkPos = 0;
  chunkPos = theMap.find("data")->second;
  if (!chunkPos)
  {
    logStream << "Required data chunk not found!" << endl;
    return false;
  }
  inFile.seekg(chunkPos, inFile.beg);
  dataWAVEchunk theData(inFile);
  if (!inFile)
  {
    logStream << "Failed to read data chunk." << endl;
    return false;
  }
  
  // show results on console
  if (verbose)
  {
    theRiff.showDetails(logStream);
    theFormat.showDetails(logStream);
    theData.showDetails(logStream);
  }
  
  // leave file pointer on start of sound data
  return true;
}

void readInputFile(istream &str)
{
  // get file path and log path from command stream
  string inPath, logPath;
  str >> inPath >> logPath;
  ofstream logFile;
  if (logPath.length()) {logFile.open(logPath);}
  ostream &logStream = logPath.length()? logFile: cout;
  if (!inPath.length())
  {
    logStream << "Required input file name missing." << endl;
    return;
  }
  
  // must be opened in binary mode
  ifstream inFile(inPath, ios::in | ios::binary);
  if (!inFile)
  {
    logStream << "Failed to open input file: " << inPath << endl;
    return;
  }
  logStream << "Reading input file," << inPath << endl;
  
  theFileType = unknownFile;
  if (readAIFFheader(inFile, logStream))
  {
    theFileType = typeAIFF;
    readInputData(inFile, logStream);
    return;
  }
  
  // reset file pointer for next operation
  inFile.clear();
  inFile.seekg(0, inFile.beg);
  if (readWAVEheader(inFile, logStream))
  {
    theFileType = typeWAVE;
    readInputData(inFile,logStream);
    return;
  }
  
  // TODO consider handling AIFC file type
  logStream << "Unexpected input file type!" << endl;
}

void readPilotFile(istream &str)
{
  // get file path and log path from command stream
  string inPath, logPath;
  str >> inPath >> logPath;
  ofstream logFile;
  if (logPath.length()) {logFile.open(logPath);}
  ostream &logStream = logPath.length()? logFile: cout;
  if (!inPath.length())
  {
    logStream << "Required pilot file name missing." << endl;
    return;
  }

  // must be opened in binary mode
  ifstream inFile(inPath, ios::in | ios::binary);
  if (!inFile)
  {
    logStream << "Failed to open pilot file: " << inPath << endl;
    return;
  }
  logStream << "Reading pilot file," << inPath << endl;
  
  theFileType = unknownFile;
  if (readAIFFheader(inFile, logStream))
  {
    theFileType = typeAIFF;
    readPilotData(inFile, logStream);
    return;
  }
  
  // reset file pointer for next operation
  inFile.clear();
  inFile.seekg(0, inFile.beg);
  if (readWAVEheader(inFile, logStream))
  {
    theFileType = typeWAVE;
    readPilotData(inFile,logStream);
    return;
  }
  
  // TODO consider handling AIFC file type
  logStream << "Unexpected pilot file type!" << endl;
}

// write a burst to output stream
void writeOutputData(ostream &outfile)
{
	int i = 0, j = 0;	// local loop indices
	short a = 0, b = 0;
	double y = 0.0;
	
	// iterate over averaging, burst interval, number of channels
	for (i = 0; i < averaging; i++)
    for (j = 0; j < burstInterval; j++)
    {
      // calculate burst waveform
      if (j < duration)
      {
        // compute the selected waveform
        switch(theBurstType)
        {
          case simpleSine:     y = sin(factor * j); break;
          case harmonicSine:   y = sin(factor * j) - sin(2.0 * factor * j) / 2.0; break;
          case simpleCosine:   y = 1.0 - cos(factor * j); break;
          case harmonicCosine: y = cos(2.0 * factor * j) - cos(factor * j); break;
          default: cout << "Unexpected burst type." << endl; return;
        }
        // normalize to +0 dB amplitude, convert to short word
        a = short(y * leftAmpl  * 32767.0);
        b = short(y * rightAmpl * 32767.0);
        if (!phase) {a = -a; b = -b;}
        if (typeAIFF == theFileType) {byteSwap(&a); byteSwap(&b);}
      }
      else {a = b = 0;}   // write silence between bursts
      
      // write left channel data, then right
      outfile.write((char *)&a, 2);
      outfile.write((char *)&b, 2);
    }
	return;
}

// get byte count for sound data in wave file
int getSize(void)
{
	// bytes/sample * num channels * (samples/burst * averaging * num bursts + delay)
	// always assumes 2 byte samples, 2 channel stereo, add one burst for pilot tone
	return (2 * 2 * (burstInterval * (averaging * numBursts + (pilot? 1:0)) + delay));
}

// write single cycle at 441 hz
void writePilotTone(ostream &outfile)
{
	int j = 0;    // local loop indices
	short a = 0, b = 0;
	double y = 0.0;
	double pilot = 2.0 * M_PI / 100.0;
  
	// iterate over burst interval, number of channels
  for (j = 0; j < burstInterval; j++)
  {
    // calculate burst waveform
    if (j < 100)
    {
      // compute the selected waveform
      switch(theBurstType)
      {
        case simpleSine:     y = sin(pilot * j);  break;
        case harmonicSine:   y = sin(pilot * j) - sin(2.0 * pilot * j) / 2.0; break;
        case simpleCosine:   y = 1.0 - cos(pilot * j); break;
        case harmonicCosine: y = cos(2.0 * pilot * j) - cos(pilot * j); break;
        default: cout << "Unexpected burst type." << endl; return;
      }
      
      // normalize to +0 dB amplitude, convert to short word
      a = short(y * leftAmpl  * 32767.0);
      b = short(y * rightAmpl * 32767.0);
      if (!phase) {a = -a; b = -b;}
      if (typeAIFF == theFileType) {byteSwap(&a); byteSwap(&b);}
    }
    else {a = b = 0;}   // write silence after pilot
    
    // write left channel data, then right
    outfile.write((char *)&a, 2);
    outfile.write((char *)&b, 2);
  }
	return;
}

void writeAIFFfile(ostream &outFile, ostream &logFile)
{
  // set up AIFF header chunks
  int theSize = getSize();
  formatAIFFchunk theFormat(theSize);
  commonAIFFchunk theCommon(theSize, sampleRate);
  soundAIFFchunk  theSound(theSize);
  
  // write AIFF header to file
  outFile.write((char *)&theFormat, sizeof(theFormat));
  outFile.write((char *)&theCommon, sizeof(theCommon));
  outFile.write((char *)&theSound,  sizeof(theSound));
  
	// show column headings here
	logFile << "numCyc,duration,nomFreq,actFreq" << endl;
  
  // write initial delay to disk
  int i = 0;
  const short a = 0;
  for (i = 0; i < delay; i++)
  {
    // no need to swap bytes, it's just zero
    outFile.write((char *)&a, 2);
    outFile.write((char *)&a, 2);
  }
  
  if (pilot) {writePilotTone(outFile);}
  
  // iterate over tone bursts while writing to disk
	for(reset(); good(); next())
	{
		// check input file before reading
		if (!outFile)
		{
			logFile << "Failed while writing WAVE data." << endl;
			return;
		}
		showBurst(logFile);
    logFile << endl;
    
    // write AIFF data to file
    writeOutputData(outFile);
	}
}

void writeWAVEfile(ostream &outFile, ostream &logFile)
{
  // set up WAVE header chunks
  int theSize = getSize();
  riffWAVEchunk   theRiff(theSize);
  formatWAVEchunk theFormat(sampleRate);
  dataWAVEchunk   theData(theSize);
  
  // write WAVE header to file
  outFile.write((char *)&theRiff,   sizeof(theRiff));
  outFile.write((char *)&theFormat, sizeof(theFormat));
  outFile.write((char *)&theData,   sizeof(theData));
  
	// show column headings here
	logFile << "numCyc,duration,nomFreq,actFreq" << endl;
  
  // write initial delay to disk
  int i = 0;
  const short a = 0;
  for (i = 0; i < delay; i++)
  {
    outFile.write((char *)&a, 2);
    outFile.write((char *)&a, 2);
  }
  
  if (pilot) {writePilotTone(outFile);}
  
  // iterate over tone bursts while writing to disk
	for(reset(); good(); next())
	{
		// check input file before reading
		if (!outFile)
		{
			logFile << "Failed while writing WAVE data." << endl;
			return;
		}
    showBurst(logFile);
    logFile << endl;
    
    // write WAVE data to file
    writeOutputData(outFile);
	}
}

void writeOutputFile(istream &str)
{
  string outPath, logPath;
  str >> outPath >> logPath;
  ofstream logFile;
  if (logPath.length()) {logFile.open(logPath);}
  ostream &logStream = logPath.length()? logFile: cout;
  if (!outPath.length())
  {
    logStream << "Required output file name missing." << endl;
    return;
  }
  
  // don't overwrite existing file
  ifstream inFile(outPath);
  if (inFile)
  {
    inFile.close();
    logStream << "File already exists: " << outPath << endl;
    return;
  }
  
  // must be opened in binary mode
  ofstream outFile(outPath, ios::out | ios::binary);
  if (!outFile)
  {
    logStream << "Failed to open output file: " << outPath << endl;
    return;
  }
  logStream << "Writing output file," << outPath << endl;
  
  switch (theFileType)
  {
    case typeAIFF:
      writeAIFFfile(outFile, logStream);
      break;
      
    case typeWAVE:
      writeWAVEfile(outFile, logStream);
      break;
      
    default:
      logStream << "Unexpected output file type!" << endl;
      break;
  }
}


