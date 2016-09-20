#include "SC_PlugIn.h"
#include <tgmath.h>

using namespace std;


// BASIC ADMINISTRATION

static InterfaceTable *ft;

struct BLOsc1section : public Unit
{
float partialOmega; // 2pi/sr : phase increment per sample while frequency = 1 
float phaseCurrent;

// For section 1
float loHarmonicsIntCurrent1;
float loHarmonicsFracCurrent1;
float hiHarmonicsIntCurrent1;
float hiHarmonicsFracCurrent1;
float spreadCurrent1;


// To check if each argument is "ar" or "kr"
int freqRate; // ar:0  kr:1
int loHarmonicsRate1; // ar:0  kr:1
int hiHarmonicsRate1; // ar:0  kr:1
int slopeRate1; // ar:0  kr:1
int evenOddRatioRate1; // ar:0  kr:1

};

extern "C"
{
  void BLOsc1section_next(BLOsc1section *unit, int inNumSamples);
  void BLOsc1section_Ctor(BLOsc1section* unit);
};

//////////////////////////////////////////////////////////////////

// CONSTRUCTOR

void BLOsc1section_Ctor(BLOsc1section* unit)
{

  SETCALC(BLOsc1section_next);

  unit->partialOmega = twopi_f / SAMPLERATE;
  unit->phaseCurrent = 0.0f;
  unit->loHarmonicsIntCurrent1 = 0.0f;
  unit->loHarmonicsFracCurrent1 = 0.0f;
  unit->hiHarmonicsIntCurrent1 = 0.0f;
  unit->hiHarmonicsFracCurrent1 = 0.0f;
  unit->spreadCurrent1 = 0.0f;

  if (INRATE(0) == calc_FullRate) unit->freqRate = 0; else unit->freqRate = 1;

  if (INRATE(1) == calc_FullRate) unit->loHarmonicsRate1 = 0; else unit->loHarmonicsRate1 = 1;
  if (INRATE(2) == calc_FullRate) unit->hiHarmonicsRate1 = 0; else unit->hiHarmonicsRate1 = 1;
  if (INRATE(3) == calc_FullRate) unit->slopeRate1 = 0; else unit->slopeRate1 = 1;
  if (INRATE(4) == calc_FullRate) unit->evenOddRatioRate1 = 0; else unit->evenOddRatioRate1 = 1;
 

  BLOsc1section_next(unit, 1);
}

//////////////////////////////////////////////////////////////////

// Function to calculate amp factor

float ampFactor_calc (float slope, float evenOddFactor, float loHarmonicsIntCurrent, float hiHarmonicsIntCurrent, float loEvenHarmonics, float hiEvenHarmonics, float numEvenHarmonics, float fundamentalAdjust, float extraHarmonicsAdjust)
{
  float result;
  result = slope == 1.0f? 
(hiHarmonicsIntCurrent - loHarmonicsIntCurrent + 1.0f) - (evenOddFactor * numEvenHarmonics) + fundamentalAdjust + extraHarmonicsAdjust
:
(pow(slope,loHarmonicsIntCurrent) - pow(slope,hiHarmonicsIntCurrent+1.0f))
/
(1.0f - slope)
- 
(evenOddFactor * 
  (pow(slope, loEvenHarmonics)
    -
    pow(slope, hiEvenHarmonics + 2.0f)
    )
  /
  (1.0f - pow(slope, 2.0f))) 
+
fundamentalAdjust * pow(slope,loHarmonicsIntCurrent-1.0f)
+
extraHarmonicsAdjust * pow(slope, hiHarmonicsIntCurrent+1.0f)
;
return result;
}
//ampFactor will be used to normalize the output amplitude. To avoid the denominator of this calculation to be 0 when slope = 1, the different formula is used when slope == 1.


// Function to calculte the first and the last fractional harmonics
complex<float> fractionalFundamental_calc (float phaseCurrent, float slope, float loHarmonicsIntCurrent, float fundamentalAdjust)
{
  complex<float> result;
  result = 
  fundamentalAdjust * 
  pow(slope,loHarmonicsIntCurrent-1.0f) * 
  exp(complex<float>(0.0f, (loHarmonicsIntCurrent-1.0f)*phaseCurrent));
  return result;
}

complex<float> fractionalExtraHarmonics_calc (float phaseCurrent, float slope, float hiHarmonicsIntCurrent, float extraHarmonicsAdjust)
{
  complex<float> result;
  result = 
  extraHarmonicsAdjust * 
  pow(slope,hiHarmonicsIntCurrent+1.0f) * 
  exp(complex<float>(0.0f, (hiHarmonicsIntCurrent+1.0f)*phaseCurrent));
  return result;
}


// Function to calculate complex signal
complex<float> complexSignal_calc (complex<float> baseOsc, complex<float> fractionalFundamental, complex<float> fractionalExtraHarmonics, float phaseCurrent, float slope, float evenOddFactor, float spreadCurrent, int spreadCompensation, float loHarmonicsIntCurrent, float hiHarmonicsIntCurrent, float loEvenHarmonics, float hiEvenHarmonics, float ampFactor)
{
  complex<float> result;
  if (slope == 1.0f && phaseCurrent == 0.0f) {
    result = complex<float>(0.0f, 1.0f);
  }
  else {
    result = 
    (
      pow(slope * baseOsc, loHarmonicsIntCurrent)
      -
      pow(slope * baseOsc, hiHarmonicsIntCurrent + 1.0f)
      )
    /
    (1.0f - slope * baseOsc)
    - 
    evenOddFactor * 
    (
      pow(slope * baseOsc, loEvenHarmonics)
      -
      pow(slope * baseOsc, hiEvenHarmonics + 2.0f)
      )
    /
    (1.0f - pow(slope * baseOsc, 2.0f))
    + fractionalFundamental
    + fractionalExtraHarmonics
    ;

    result = pow(result/ampFactor, spreadCurrent);
  };

/* FREQUENCY SHIFTER */
  if (spreadCompensation != 0) {
   result = polar(abs(result), arg(result) - (phaseCurrent * loHarmonicsIntCurrent * (spreadCurrent - 1.0f)));
 }
  // if spreadCompensation is not 0, frequency shifts down according to given spread value (in other words, it always keeps the same lowest harmonic number)
  // if spread is 0, frequency shift does NOT happen (Therefore it keeps the original result of power function)

  return result;
}


// Function to calculate one block as real signal (imag = sine phase)

float oneBlock_calc (complex<float> baseOsc, float phaseCurrent, float slope, float evenOddRatio, float spreadCurrent, int spreadCompensation, float loHarmonicsIntCurrent, float hiHarmonicsIntCurrent, float loHarmonicsFracCurrent, float hiHarmonicsFracCurrent)
{
  float loEvenHarmonics = static_cast<int>(loHarmonicsIntCurrent)%2 == 0? loHarmonicsIntCurrent : loHarmonicsIntCurrent + 1.0f; // The lowest even harmonic index
float hiEvenHarmonics = static_cast<int>(hiHarmonicsIntCurrent)%2 == 0? hiHarmonicsIntCurrent : hiHarmonicsIntCurrent - 1.0f; // The highest even harmonic index
float numEvenHarmonics = (hiEvenHarmonics - loEvenHarmonics) / 2.0f + 1.0f; //The total number of even harmonics

float evenOddFactor = 1.0f - evenOddRatio;

  float fundamentalAdjust = static_cast<int>(loHarmonicsIntCurrent)%2 == 0? loHarmonicsFracCurrent : evenOddRatio * loHarmonicsFracCurrent;
float extraHarmonicsAdjust = static_cast<int>(hiHarmonicsIntCurrent)%2 == 0? hiHarmonicsFracCurrent : evenOddRatio * hiHarmonicsFracCurrent;

  float ampFactor = ampFactor_calc(slope, evenOddFactor, loHarmonicsIntCurrent, hiHarmonicsIntCurrent, loEvenHarmonics, hiEvenHarmonics, numEvenHarmonics, fundamentalAdjust, extraHarmonicsAdjust);

complex<float> fractionalFundamental = fractionalFundamental_calc(phaseCurrent, slope, loHarmonicsIntCurrent, fundamentalAdjust);
complex<float> fractionalExtraHarmonics = fractionalExtraHarmonics_calc(phaseCurrent, slope, hiHarmonicsIntCurrent, extraHarmonicsAdjust);

complex<float> signalC = 
complexSignal_calc(baseOsc, fractionalFundamental, fractionalExtraHarmonics, phaseCurrent, slope, evenOddFactor, spreadCurrent, spreadCompensation, loHarmonicsIntCurrent, hiHarmonicsIntCurrent, loEvenHarmonics, hiEvenHarmonics, ampFactor); 
// signal as complex number

float signalR = signalC.imag(); // signal as real number
return signalR;
}

//////////////////////////////////////////////////////////////////

// UGEN CALCULATION

void BLOsc1section_next(BLOsc1section *unit, int inNumSamples)
{

float *out = ZOUT(0);

float *freqIn = ZIN(0);
float freq = ZIN0(0);


/* Arguments for section 1 */
float *loHarmonicsIn1 = ZIN(1);
float loHarmonics1 = ZIN0(1);

float *hiHarmonicsIn1 = ZIN(2);
float hiHarmonics1 = ZIN0(2);

float *slopeIn1 = ZIN(3);
float slope1 = ZIN0(3);

float *evenOddRatioIn1 = ZIN(4);
float evenOddRatio1 = ZIN0(4);

int spreadIn1 = ZIN0(5);
float spreadCurrent1 = unit->spreadCurrent1 == 0.0f? static_cast<float>(spreadIn1) : unit->spreadCurrent1;
// initial value = ZIN0(5), value after running = unit->spreadCurrent


int spreadCompensation = ZIN0(6);

float phaseCurrent = unit->phaseCurrent;

/* Section 1 */
float loHarmonicsIntCurrent1 = unit->loHarmonicsIntCurrent1 == 0.0f? ceil(loHarmonics1) : unit->loHarmonicsIntCurrent1;
float hiHarmonicsIntCurrent1 = unit->hiHarmonicsIntCurrent1 == 0.0f? floor(hiHarmonics1) : unit->hiHarmonicsIntCurrent1;
// unit->... == 0? means set up initial value when starting the UGen

float loHarmonicsFracCurrent1 = unit->loHarmonicsFracCurrent1 == 0.0f? loHarmonicsIntCurrent1 - loHarmonics1 : unit->loHarmonicsFracCurrent1;
float hiHarmonicsFracCurrent1 = unit->hiHarmonicsFracCurrent1 == 0.0f? hiHarmonics1 - hiHarmonicsIntCurrent1 : unit->hiHarmonicsFracCurrent1;
// unit->... == 0? means set up initial value when starting the UGen


LOOP(inNumSamples,
  
complex<float> baseOsc = exp(complex<float>(0.0f, phaseCurrent));


/* Section 1 */
float loHarmonicsInt1;
float loHarmonicsFrac1;
float hiHarmonicsInt1;
float hiHarmonicsFrac1;
float signalR1;


if (unit->loHarmonicsRate1 == 0) loHarmonics1 = ZXP(loHarmonicsIn1); 
if (unit->hiHarmonicsRate1 == 0) hiHarmonics1 = ZXP(hiHarmonicsIn1); 
if (unit->slopeRate1 == 0) slope1 = ZXP(slopeIn1); 
if (unit->evenOddRatioRate1 == 0) evenOddRatio1 = ZXP(evenOddRatioIn1); 

loHarmonicsInt1 = ceil(loHarmonics1);
loHarmonicsFrac1 = loHarmonicsInt1 - loHarmonics1;
if (loHarmonicsIntCurrent1 == loHarmonicsInt1) loHarmonicsFracCurrent1 = loHarmonicsFrac1;

hiHarmonicsInt1 = floor(hiHarmonics1);
hiHarmonicsFrac1 = hiHarmonics1 - hiHarmonicsInt1;
if (hiHarmonicsIntCurrent1 == hiHarmonicsInt1) hiHarmonicsFracCurrent1 = hiHarmonicsFrac1;

signalR1 = oneBlock_calc(baseOsc, phaseCurrent, slope1, evenOddRatio1, spreadCurrent1, spreadCompensation, loHarmonicsIntCurrent1, hiHarmonicsIntCurrent1, loHarmonicsFracCurrent1, hiHarmonicsFracCurrent1);
;

float omega;
if (unit->freqRate == 0) freq = ZXP(freqIn);
omega = freq * unit->partialOmega;
// phase increment per sample at given frequency

  phaseCurrent += omega;
  while (phaseCurrent >= twopi_f){
    phaseCurrent -= twopi_f;

   /* Section 1 */
    if (loHarmonicsInt1 != loHarmonicsIntCurrent1) loHarmonicsIntCurrent1 = loHarmonicsInt1;
    // Change loHarmonicsInt value when phase crosses 2pi 

    if (hiHarmonicsInt1 != hiHarmonicsIntCurrent1) hiHarmonicsIntCurrent1 = hiHarmonicsInt1;
    // Change hiHarmonicsInt value when phase crosses 2pi 

    if (static_cast<float>(spreadIn1) != spreadCurrent1) spreadCurrent1 = static_cast<float>(spreadIn1);
    // Change spread value when phase crosses 2pi

}
  ZXP(out) = signalR1;
  )

unit->phaseCurrent = phaseCurrent;

unit->loHarmonicsIntCurrent1 = loHarmonicsIntCurrent1;
unit->hiHarmonicsIntCurrent1 = hiHarmonicsIntCurrent1;
unit->loHarmonicsFracCurrent1 = loHarmonicsFracCurrent1;
unit->hiHarmonicsFracCurrent1 = hiHarmonicsFracCurrent1;
unit->spreadCurrent1 = spreadCurrent1;

}

////////////////////////////////////////////////////////////////////

// LOAD FUNCTION

PluginLoad(BLOsc1section)
{
  ft = inTable;

  DefineSimpleUnit(BLOsc1section);
}

////////////////////////////////////////////////////////////////////

