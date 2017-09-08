//
//  water.cpp
//  pbrt
//
//  Created by Trisha Lian on 5/26/16.
//
//

#include "stdafx.h"
#include "water.h"
#include "paramset.h"




// WaterVolumeDensity Method Definitions

WaterVolumeDensity::WaterVolumeDensity(Spectrum absorption, Spectrum scattering, string phaseFile, const BBox &e, const Transform &v2w) {
    

    // --------- Read in phase file --------------
    ReadPhaseFile(phaseFile, &phaseFunction);

    
    sig_a = absorption/1000.f; // mm^-1
    sig_s = scattering/1000.f; // mm^-1

    // ------- PBRT Transforms -----------
    WorldToVolume = Inverse(v2w);
    extent = e;
    
}


void WaterVolumeDensity::ReadPhaseFile(const string &filename, vector<Spectrum> *phase){
    
    // Check to see if the VSF filename is valid
    string phaseFile = AbsolutePath(ResolveFilename(filename));
    vector<float> phaseSpectrums;
    if (!ReadFloatFile(phaseFile.c_str(), &phaseSpectrums)) {
        Error("Unable to read phase function file!");
        return;
    }
    
    // Find the number of samples in each spectrum. This should be saved as the first value in the file.
    int numWlsSamples = (int)phaseSpectrums[0];
    if ((phaseSpectrums.size()-1) % numWlsSamples != 0){
        Error("The number of wavelength samples specified in the phase function spectrum is incorrect!");
        return;
    }
    
    // Find the wavelengths. This should be the first row in the file.
    vector<float>::const_iterator first = phaseSpectrums.begin() + 1;
    vector<float>::const_iterator last = phaseSpectrums.begin() + numWlsSamples + 1;
    vector<float> wavelengths(first, last);
    
    
    // Read in each row of float as a spectrum and place in phase function vector
    Spectrum currSpectrum;
    vector<float>::iterator it;
    for (it = phaseSpectrums.begin() + numWlsSamples + 1; it < phaseSpectrums.end(); it += numWlsSamples){
        
        vector<float> currRow(it, it+numWlsSamples);
        currSpectrum = Spectrum::FromSampled(&wavelengths[0], &currRow[0], numWlsSamples); // Pass pointer to first element in vector, since FromSampled takes an array
        phase->push_back(currSpectrum);
        
        currRow.clear();
    }
    
    return;
}

Spectrum WaterVolumeDensity::p(const Point &p, const Vector &wi, const Vector &wo, float) const{
    
    if (!extent.Inside(WorldToVolume(p))) return 0.;
    if (sig_s.IsBlack()) return 0.;
    
    float costheta = Dot(wi, wo);
    float theta = acosf(costheta);
    
    // TODO: Does this ever happen?
    if (isnan(theta)) {
        //theta is probably 180
        theta = 2*3.14159f;
    }
    
    // Find the corresponding phase spectrum for this angle
    
    float thetaDegrees = theta*rad2deg;
    
    Spectrum phase_thisAngle;
    bool foundSpectrum = false; // For debugging
    
    // We linearly interpolate between spectra
    for(int angle = 0; angle < (181-1); angle++){
        if(thetaDegrees >= angle & thetaDegrees < (angle+1)){
            float t = thetaDegrees - angle;
            phase_thisAngle = Lerp(t, phaseFunction[angle], phaseFunction[angle+1]);
            foundSpectrum = true;
            break;
        }
    }
    
    // Error check
    if (!foundSpectrum) {
        std::cout << "thetaDegrees = " << thetaDegrees << std::endl;
        Error("Couldn't find phase function for this angle!");
    }
    
    return phase_thisAngle;
}

WaterVolumeDensity *CreateWaterVolumeDensity(const Transform &volume2world,
                                             const ParamSet &params) {
    
    Spectrum absorption = params.FindOneSpectrum("absorptionCurveFile", 0.f); // m^-1
    Spectrum scattering = params.FindOneSpectrum("scatteringCurveFile", 0.f); // m^-1
    string phaseFile = params.FindOneString("phaseFunctionFile",""); // sr^-1
    
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
    
    return new WaterVolumeDensity(absorption, scattering, phaseFile, BBox(p0, p1),
                                  volume2world);
}