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

WaterVolumeDensity::WaterVolumeDensity(Spectrum absorption, string vsfFile, const BBox &e, const Transform &v2w) {
    

    // --------- Read in VSF--------------
    ReadVSFFile(vsfFile, &VSF);

    
    // --------- Calculate sig_a ---------
    // Absorption and scattering coefficients are in m^-1, but PBRT-spectral expects units as mm.
    // We adjust this difference here.
    sig_a = absorption/(1000.f);
    

    // --------- Calculate sig_s ---------
    // To get sig_s we integrate the VSF over the 0 to pi
    // and multiply by 2*pi (see Sedlazech et al 2011, eq 7)
    sig_s = 0.f;
//    vector<Spectrum>::iterator angle_i;
    float currAngle = 0;
    for (int angle = 0; angle < 181; angle++) {
        // Sum up VSF's over 0 to 180 degrees
        sig_s += VSF[angle]*sinf(angle*deg2rad);
    }
    
    //This deg2rad is here is because of the way we are integrating above
    sig_s *= 2*3.14159f*deg2rad; //2*pi * deg2rad(1)
    
    // Absorption and scattering coefficients are in m^-1, but PBRT-spectral expects units as mm.
    // We adjust this difference here.
    sig_s = sig_s/(1000.f);
    

    // ------- PBRT Transforms -----------
    WorldToVolume = Inverse(v2w);
    extent = e;
    
}


void WaterVolumeDensity::ReadVSFFile(const string &filename, vector<Spectrum> *vsf){
    
    // Check to see if the VSF filename is valid
    string vsfFile = AbsolutePath(ResolveFilename(filename));
    vector<float> vsfSpectrums;
    if (!ReadFloatFile(vsfFile.c_str(), &vsfSpectrums)) {
        Error("Unable to read VSF file!");
        return;
    }
    
    // Find the number of samples in each spectrum. This should be saved as the first value in the file.
    int numWlsSamples = (int)vsfSpectrums[0];
    if ((vsfSpectrums.size()-1) % numWlsSamples != 0){
        Error("The number of wavelength samples specified in the VSF spectrum is incorrect!");
        return;
    }
    
    // Find the wavelengths. This should be the first row in the file.
    vector<float>::const_iterator first = vsfSpectrums.begin() + 1;
    vector<float>::const_iterator last = vsfSpectrums.begin() + numWlsSamples + 1;
    vector<float> wavelengths(first, last);
    
    
    // Read in each row of float as a spectrum and place in VSF vector
    Spectrum currSpectrum;
    vector<float>::iterator it;
    for (it = vsfSpectrums.begin() + numWlsSamples + 1; it < vsfSpectrums.end(); it += numWlsSamples){
        
        vector<float> currRow(it, it+numWlsSamples);
        currSpectrum = Spectrum::FromSampled(&wavelengths[0], &currRow[0], numWlsSamples); // Pass pointer to first element in vector, since FromSampled takes an array
        vsf->push_back(currSpectrum);
        
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
    
    // Find the corresponding VSF spectrum for this angle
    
    float thetaDegrees = theta*rad2deg;
    Spectrum vsf_thisAngle;
    bool foundSpectrum = false; // For debugging
    
    // We linearly interpolate between spectra
    for(int angle = 0; angle < (181-1); angle++){
        if(thetaDegrees > angle & thetaDegrees <= (angle+1)){
            float t = thetaDegrees - angle;
            vsf_thisAngle = Lerp(t, VSF[angle], VSF[angle+1]);
            foundSpectrum = true;
            break;
        }
    }
    
    // Error check
    if (!foundSpectrum) {
        std::cout << "thetaDegrees = " << thetaDegrees << std::endl;
        Error("Couldn't find phase function for this angle!");
    }
    
    Spectrum test = (vsf_thisAngle/1000.f)/sig_s;
    
    return ((vsf_thisAngle/1000.f)/sig_s);
}

WaterVolumeDensity *CreateWaterVolumeDensity(const Transform &volume2world,
                                             const ParamSet &params) {
    
    Spectrum absorption = params.FindOneSpectrum("absorptionCurveFile", 0.f);
    string vsfFile = params.FindOneString("vsfFile","");
    
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
    
    return new WaterVolumeDensity(absorption, vsfFile, BBox(p0, p1),
                                  volume2world);
}