//
//  water.h
//  pbrt
//
//  Created by Trisha Lian on 5/26/16.
//
//

#ifndef __pbrt__water__
#define __pbrt__water__

#include <stdio.h>
#include "volume.h"
#include <iostream>
#include <fstream> // For DEBUGGING

// WaterVolumeDensity Declarations
class WaterVolumeDensity : public VolumeRegion {
public:
    // WaterVolumeDensity Public Methods
    WaterVolumeDensity(float spc, float lpc, float cc, float domc, Spectrum plankton, const BBox &e, const Transform &v2w) {

        largePartConc = lpc;
        smallPartConc = spc;
        chlConc = cc;
        domConc = domc;
        planktonAbsorption = plankton;
        
        // -----------------------------------------------
        // --------- Print out water parameters  ---------
        // -----------------------------------------------
        std::cout << "Large particle concentration = " << largePartConc << std::endl;
        std::cout << "Small particle concentration = " << smallPartConc << std::endl;
        std::cout << "Chlorophyll concentration = " << chlConc << std::endl;
        std::cout << "Dissolved organic matter concentration = " << domConc << std::endl;
        
        // ---------------------------------------------
        // --------- Load/Calculate Constants ----------
        // ---------------------------------------------
        
        // A spectrum that holds pureWaterAbsorption
        pureWaterSpectrum = Spectrum::FromSampled(wavelengths, pw, waveSamples);
        
        // A spectrum that holds the wavelength numbers to be used in calculations below
        waveSpectrum = Spectrum::FromSampled(wavelengths, wavelengths, waveSamples);
        
        // TODO: We can precalculate these to save a little bit of time?
        lam0_over_lam = Pow((Spectrum(530.f)/waveSpectrum),4.3); //(see p.103 in LaW)
        VSFsmallConst =Pow((Spectrum(550.f)/waveSpectrum),1.7);
        VSFlargeConst =Pow((Spectrum(550.f)/waveSpectrum),0.3);
        
        // -----------------------------------
        // --------- Calculate sig_a ---------
        // -----------------------------------
        
        // Absorption from organic detritus
        if (chlConc != 0.f) {
            c_nap = 0.0124*(pow(chlConc,0.724))*Exp(-0.011*(waveSpectrum-440.f));
        }else{
            c_nap = 0.f;
        }
        
        // Absorption from dissolved organic matter
        if(domConc != 0.f){
            c_dom = domConc*Exp(-0.014*(waveSpectrum - 440.f));
        }else{
            c_dom = 0.f;
        }
        
        // Absorption from phytoplankton
        if (chlConc != 0.f) {
            // Get plankton absorption at 440 in order to scale it by the cholorophyll concentration
            float plankton440;
            planktonAbsorption.GetSpectrumAtWavelength(440, &plankton440);
            float scale = (0.0378*pow(chlConc,0.627))/plankton440;
            planktonAbsorption = planktonAbsorption*scale;
        }
        
        // Total absorption
        sig_a = pureWaterSpectrum + c_nap + c_dom + planktonAbsorption;
        
        // Absorption and scattering coefficients are in m^-1, but PBRT-spectral expects units as mm.
        // We adjust this difference here.
        sig_a = sig_a/(1000.f);
        
        // -----------------------------------
        // --------- Calculate sig_s ---------
        // -----------------------------------
        
        // To get sig_s we integrate the VSF over the 0 to pi
        // and multiply by 2*pi (see Sedlazech et al 2011)
        sig_s = 0.f;
        if(largePartConc != 0 && smallPartConc != 0){
            Spectrum currVSF;
            for (int angle = 0; angle < 181; angle++) {
                
                // Calculate VSF for this angle
                currVSF = sinf(angle*0.0174533f)*PhaseKopelevic(angle);
                
                // Sum up VSF's
                sig_s += currVSF;
            }
            sig_s *= 2*3.14159f*deg2rad; //2*pi * deg2rad(1)
            
            // Absorption and scattering coefficients are in m^-1, but PBRT-spectral expects units as mm.
            // We adjust this difference here.
            sig_s = sig_s/(1000.f);
        }
        
        
        WorldToVolume = Inverse(v2w);
        extent = e;
        
    }
    
    // PBRT (p.590): "Because the bound is maintained internally in the volume's object space, it must be transformed to world space for the WorldBound() method."
    BBox WorldBound() const {
        return Inverse(WorldToVolume)(extent);
    }
    
    // PBRT (p.590): "If the region's world-to-object-space transformation includes a rotation such that the volume isn't axis aligned in world space, it is possible to compute a tighter segment of the ray-volume overlap by transforming the ray to the volume's object space and doing the overlap test there." See Figure in book.
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    
    // If the renderer wants these values, we return it.
    Spectrum sigma_a(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? sig_a : 0.;
    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? sig_s : 0.;
    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? (sig_a + sig_s) : 0.;
    }
    
    // We don't have any emission for water.
    Spectrum Lve(const Point &p, const Vector &, float) const {
        return 0.;
    }
    
    // "Given a pair of directions, the VolumeRegion::p() returns the value of the phase function at the given point and the given time."
    // We'll want to use our own phase function (e.g. Kopelevic model)
    Spectrum p(const Point &p, const Vector &wi, const Vector &wo, float) const{
        
        if (!extent.Inside(WorldToVolume(p))) return 0.;
        if (largePartConc == 0 && smallPartConc == 0) return 0.;
        
        float costheta = Dot(wi, wo);
        float theta = acosf(costheta);
        
        if (isnan(theta)) {
            //Input vector is probably the same as the output vector.
            theta = 180;
        }
    
        return (PhaseKopelevic(theta)/sig_s);
    }
    
    
    // TODO: Move this into volume.cpp with all the other phase functions at some point.
    Spectrum PhaseKopelevic(float angle) const{

        float costheta = cosf(angle*0.0174533f);
        
        float coeff = (1.21e-4)*(1+0.835*costheta*costheta);
        Spectrum p = coeff*lam0_over_lam;
        
        p += smallPartConc*getVSFsmall(angle)*VSFsmallConst;
        p += largePartConc*getVSFlarge(angle)*VSFlargeConst;
        
        return p;
        
    }
    
    // This should return the optical thickness. \sigma_a and \sigma_s should be constant throughout the volume, so we can use Beer's Law (see p. 591 in the textbook)
    // Note: It seems like the volume renderer calls this method primarily.
    Spectrum tau(const Ray &ray, float, float) const {
        float t0, t1;
        if (!IntersectP(ray, &t0, &t1)) return 0.;
        return Distance(ray(t0), ray(t1)) * (sig_a + sig_s);
    }
    
    
    // Get the corresponding VSF value for a given angle with VSFsmall
    float getVSFsmall(float &angle) const{
        
        float value = 0;
        
        // Symmetric VSF
        if (angle < 0) {
            angle = -1.*angle;
        }
        
        //TODO: Binary search would be faster
        bool DEBUG_foundValue = false;
        float a0; float a1; float t;
        for (int i = 0; i < (numAngleSamples-1); i++) {
            a0 = angles[i];
            a1 = angles[i+1];
            
            if (angle == a0) {
                value = VSFsmall[i];
                DEBUG_foundValue = true;
                break;
            }
            
            if ((angle <= a1) && (angle > a0)) {
                DEBUG_foundValue = true;
                
                // Interpolate
                t = (angle - a0)/(a1-a0);
                value = Lerp(t, VSFsmall[i], VSFsmall[i+1]);
                break;
            }
        }
        
        // DEBUG
        if (!DEBUG_foundValue) {
            Error("Couldn't find a value for VSFsmall!!");
            std::cout << "Angle is " << angle << std::endl;
        }
        return value;
    }
    
    // Get the corresponding VSF value for a given angle with VSFsmall
    float getVSFlarge(float &angle) const{
        
        float value = 0;
        
        // Symmetric VSF
        if (angle < 0) {
            angle = -1.*angle;
        }
        
        //TODO: Binary search would be faster
        bool DEBUG_foundValue = false;
        float a0; float a1; float t;
        for (int i = 0; i < (numAngleSamples-1); i++) {
            a0 = angles[i];
            a1 = angles[i+1];
            
            if (angle == a0) {
                value = VSFlarge[i];
                DEBUG_foundValue = true;
                break;
            }
            
            if ((angle <= a1) && (angle > a0)) {
                DEBUG_foundValue = true;
                
                // Interpolate
                t = (angle - a0)/(a1-a0);
                value = Lerp(t, VSFlarge[i], VSFlarge[i+1]);
                break;
            }
        }
        
        // DEBUG
        if (!DEBUG_foundValue) {
            Error("Couldn't find a value for VSFsmall!!");
            std::cout << "Angle is " << angle << std::endl;
        }
        return value;
    }
    
private:
    // WaterVolumeDensity Private Data
                                
    // User Inputs
    float smallPartConc; // Small particle concentration
    float largePartConc; // Large particle concentration
    float chlConc; // Chlorophyll concentration
    float domConc; // Dissolved organic matter concentration
    Spectrum planktonAbsorption; //Plankton absorption curve
    BBox extent;
    
    // We will calculate these values using the above inputs
    Spectrum sig_a, sig_s, le;
    float g;
    Transform WorldToVolume;
    
    // Calculated data
    Spectrum c_nap;
    Spectrum c_dom;
    
    // Pre-calculated values
    Spectrum waveSpectrum;
    Spectrum pureWaterSpectrum;
    Spectrum lam0_over_lam; // For use in the volume scattering function of pure sea water (p.103 LaW)
    Spectrum VSFsmallConst; // For calculating the VSF for small/large particles
    Spectrum VSFlargeConst;
    
    float deg2rad = 0.0175f;
    
    // Wavelengths
    static const int waveSamples = 30;
    const float wavelengths[waveSamples] = {
        405,415,425,435,445,455,465,475,485,495,
        505,515,525,535,545,555,565,575,585,595,
        605,615,625,635,645,655,665,675,685,695};
    
    // Pure water absorption spectrum from literature
    // Wavelength from 405:10:695 (see wavelength)
    const float pw[waveSamples] = {
        0.0331,0.0286,0.0249,0.0214,
        0.0189,0.0179,0.0183,0.0190,
        0.0208,0.0235,0.0292,0.0371,
        0.0421,0.0467, 0.0528,0.0610,
        0.0708,0.0859,0.1216,0.1809,
        0.2431,0.2675,0.2854,0.2995,
        0.3146,0.3490,0.3931,0.4225,
        0.4708,0.5575};
    
    // VSF for small and large particles
    static const int numAngleSamples = 19;
    const float angles[numAngleSamples] = {
        0,0.5000,1.0000,1.5000,2.0000,4.0000,
        6.0000,10.0000,15.0000,30.0000,45.0000,
        60.0000,75.0000,90.0000,105.0000,120.0000,
        135.0000,150.0000, 180.0000
    };
    const float VSFsmall[numAngleSamples] = {
        5.300000,5.300000,5.200000,5.200000,
        5.100000,4.600000,3.900000,2.500000,
        1.300000,0.290000,0.098000,0.041000,
        0.020000,0.012000,0.008600,0.007400,
        0.007400,0.007500,0.008100};
    
    const float VSFlarge[numAngleSamples] = {
        140.000000,98.000000,46.000000,26.000000,
        15.000000,3.600000,1.100000,0.200000,
        0.050000,0.002800,0.000620,0.000380,
        0.000200,0.000063,0.000044,0.000029,
        0.000020,0.000020,0.000070};
    
};


WaterVolumeDensity *CreateWaterVolumeDensity(const Transform &volume2world,
                                                               const ParamSet &params);

#endif /* defined(__pbrt__water__) */
