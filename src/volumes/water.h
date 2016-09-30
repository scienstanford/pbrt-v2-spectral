//
//  water.h
//  pbrt
//
//  Created by Trisha Lian on 5/26/16.
//
//

#ifndef __pbrt__water__
#define __pbrt__water__

#include "floatfile.h"// For ReadFloatFile
#include "fileutil.h" // For ResolveFilename
#include <stdio.h>
#include "volume.h"
#include <iostream>
#include <fstream> // For DEBUGGING

// WaterVolumeDensity Declarations
class WaterVolumeDensity : public VolumeRegion {
public:
    // WaterVolumeDensity Public Methods
    WaterVolumeDensity(Spectrum absorption, string vsfFile, const BBox &e, const Transform &v2w);
    
    // PBRT (p.590): "Because the bound is maintained internally in the volume's object space, it must be transformed to world space for the WorldBound() method."
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent);}
    
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
    Spectrum Lve(const Point &p, const Vector &, float) const {return 0.;}
    
    // "Given a pair of directions, the VolumeRegion::p() returns the value of the phase function at the given point and the given time."
    // We'll want to use our own phase function (e.g. Kopelevic model)
    Spectrum p(const Point &p, const Vector &wi, const Vector &wo, float) const;
    
    // This should return the optical thickness. \sigma_a and \sigma_s should be constant throughout the volume, so we can use Beer's Law (see p. 591 in the textbook)
    // Note: It seems like the volume renderer calls this method primarily.
    Spectrum tau(const Ray &ray, float, float) const {
        float t0, t1;
        if (!IntersectP(ray, &t0, &t1)) return 0.;
        return Distance(ray(t0), ray(t1)) * (sig_a + sig_s);
    }
    
    /*
     ------- READ VSF FILE ------
     The VSF is a 2D matrix with the rows being the angles (0 to 180 degrees) and the columns being wavelength samples.
     We read this in from a text file outputted by MATLAB code.
     This first value should be the number of wavelengths in each spectrum.
     This first row (after the first value) should be the wavelength samples. We need this to create the spectrum from sampled values.
     E.g.
     
     30
     400 410 420 ... 700
     0_1  0_2  0_3  ... 0_30
     1_1  1_2  1_3  ... 1_30
     ...
     ...
     180_1 180_2 180_3 ... 180_30
     
     WARNING: VSF rows MUST be in degrees, at 1 degree intervals (i.e. 0,1,2,3...,178,179,180)
     In other words there MUST be 181+1+1 = 183 rows in the text file
     */
    void ReadVSFFile(const string &filename, vector<Spectrum> *vsf);
    

private:
    // WaterVolumeDensity Private Data
    BBox extent;
    
    // Read this in from user
    Spectrum sig_a;
    vector<Spectrum> VSF;
    
    // Calculated
    Spectrum sig_s;
    Transform WorldToVolume;
    
    float deg2rad = 0.0174533f;
    float rad2deg = 57.2958f;
    
};


WaterVolumeDensity *CreateWaterVolumeDensity(const Transform &volume2world,
                                                               const ParamSet &params);

#endif /* defined(__pbrt__water__) */
