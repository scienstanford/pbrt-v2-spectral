#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_EYE_H
#define PBRT_CAMERAS_EYE_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <gsl/gsl_randist.h>
#include "spectrum.h" // This is necessary to declare Spectrum class in this header file. However, it doesn't seem right to have to include this here...maybe we need to think about reorganizing this class?

// This is slightly different than our normal lens class. We have more variables to be able to describe the elements of the eye including biconic surfaces and unique ocular mediums.
struct LensElementEye{
    float radiusX;
    float radiusY;
    float thickness;
    float mediumIndex; // Corresponds to the mediumElement. Describes media directly "behind" (-z direction, toward the retina) the surface.
    float semiDiameter;
    float conicConstantX;
    float conicConstantY;
};

struct GRINelement{
//    float minWavelength, maxWavelength, refWavelength;
//    int K_max, L_max;
    float K1, K2, K3; // We precalculate these using the Sellmeier formula
    float L1, L2, L3;
    
};

// TEMP: Don't use this for now.
/*
struct MediumElement{
    float indexLabel; // Corresponds to the mediumIndex in LensElementEye
    bool GRINflag; // Whether or not this is a media directly "behind" (-z direction, toward the retina) of a GRIN lens
    vector<float> a; // For media described by Schott dispersion coefficients
    vector<float> k; // For media described by Sellmeier dispersion coefficients (e.g. the crystalline lens)
    vector<float> L; // For media described by Sellmeier dispersion coefficients (e.g. the crystalline lens)
};
*/

class RealisticEyeCamera : public Camera {
public:
    RealisticEyeCamera(const AnimatedTransform &cam2world,
                       Film *f,
                       float hither, float yon,
                       float sopen, float sclose,
                       string specfile,
                       float pupDiam,
                       float lensDecenterX,
                       float lensDecnterY,
                       float lensTiltX,
                       float lensTiltY,
                       float retinaDistance,
                       float retinaRadius,
                       float retinaSemiDiam,
                       bool chromaticFlag,
                       bool flipRad,
                       vector<Spectrum> iorSpectra,
                       int grinSurfaceIndex,
                       string grinFile);
    
    
    ~RealisticEyeCamera();
    float GenerateRay(const CameraSample &sample, Ray *) const;
    float GenerateCameraRay(const CameraSample &sample, Ray *) const;
    void runLensFlare(const Scene * scene, const Renderer * renderer) const;
    float getFStop();
    float getFocalLength();
    float getSensorWidth();
    
private:
    
    // General parameters for the camera class
    float ShutterOpen;
    float ShutterClose;
    Film * film;
    float filmDiag;
    
    // Lens information
    vector<LensElementEye> lensEls;
//    vector<MediumElement> lensMedia;
    float effectiveFocalLength;
    
    // Specific parameters for the human eye
    float pupilDiameter;
    float lensDecenterX;
    float lensDecenterY;
    float lensTiltX;
    float lensTiltY;
    float retinaDistance;
    float retinaRadius;
    float retinaSemiDiam;
    vector<Spectrum> iorSpectra;
    
    // GRIN lens
    bool grinLensFlag;
    int grinSurfaceIndex;
    vector<float> wavelengths;
    vector<float> radialSamples ;
    vector<Spectrum> grinIORsamples;

    // Flags for speed
    bool chromaticAberrationEnabled;
    bool IORforEyeEnabled;
    
    // Flags for conventions
    bool flipLensRadius;
    
    // Private methods for tracing through lens
    bool IntersectLensEl(const Ray &r, float *tHit, float radius, Vector dist, Vector & normalVec);
    bool IntersectLensElAspheric(const Ray &r, float *tHit, LensElementEye currElement, float zShift, Vector *n) const;
    void applySnellsLaw(float n1, float n2, float lensRadius, Vector &normalVec, Ray * ray ) const;
    float lookUpIOR(int mediumIndex, const Ray &ray) const;
    
    // Handy method to explicity solve for the z(x,y) at a given point (x,y), for the biconic SAG
    float BiconicZ(float x, float y, LensElementEye currElement) const;
    
};

RealisticEyeCamera *CreateRealisticEyeCamera(const ParamSet &params,
                                             const AnimatedTransform &cam2world, Film *film);

#endif
