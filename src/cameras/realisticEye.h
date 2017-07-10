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
                       bool IORforEyeFlagIn,
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
    
    // --------------------------------------------------------------
    // ---- Realistic index of refraction for each ocular medium ----
    // --------------------------------------------------------------
    
    // These are calculated using equations in Navarro 1985 and Atchison 2005.
    
    // R. Navarro, J. Santamaría, and J. Bescós, "Accommodation-dependent model of the human eye with aspherics," J. Opt. Soc. Am. A, vol. 2, no. 8, p.1273, Aug. 1985.
    // David A. Atchison and George Smith, "Chromatic dispersions of the ocular media of human eyes," J. Opt. Soc. Am. A 22, 29-37 (2005)
    
    int numEyeWaveSamples = 130;
    const float eyeWaveSamples[130] = {365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600,605,610,615,620,625,630,635,640,645,650,655,660,665,670,675,680,685,690,695,700,705,710,715,720,725,730,735,740,745,750,755,760,765,770,775,780,785,790,795,800,805,810,815,820,825,830,835,840,845,850,855,860,865,870,875,880,885,890,895,900,905,910,915,920,925,930,935,940,945,950,955,960,965,970,975,980,985,990,995,1000,1005,1010};
    
    const float corneaIORraw[130] = {1.398,1.396,1.395,1.394,1.393,1.392,1.391,1.390,1.389,1.388,1.388,1.387,1.386,1.386,1.385,1.385,1.384,1.384,1.383,1.383,1.382,1.382,1.381,1.381,1.381,1.380,1.380,1.380,1.380,1.379,1.379,1.379,1.379,1.378,1.378,1.378,1.378,1.377,1.377,1.377,1.377,1.377,1.376,1.376,1.376,1.376,1.376,1.376,1.375,1.375,1.375,1.375,1.375,1.375,1.375,1.374,1.374,1.374,1.374,1.374,1.374,1.374,1.374,1.373,1.373,1.373,1.373,1.373,1.373,1.373,1.373,1.373,1.372,1.372,1.372,1.372,1.372,1.372,1.372,1.372,1.372,1.372,1.371,1.371,1.371,1.371,1.371,1.371,1.371,1.371,1.371,1.371,1.370,1.370,1.370,1.370,1.370,1.370,1.370,1.370,1.370,1.370,1.369,1.369,1.369,1.369,1.369,1.369,1.369,1.369,1.369,1.369,1.369,1.368,1.368,1.368,1.368,1.368,1.368,1.368,1.368,1.368,1.368,1.367,1.367,1.367,1.367,1.367,1.367,1.367};
    
    const float aqueousIORraw[130] = {1.359,1.358,1.357,1.355,1.354,1.353,1.352,1.351,1.351,1.350,1.349,1.348,1.348,1.347,1.347,1.346,1.346,1.345,1.345,1.344,1.344,1.343,1.343,1.343,1.342,1.342,1.342,1.341,1.341,1.341,1.340,1.340,1.340,1.340,1.339,1.339,1.339,1.339,1.339,1.338,1.338,1.338,1.338,1.338,1.338,1.337,1.337,1.337,1.337,1.337,1.337,1.336,1.336,1.336,1.336,1.336,1.336,1.336,1.335,1.335,1.335,1.335,1.335,1.335,1.335,1.335,1.334,1.334,1.334,1.334,1.334,1.334,1.334,1.334,1.334,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.328,1.328,1.328,1.328,1.328,1.328};
    
    const float lensIORraw[130] = {1.449,1.447,1.446,1.444,1.443,1.441,1.440,1.439,1.438,1.437,1.436,1.435,1.434,1.433,1.432,1.431,1.431,1.430,1.429,1.429,1.428,1.428,1.427,1.427,1.426,1.426,1.425,1.425,1.425,1.424,1.424,1.424,1.423,1.423,1.423,1.422,1.422,1.422,1.422,1.421,1.421,1.421,1.421,1.420,1.420,1.420,1.420,1.419,1.419,1.419,1.419,1.419,1.419,1.418,1.418,1.418,1.418,1.418,1.418,1.417,1.417,1.417,1.417,1.417,1.417,1.417,1.416,1.416,1.416,1.416,1.416,1.416,1.416,1.415,1.415,1.415,1.415,1.415,1.415,1.415,1.415,1.415,1.414,1.414,1.414,1.414,1.414,1.414,1.414,1.414,1.414,1.413,1.413,1.413,1.413,1.413,1.413,1.413,1.413,1.413,1.413,1.412,1.412,1.412,1.412,1.412,1.412,1.412,1.412,1.412,1.412,1.411,1.411,1.411,1.411,1.411,1.411,1.411,1.411,1.411,1.411,1.411,1.410,1.410,1.410,1.410,1.410,1.410,1.410,1.410};
    
    const float vitreousIORraw[130] = {1.357,1.355,1.354,1.353,1.352,1.351,1.350,1.349,1.349,1.348,1.347,1.347,1.346,1.345,1.345,1.344,1.344,1.343,1.343,1.343,1.342,1.342,1.341,1.341,1.341,1.340,1.340,1.340,1.340,1.339,1.339,1.339,1.339,1.338,1.338,1.338,1.338,1.337,1.337,1.337,1.337,1.337,1.336,1.336,1.336,1.336,1.336,1.336,1.336,1.335,1.335,1.335,1.335,1.335,1.335,1.335,1.334,1.334,1.334,1.334,1.334,1.334,1.334,1.334,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.333,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.331,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.329,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.327,1.327};
    
    // We will load the above values into the following spectrums when the eye model is being set up in the constructor
    Spectrum corneaIOR;
    Spectrum aqueousIOR;
    Spectrum eyeLensIOR;
    Spectrum vitreousIOR;
    
};

RealisticEyeCamera *CreateRealisticEyeCamera(const ParamSet &params,
                                             const AnimatedTransform &cam2world, Film *film);

#endif
