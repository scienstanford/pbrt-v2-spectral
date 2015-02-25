#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <gsl/gsl_randist.h>

struct LensElement{
    float radius;
    float separation;
    float n;
    float aperture;
};

class RealisticDiffractionCamera : public Camera {
public:
   RealisticDiffractionCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
      float filmdiag,
      float curveradius,
	  Film *film,
      bool diffractionFlag,
      bool chromaticFlag,
      float xOffset,
      float yOffset,
      float pinholeExitXIn,
      float pinholeExitYIn,
      float pinholeExitZIn,
      float filmCenterXIn,
      float filmCenterYIn);

      
   ~RealisticDiffractionCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   float GenerateCameraRay(const CameraSample &sample, Ray *) const;
    void runLensFlare(const Scene * scene, const Renderer * renderer) const;
   float getFStop();
   float getFocalLength();
   float getSensorWidth();
private:
   float ShutterOpen;
   float ShutterClose;
   float filmDiag;
   float curveRadius;
   float filmDistance;
   float xApertureOffset;
   float yApertureOffset;
   float pinholeExitX;
   float pinholeExitY;
   float pinholeExitZ;
   float filmCenterX;
   float filmCenterY;
   gsl_rng * r;
   Film * film;
   vector<LensElement> lensEls;
   bool IntersectLensEl(const Ray &r, float *tHit, float radius, float dist,  Vector & normalVec) const;
   bool diffractionEnabled;
   bool chromaticAberrationEnabled;
   float fstop;
   float focalLength;

;
};

RealisticDiffractionCamera *CreateRealisticDiffractionCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
