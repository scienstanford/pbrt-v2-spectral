
#include "stdafx.h"
#include "cameras/realisticDiffraction.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/spectralImage.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"
#include "floatfile.h"
#include "camera.h"
#include "scene.h"
#include "lights/point.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <time.h>
#include "../vdb.h"
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

using namespace std;

RealisticDiffractionCamera *CreateRealisticDiffractionCamera(const ParamSet &params,
                                                             const AnimatedTransform &cam2world, Film *film) {
    // -----------------------------------
    // ----- Common Camera Parameters ----
    // -----------------------------------
    
    float hither = params.FindOneFloat("hither", -1); // Do we need this?
    float yon = params.FindOneFloat("yon", -1); // Do we need this?
    float shutteropen = params.FindOneFloat("shutteropen", -1);
    float shutterclose = params.FindOneFloat("shutterclose", -1);

    // ---------------------------------------------------
    // ----- Realistic Diffraction Camera Parameters -----
    // ---------------------------------------------------
    
    // Check for lens file
    string specfile = params.FindOneString("specfile", "");
    if (specfile == "") {
        Severe( "No lens spec file supplied!\n" );
    }
    float filmdistance = params.FindOneFloat("filmdistance", 70.0);
    float apdiameter = params.FindOneFloat("aperture_diameter", 1.0);
    float filmdiag = params.FindOneFloat("filmdiag", 35.0);

    float curveradius = params.FindOneFloat("curveRadius", 0);

    float xOffset = params.FindOneFloat("x_aperture_offset", 0);
    float yOffset = params.FindOneFloat("y_aperture_offset", 0);
    
    bool diffractFlag = params.FindOneBool("diffractionEnabled", 1.0);
    bool chromaticFlag = params.FindOneBool("chromaticAberrationEnabled", 0.0);
    
    // -----------------------------------------
    // ----- Light Field Camera Parameters -----
    // -----------------------------------------
    
    float pinholeExitApX = params.FindOneFloat("pinhole_exit_x", -1);
    float pinholeExitApY = params.FindOneFloat("pinhole_exit_y", -1);
    float pinholeExitApZ = params.FindOneFloat("pinhole_exit_z", -1);
    float filmcenterX = params.FindOneFloat("film_center_x", 0);
    float filmcenterY = params.FindOneFloat("film_center_y", 0);
    
    int numPinholesW = (int) params.FindOneFloat("num_pinholes_w", -1);
    int numPinholesH = (int) params.FindOneFloat("num_pinholes_h", -1);
    bool microlensFlag = params.FindOneFloat("microlens_enabled", 0);
    
    // -----------------------------------------
    // ----- Gullstrand Eye Parameters ---------
    // -----------------------------------------
    // Whether to use realistic index of refractions for tracing through the Gullstrand eye.
    bool IORforEyeFlag = params.FindOneBool("IORforEyeEnabled", 0.0);
    
	   
    // TL: We don't use hither, yon, or the shutter parameters in our code, so no need to assert.
    /*	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
    	      shutterclose != -1 && filmdistance!= -1);
     */

    
    return new RealisticDiffractionCamera(cam2world, hither, yon,
                                          shutteropen, shutterclose, filmdistance, apdiameter,
                                          specfile, filmdiag, curveradius, film, diffractFlag, chromaticFlag, xOffset, yOffset, pinholeExitApX, pinholeExitApY, pinholeExitApZ, filmcenterX, filmcenterY, numPinholesW, numPinholesH, microlensFlag,IORforEyeFlag);
}




RealisticDiffractionCamera::RealisticDiffractionCamera(const AnimatedTransform &cam2world,
                                                       float hither, float yon,
                                                       float sopen, float sclose,
                                                       float filmdistance, float aperturediameter,
                                                       const string &specfile,
                                                       float filmdiag,
                                                       float curveradius,
                                                       Film *f,
                                                       bool diffractFlag,
                                                       bool chromaticFlag,
                                                       float xOffset,
                                                       float yOffset,
                                                       float pinholeExitXIn,
                                                       float pinholeExitYIn,
                                                       float pinholeExitZIn,
                                                       float filmCenterXIn,
                                                       float filmCenterYIn,
                                                       int numPinholesWIn,
                                                       int numPinholesHIn,
                                                       bool microlensFlagIn,
                                                       bool IORforEyeFlagIn)
: Camera(cam2world, sopen, sclose, f), ShutterOpen(sopen), ShutterClose(sclose),film(f)
{
    

    // Store passed in values into private variables
    filmDiag = filmdiag;
    filmDistance = filmdistance;
    string fn = AbsolutePath(ResolveFilename(specfile)); //parse the specfile
    diffractionEnabled = diffractFlag;
    chromaticAberrationEnabled = chromaticFlag;
    xApertureOffset = xOffset;
    yApertureOffset = yOffset;
    curveRadius = curveradius;
    pinholeExitX = pinholeExitXIn;
    pinholeExitY = pinholeExitYIn;
    pinholeExitZ = pinholeExitZIn;
    filmCenterX = filmCenterXIn;
    filmCenterY = filmCenterYIn;
    numPinholesW = numPinholesWIn;
    numPinholesH = numPinholesHIn;
    microlensFlag = microlensFlagIn;
    IORforEyeFlag = IORforEyeFlagIn;
    
    // Display Light Field parameters
    /*
    std::cout <<"xApertureOffset: " << xApertureOffset << "\n";
    std::cout <<"yApertureOffset: " << yApertureOffset << "\n";
    std::cout <<"DiffractionEnabled: " << diffractionEnabled << "\n";
    std::cout <<"chromaticAberrationEnabled: " << chromaticAberrationEnabled << "\n";
    std::cout <<"pinholeExitApX: " << pinholeExitX << "\n";
    std::cout <<"pinholeExitApY: " << pinholeExitY << "\n";
    std::cout <<"pinholeExitApX: " << pinholeExitX << "\n";
    std::cout <<"filmCenterX: " << filmCenterX << "\n";
    std::cout <<"filmCenterY: " << filmCenterY << "\n";
    std::cout <<"numPinholesW: " << numPinholesW << "\n";
    std::cout <<"numPinholesH: " << numPinholesH << "\n";
    std::cout <<"microlensFlag: " << microlensFlag << "\n";
    */
    
    // -------------------------
    // --- Read in lens file ---
    // -------------------------

    vector<float> vals;
    float lastAperture = 0; // This holdes the diameter of the last lens. What is it for? (TL)
    
    // Check to see if there is valid input in the lens file.
    if (!ReadFloatFile(fn.c_str(), &vals)) {
        Warning("Unable to read lens file!");
        return;
    }

    // Check to see if the number of floats (not including focal length) is a multiple of 4
    if ((vals.size()-1) % 4 != 0)
    {
        Warning("Wrong number of float values in lens file!  Check file format!  Did you forget to specify the focal length?");
        return;
    }
    
    focalLength = vals[0];   // Read the focal length
    fstop = focalLength/aperturediameter;
    
    for (int i = 1; i < vals.size(); i+=4)
    {
        LensElement currentLensEl;
        currentLensEl.radius = vals[i];
        currentLensEl.separation = vals[i+1];
        currentLensEl.n = vals[i+2];
        currentLensEl.aperture = vals[i+3];
        if (currentLensEl.radius == 0)
            currentLensEl.aperture = aperturediameter;
        lensEls.push_back(currentLensEl);
        lastAperture = currentLensEl.aperture;
    }
    
    // -------------------------
    // --- Set Up Eye Model ----
    // -------------------------
    
    // If we are ray tracing through the eye and want chromatic aberration, we set up the IOR curves to be used, i.e. n(lambda)
    // Load values into a spectrum class. These values are interpolated to the wave samples specified in "spectrum.h"
    if(IORforEyeFlag){
        corneaIOR = Spectrum::FromSampled(eyeWaveSamples, corneaIORraw, numEyeWaveSamples);
        aqueousIOR = Spectrum::FromSampled(eyeWaveSamples, aqueousIORraw, numEyeWaveSamples);
        eyeLensIOR = Spectrum::FromSampled(eyeWaveSamples, lensIORraw, numEyeWaveSamples);
        vitreousIOR = Spectrum::FromSampled(eyeWaveSamples, vitreousIORraw, numEyeWaveSamples);
        
        //DEBUG
        std::cout << "Eye IOR curves have been set." << std::endl;
    }
    
    // Display Camera parameters
    
     /*
     std::cout << "focalLength :" << focalLength << "\n";
     std::cout << "apertureDiameter: " << aperturediameter << "\n";
     std::cout << "fstop: f/" << fstop << "\n";
     std:: cout << "curveRadius: " << curveRadius << "\n";
     
     // Note: This is read out the way it's traced through.
     for (int i = lensEls.size()-1; i >=0 ; i--)
     {
     std::cout << "radius: " << lensEls[i].radius << " ";
     std::cout << "separation: " << lensEls[i].separation << " ";
     std::cout << "n: " << lensEls[i].n << " ";
     std:: cout << "aperture: " << lensEls[i].aperture << "\n";
     }
     */
     
    // --- Random Number Generator ---
    // I think this was originally intended for diffraction.
    //code borrowed from http://www.gnu.org/software/gsl/manual/gsl-ref.html#Random-Number-Generation  18.13
    const gsl_rng_type * T;
    
    int i, n = 10;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    // ---------------------------------
    // --- Set Up Light Field Camera ---
    // ---------------------------------
    
    // TL: This section has not been cleaned up yet.
    
    // AL: Figure out pinhole locations film distance.
    // we solved for this algebraically/using similar triangles
	// we are assuming a square sensor for now
    
	   if (numPinholesW > 0 && numPinholesH > 0)
       {
           
           float sPixPitch = getSensorWidth()/((float)numPinholesW);
           float pinholeArrayDistance = sPixPitch * filmDistance/(lastAperture + sPixPitch);
           std::cout << "superpixel pitch: " << sPixPitch << "\n";
           std::cout << "pinholeArrayDistance: " << pinholeArrayDistance << "\n";
           
           //loop through all the pinhole locations and figure out the coordinates
           //store the pinhole locations in a matrix
           
           float pinholePosition = -filmDistance + pinholeArrayDistance;
           pinholeArray = new Point*[numPinholesW];
           
           for (int i = 0; i < numPinholesW; i++)
           {
               pinholeArray[i] = new Point[numPinholesH];
               for (int j = 0; j < numPinholesH; j++)
               {
                   Point centerPos;
                   centerPos.x = -(i - numPinholesW/2.0 + .5) * sPixPitch;  //potential rounding problems here?
                   centerPos.y = (j- numPinholesH/2.0 + .5) * sPixPitch;
                   centerPos.z = -filmDistance;  //pbrt handedness is film is negative, scene positive
                   
                   Point lensCenter(0,0,0); //this is assuming the lens is centered here!
                   
                   Vector chiefRayDir;
                   chiefRayDir = lensCenter -  centerPos;
                   chiefRayDir = Normalize(chiefRayDir);
                   
                   float tHit = pinholePosition/chiefRayDir.z;
                   
                   Point currentPinhole;
                   currentPinhole.x = tHit * chiefRayDir.x;   //ray originates at the lens center, and is traced back
                   currentPinhole.y = tHit * chiefRayDir.y;
                   
                   //currentPinhole.x = tHit * chiefRayDir.x + centerPos.x;
                   //currentPinhole.y = tHit * chiefRayDir.y + centerPos.y;
                   
                   currentPinhole.z = pinholePosition;
                   pinholeArray[i][j] = currentPinhole;
                   
                   std::cout << "x: " << pinholeArray[i][j].x << " ";
                   std::cout << "y: " << pinholeArray[i][j].y << " ";
                   std::cout << "z: " << pinholeArray[i][j].z << "\t";
                   
                   //check for microlenses...
                   //calculate microlens specs
                   
               }
               std::cout <<"\n";
           }
       }
	   
    
    
}


RealisticDiffractionCamera::~RealisticDiffractionCamera()
{
}


float RealisticDiffractionCamera::getFStop()
{
    return fstop;
}

float RealisticDiffractionCamera::getFocalLength()
{
    return focalLength;
}

// TL: I don't think Andy got this to work.
void RealisticDiffractionCamera::runLensFlare(const Scene * scene, const Renderer * renderer) const
{
    /*    std::cout << "in runLensFlare()!\n";
     
     for (int i = 0; i < scene->lights.size(); i++)
     {
     std::cout << "isPointLight?: " << scene->lights[i]->isPointLight() << "\n";
     if (!scene->lights[i]->isPointLight())
     continue;
     
     //if it is a point light, make lens flare!  shoot rays into the lens!
     Point pointLight = ((const PointLight*)(scene->lights[i]))->getLightPos();
     std::cout << "pointLight position: " << pointLight.x << ", " << pointLight.y << ", " << pointLight.z << "\n";
     
     //now we want to shoot rays uniformly from the light source to the front of the lens
     
     }*/
    
}

void RealisticDiffractionCamera::applySnellsLaw(float n1, float n2, float lensRadius, Vector &normalVec, Ray * ray ) const
{
    
    // n1 = IOR for medium closer to sensor
    // n2 = IOR for medium closer to scene
    
    // scene....n2 <----> Lens Element <----> n1 ....sensor
    
    
    //Andy: add chromatic aberration effect (changing index of refraction) here - basic for now
    // TODO: There probably doesn't need to be a flag for chromatic aberration and then another flag for the eye IOR. What's the best way to combine these? Anyway for now, let's just make sure the IOR works.
    
    // In the case that the user made a mistake and specified the IORforEyeFlag but not the spectralrenderer, the rays don't have a wavelength. They seem to be assigned random values that are outside the visible wavelength. Here I'm checking to make sure they're within the sampledLambda ranges in spectrum.h. This is not a very safe fix, so let's think about how to do this in the future (TODO: Make this less hacky?).
    
    if (IORforEyeFlag && (ray->wavelength != 0)){
        
        // Look up the accurate n1 and n2 values according to the ray wavelength.
        // How do we determine which ocular medium we're in right now? Right now I'm just going to check the n value loaded from the lens file, but this might not be consistent long term.
        // TODO: What's the best way to organize this? Maybe we should just make subclasses of realisticDiffraction somehow.
        if (abs(n1-1.336) < 0.001) {
            vitreousIOR.GetValueAtWavelength(ray->wavelength,&n1);
            eyeLensIOR.GetValueAtWavelength(ray->wavelength, &n2);
        }else if(abs(n1-1.42) < 0.001){
            eyeLensIOR.GetValueAtWavelength(ray->wavelength, &n1);
            aqueousIOR.GetValueAtWavelength(ray->wavelength, &n2);
        }else if(abs(n1-1.3374) < 0.001){
            aqueousIOR.GetValueAtWavelength(ray->wavelength, &n1);
            corneaIOR.GetValueAtWavelength(ray->wavelength, &n2);
        }else if(abs(n1-1.3771) < 0.001){
            corneaIOR.GetValueAtWavelength(ray->wavelength, &n1);
            n2 = 1;
        }
    
    }else if(chromaticAberrationEnabled)
    {
        if (n1 != 1)
            n1 = (ray->wavelength - 550) * -.04/(300)  +  n1;              //should be .04
        if (n2 != 1)
            n2 = (ray->wavelength - 550) * -.04/(300)  +  n2;
    }

    // ------ DEBUG
    // Verify that the IOR's have changed:
    /*
    std::cout << "Wavelength: " << ray->wavelength << std::endl;
    std::cout << "n1 = " << n1 << std::endl;
    std::cout << "n2 = " << n2 << std::endl;
     */
    // ------
    
    Vector s1 = ray->d;
    if (lensRadius >0)
        normalVec = -normalVec;

    float radicand = 1 - (n1/n2) * (n1/n2) * Dot(Cross(normalVec, s1), Cross(normalVec, s1));
    if (radicand < 0){
        ray->d = Vector(0, 0,0);
        // Tlian: I put this return back in so the assert doesn't catch.
        return;   //reflection, no refraction - might want to change for lens flare!
    }
    
    Vector s2 = n1/n2 * (Cross(normalVec, Cross(-1 * normalVec, s1))) - normalVec * sqrt(radicand);
    ray->d = Normalize(s2);  //reassign the direction to the new direction
}

bool RealisticDiffractionCamera::IntersectLensEl(const Ray &r, float *tHit, float radius, Vector dist, Vector & normalVec) const{
    
    // TL: This is the same code that is in sphere.cpp
    // We have to be careful about where the zero axis is in all of this.
    // TODO: Double check Andy's calculations. For now I'm assuming its correct.
    
    float phi;
    Point phit;
    
    // Transform _Ray_ to object space
    Transform shiftZ =  Translate(dist);
    //Transform shiftZ =  Translate(Vector(0,0,radius - dist));
    
    Ray ray = r;
    (shiftZ)(r, &ray);
    //(*WorldToObject)(r, &ray);
    
    
    if (radius < 0)
        radius = -radius;
    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
    ray.o.z*ray.o.z - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;

    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }
    
    // Compute sphere hit position and $\phi$
    phit = ray(thit);
    if (phit.x == 0.f && phit.y == 0.f) phit.x = 1e-5f * radius;
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;
    
    
    // Update _tHit_ for quadric intersection
    *tHit = thit;
    
    // Compute _rayEpsilon_ for quadric intersection
    //*rayEpsilon = 5e-4f * *tHit;
    normalVec = Normalize(Vector(ray.d.x * thit + ray.o.x, ray.d.y*thit + ray.o.y, ray.d.z*thit + ray.o.z));
    
    return true;
}

float RealisticDiffractionCamera::getSensorWidth()
{
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    return width;
}

float RealisticDiffractionCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
    
    // TL:                     z=0         z= -filmDistance
    //                  |  |    ||           |
    //                  |  |    ||           |
    //  Scene <---------|--|----||---------> |
    //            +z    |  |    ||    -z     |
    //                  |  |    ||           |
    //                Lens Elements        Sensor
    //
    
    
    
    
    // use sample->imageX and sample->imageY to get raster-space coordinates
    // of the sample point on the film.
    // use sample->lensU and sample->lensV to get a sample position on the lens
    
    // GenerateRay() should return the weight of the generated ray
    
    Point startingPoint;
    //std::cout << "imageX: " << sample.imageX << "imageY: " << sample.imageY << "\n";
    // startingPoint.x = -(sample.imageX -2 - film->xResolution/2.f + .5)/(film->xResolution/2.f);
    // startingPoint.y = (sample.imageY-2 - film->yResolution/2.f + .5)/(film->yResolution/2.f);
    
    startingPoint.x = -((sample.imageX) - film->xResolution/2.f - .25)/(film->xResolution/2.f);  //potential rounding/alignment problem here
    startingPoint.y = ((sample.imageY) - film->yResolution/2.f - .25)/(film->yResolution/2.f);
    startingPoint.z = -filmDistance;
    
    //cout << "startingPoint.z: " << startingPoint.z;
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    //float width = filmDiag /sqrt((1.f + aspectRatio * aspectRatio));
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    float height = width/aspectRatio;
    
    startingPoint.x = startingPoint.x * width/2.f + filmCenterX;
    startingPoint.y = startingPoint.y * height/2.f + filmCenterY;
    
    //curved Sensor stuff - sensor offset with sperical sensors still needs to be evaluated
    if (curveRadius != 0)
    {
        //convert to angle
        float startTheta = startingPoint.x/curveRadius;
        float startPhi = startingPoint.y/curveRadius;
        
        //convert to x,y,z, on sphere
        startingPoint.x = curveRadius * cos(startPhi) * sin(startTheta);   //sign convention needs to be checked
        startingPoint.z = curveRadius * cos(startPhi) * cos(startTheta);
        startingPoint.y = curveRadius * sin(startPhi);
        float sphereCenter =  (-filmDistance - curveRadius) ;
        startingPoint.z = sphereCenter + startingPoint.z;  //check the sign convention here...
        
    }
    
    float lensU, lensV;
    float prevN = 1;
    
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    
    //use geometry to solve for z intercept of disc on sphere
    float firstAperture = lensEls[lensEls.size()-1].aperture/2;
    float firstRadius = lensEls[lensEls.size()-1].radius;
    
    //special case for when aperture is 1st element
    float zIntercept = 0;
    if (firstRadius ==0)
        zIntercept = 0;
    else
        zIntercept = (-firstRadius - sqrt(firstRadius * firstRadius - firstAperture * firstAperture));
    
    //experiment
    //zIntercept = 0;
    
    //std::cout << " firstRadius: " << firstRadius;
    //std::cout << " firstAperture: " << firstAperture;
    //std::cout << " zIntercept: " << zIntercept;
    
    //TODO: adjust first aperture depending on if we are using depth-map mode or NOT
    
    float lensUNoScale = lensU;
    float lensVNoScale = lensV;
    float microlensAperture = 1;  //the aperture for the microlens.  Should be the same size as superpixel pitch.
    float superpixelPitch = 0;    //pitch in mm of the superpixels (a super pixel constitutes the collection of pixels behind 1 pinhole/microlens)
    int xPinhole = 0;
    int yPinhole = 0;
    
    lensU *= firstAperture;
    lensV *= firstAperture;
    
    Point pointOnLens = Point(lensU, lensV, zIntercept);   //can we even assume that lens is a flat disk, maybe this has problems later?
    
    //special case when we have an offset pinhole
    if(pinholeExitX != -1 && pinholeExitY != -1 && pinholeExitZ != -1)
    {
        pointOnLens = Point(pinholeExitX, pinholeExitY, pinholeExitZ);
    }
    
    //special case when we have a pinhole array
    else if(numPinholesW >0 && numPinholesH >0)
    {
        //determine which pinhole to trace to
        
        int numPixelsPerPinholeW = film->xResolution/numPinholesW;  //figure out number of pixels per superpixel in both dirs
        int numPixelsPerPinholeH = film->yResolution/numPinholesH;
        
        xPinhole = (int)(((sample.imageX - .25 )/numPixelsPerPinholeW));    //figure out which pinhole we are under right now
        yPinhole = (int)(((sample.imageY - .25)/numPixelsPerPinholeH));    //same with y
        
        if (xPinhole > (numPinholesW - 1))
        {
            //for some weird reason... imageX starts at about 1 -> 721 ... why?  is this a pbrt bug?
            std::cout << "warning: xPinhole is out of range (too big)! " << xPinhole << " " << sample.imageX << "\n";
            xPinhole = numPinholesW - 1;
        }
        else if (xPinhole < 0)
        {
            std::cout << "warning: xPinhole is out of range (too small)!" << xPinhole << " " << sample.imageX << "\n";
            xPinhole = 0;
        }
        if (yPinhole > (numPinholesH - 1))
        {
            std::cout << "warning: yPinhole is out of range (too big)!\n";
            yPinhole = numPinholesH - 1;
        }
        else if (yPinhole < 0)
        {
            std::cout << "warning: yPinhole is out of range (too small)!\n";
            yPinhole = 0;
        }
        //std::cout <<"numPixsPerPinholeW" << numPixelsPerPinholeW << "\n";
        //std::cout << "numPixsPerPinholeH" << numPixelsPerPinholeH << "\n";
        //std::cout <<"xPinhole" << xPinhole << "\t";
        //	std::cout << "yPinhole" << yPinhole << "\n";
        //std::cout << "film->xResolution: " << film->xResolution << "\t" << film->yResolution << "\n";
        //std::cout <<"xsample: " << sample.imageX << " ysample" << sample.imageY << "\n";
        
        Point pinholeLocation = pinholeArray[xPinhole][yPinhole];
        
        //std::cout << "pinholeLocation: " << pinholeLocation.x << " " << pinholeLocation.y << " " << pinholeLocation.z << 			"\n";
        //std::cout << "startingPoint: " << startingPoint.x << " " << startingPoint.y << " " << startingPoint.z << "\n";
        
        
        //check for microlens boolean
        //change pointOnLens to be a rectangle at the entrance aperture of microlens
        if(microlensFlag)
        {
            
            //shoot rays through at the superpixel aperture
            superpixelPitch = width/numPinholesW;
            //float superpixelPitchH = height/numPinholesH;
            
            //this isn't right... we need the pinhole pitch instead... but will do for now... TODO: fix!!...but it's close
            pointOnLens = Point(lensUNoScale * superpixelPitch/2.f + pinholeLocation.x, lensVNoScale * superpixelPitch/2.f+ pinholeLocation.y, pinholeLocation.z);
            
            //pointOnLens = pinholeLocation;
        }
        else
        {
            pointOnLens = pinholeLocation;
        }
    }
    
    //wanted focal length of microlenses right now .23739
    //superpixel pitch: 0.106066
    //assume these specs for now.  do on the fly calculations later!
    
    
    
    //SPECIAL CASE WHEN WE HAVE A MICROLENS ARRAY
    //take sampleX and sampleY and scale them so they make up the square of microlens
    //then, trace rays through the microlenses (should we use a function of some sort?)
    
    float tempWavelength = ray->wavelength;
    *ray = Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
    ray->o = startingPoint;    //initialize ray origin
    ray->d = Normalize(pointOnLens - ray->o);
    ray->wavelength = tempWavelength;  //so that wavelength information is retained.
    
    // vdb_color(.7, 0, 0);
    // vdb_point(pointOnLens.x, pointOnLens.y, pointOnLens.z);
    
    
    //trace through microlenses! (if necessary)
    //the microlens has 2 spherical lens elements
    if (microlensFlag && numPinholesW >0 && numPinholesH >0)
    {
        //we assume an R of 0.3143 for a focal length of f = .23739.  calculate for real in further mod
        
        //float lensRadius = .3558;//0.3762; //.3191;  //calculated from matlab script s_s3dLFMicroLensCalculate.m
        float currentAperture = superpixelPitch;
        
        //these are example values for now.  we will compute them in a subsequent version
        float microlensThickness = .01;   //mm
        float microlensFilmDistance = filmDistance+pointOnLens.z;//.232; //.237; //mm
        float microlensFocalLength = microlensFilmDistance + microlensThickness/2; //.237;//mm   //this is to be used for for ideal microlenses
        
        float microlensN = 1.67; //index of refraction of microlens
        float oneOverR = (-2 * (microlensN-1) + sqrt(4 * pow(microlensN-1,2) +
                                                     4 * pow(microlensN-1, 2) * microlensThickness/(microlensN * microlensFocalLength)))/(2 * pow((microlensN-1), 2) * microlensThickness/microlensN);
        float lensRadius = 1/oneOverR; //solved for lens radius given focal length and mirolensN
        float lensDistance = -filmDistance + microlensFilmDistance;
        
        /*
         std::cout <<"microlensFilmDistance: " << microlensFilmDistance <<"\n";
         std::cout <<"microlensFocalLength: " << microlensFocalLength << "\n";
         std::cout << "oneOverR: " << oneOverR << "\n";
         std::cout << "lensRadius: " << lensRadius << "\n";*/
        
        Vector normalVec(0,0,1);
        Point intersectPoint(0,0,0);
        
        ray->o = startingPoint;
        
        
        /*
         //use ideal microlenses for now to see if the problem is with lens aberrations, the code below, or even something else.  this is a sanity check
         
         //trace chief ray
         Ray * centerRay = new Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
         Point centerOfLens = Point(pinholeArray[xPinhole][yPinhole].x,pinholeArray[xPinhole][yPinhole].y,pinholeArray[xPinhole][yPinhole].z);
         centerRay->o = startingPoint;
         centerRay->d = Normalize(centerOfLens - startingPoint);
         
         centerRay->o = centerOfLens;  //put origin at origin (center of lens)
         
         //calculate chief ray intersection
         
         //use thin lens equation to compute new point of focus
         //this is the distance to a point that is in focus.  all rays from a single point must converge at this distance.
         float focalDistance = 0;
         float denominator = 1.f/microlensFocalLength - 1.f/microlensFilmDistance;
         if (abs(microlensFocalLength - microlensFilmDistance) < .00000000001f)   //protect against division by 0
         focalDistance = 100000000000 * microlensFocalLength;
         else
         focalDistance = 1.f/denominator;
         
         focalDistance -= filmDistance;
         
         //compute point of focus in 3 space (world coordinates).
         //first, compute t intersection with in focus plane - then position
         float ft = focalDistance / centerRay->d.z;
         //Pfocus is the object point in world coordinates
         Point Pfocus = (*centerRay)(ft);  //PBRT notation assigns point on a ray, if given the desired time
         
         //intersectPoint is a point in world coordinates corresponding to where the ray hits the pupil plane
         //Point intersectPoint;
         //intersectPoint.x = lensU;
         //intersectPoint.y = lensV;
         //intersectPoint.z = 0;
         
         ray->o = pointOnLens;
         ray->d = Normalize(Pfocus - pointOnLens);
         
         startingPoint = pointOnLens;
         */
        
        
        //this loop is extremely similar to the one below.  However, to prevent me breaking the loop below, I made a separate one to decouple them for better debugging.  They can be combined in the future if you wish
        //i = 0 case: first microlens element (closest to sensor;  i = 1 case: 2nd microlens element
        
        
        for (int i = 0; i < 2; i++)
        {
            float dummy = 0;
            float *tHit;
            tHit = &dummy;
            bool intersected = false;
            Vector normalVec(0,0,1);
            Point intersectPoint(0,0,0);
            Point microlensCenter(pinholeArray[xPinhole][yPinhole].x, pinholeArray[xPinhole][yPinhole].y,lensDistance);
            
            lensRadius = -lensRadius;  //to account for the 2 opposite polarity radii
            
            
            ray->o = startingPoint;
            
            //we should change the lens position to the correct one (with the correct offset)
            //we might consider making the next loop a while loop instead and incorporating all of this in there so that we don't need to remake everything...
            intersected = IntersectLensEl(*ray, tHit, lensRadius, Vector(-microlensCenter.x, -microlensCenter.y, lensRadius - lensDistance), normalVec);   //
            
            if (intersected)
            {
                intersectPoint.x = *tHit * ray->d.x + ray->o.x;
                intersectPoint.y = *tHit * ray->d.y + ray->o.y;
                intersectPoint.z = *tHit * ray->d.z + ray->o.z;
                
                //aperture test
                if ((intersectPoint.x - microlensCenter.x) * (intersectPoint.x - microlensCenter.x) +
                    (intersectPoint.y - microlensCenter.y) * (intersectPoint.y - microlensCenter.y) >= (currentAperture * currentAperture) / (2 * 2))
                    return 0.f;
                
                //call Snell's Law Function
                float n1 = 1;
                float n2 = 1;
                
                if (i == 0)
                {
                    n2 = microlensN;
                    lensDistance += microlensThickness;
                }
                else
                    n1 = microlensN;
                
                applySnellsLaw( n1,  n2,  lensRadius, normalVec, ray );
                
                //update starting point to current point
                startingPoint.x = intersectPoint.x;
                startingPoint.y = intersectPoint.y;
                startingPoint.z = intersectPoint.z;
            }
            
            
            // --------------------add effect of diffraction----------------
            
            if (diffractionEnabled)
            {
                
                double initx = 0;
                double inity = 0;
                double * x = &initx;
                double * y = &inity;
                //double currentAperture = lensRadius * 2;
                
                //calculate min distance direction
                
                //calculate radius
                //Point intersectPoint(lensU, lensV, 0.f);
                double radius = sqrt((intersectPoint.x - microlensCenter.x) * (intersectPoint.x - microlensCenter.x) + (intersectPoint.y - microlensCenter.y) * (intersectPoint.y - microlensCenter.y) );
                
                //cout << "radius: " << radius << "\n";
                //std::cout << "rayWavelength: " << ray->wavelength << "\n";
                
                //calculate direction
                Vector direction(intersectPoint.x, intersectPoint.y, 0);
                
                Vector orthoDirection(-intersectPoint.y, intersectPoint.x, 0);
                // double direction = atan (intersectPoint.y/intersectPoint.x);
                
                double a = currentAperture/2 - radius;
                double b = sqrt((currentAperture/2 * currentAperture/2) - radius * radius);
                
                a = a;
                b = b;
                double pi = 3.14159265359;
                //double lambda = .000000550;  //550 nanometers for now
                double lambda = ray->wavelength * 1e-9;
                double sigma_x = atan(1/(sqrt(2) * a *.001 * 2 * pi/lambda));
                double sigma_y = atan(1/(sqrt(2) *b * .001 * 2 * pi/lambda));
                
                //gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)
                gsl_ran_bivariate_gaussian (r, sigma_x, sigma_y, 0, x, y);    //experiment for now
                
                //add r.v. in directions of direction and orthoDirection
                
                //calculate component of these vectors based on 2 random degrees
                direction = Normalize(direction);
                orthoDirection = Normalize(orthoDirection);
                
                float noiseA = (float)(*x);
                float noiseB = (float)(*y);
                
                //project the original ray onto the new bases
                double projA = (ray->d.x * direction.x + ray->d.y * direction.y)/sqrt(direction.x * direction.x + direction.y * direction.y);
                double projB = (ray->d.x * orthoDirection.x + ray->d.y * orthoDirection.y)/sqrt(orthoDirection.x * orthoDirection.x + orthoDirection.y * orthoDirection.y);
                double projC = ray->d.z;
                
                double rA = sqrt(projA * projA + projC * projC);
                double rB = sqrt(projB * projB + projC * projC);
                double thetaA = acos(projA/rA);
                double thetaB = acos(projB/rB);
                
                //add uncertainty
                thetaA = thetaA + noiseA;   //removed for now
                thetaB = thetaB + noiseB;
                
                //convert angles back into cartesian coordinates, but in a,b space
                double newProjA = cos(thetaA) * rA;
                
                ray->d.z = sin(thetaA) * rA;
                
                projC = ray->d.z;
                rB = sqrt(projB * projB + projC * projC);
                
                //rB = sqrt(ray->d.y * ray->d.y + ray->d.z * ray->d.z); // THIS LINE IS WRONG!!!!!
                
                //need to recalculate thetaB after you modify z (this is a new addition, also questionable)
                thetaB = acos(projB/rB);
                
                double newProjB = cos(thetaB) * rB;
                ray->d.z = sin(thetaB) * rB;
                
                //convert from ab space back to x,y space
                ray->d.x = direction.x * newProjA + orthoDirection.x * newProjB;
                ray->d.y = direction.y * newProjA + orthoDirection.y * newProjB;
                
                ray->d = Normalize(ray->d);
            }
            //------------end diffraction-----------------
            
            
        }
        
    }
    
    // --------------------------------------------------------
    // --- Trace through the lens elements of the main lens ---
    // --------------------------------------------------------
    
    
    float lensDistance = 0; // How far we are from the "0" on the z-axis, i.e. the position of the lens closest to the sensor
    
    for (int i = lensEls.size()-1; i>=0 ; i--)
    {
        
        float lensRadius = lensEls[i].radius;
        lensDistance += lensEls[i].separation;
        float currentAperture = lensEls[i].aperture;
        
        ray->o = startingPoint;

        float dummy = 0;
        float *tHit;
        tHit = &dummy;
        
        bool intersected = false;
        Vector normalVec(0,0,1);
        Point intersectPoint(0,0,0);
        
        if (lensRadius ==0)
        {
            // ---------------------
            // --- APERTURE CASE ---
            // ---------------------
            
            float tAperture = 0;
            if (i == lensEls.size()-1)
                tAperture = filmDistance/ray->d.z;   //special case for when aperture is the 1st element
            else
                tAperture = (lensDistance - ray->o.z)/(ray->d.z);
            
            // Point where the ray intersects the aperture plane
            Point apertureIntersect = Point(ray->o.x + ray->d.x * tAperture,ray->o.y + ray->d.y * tAperture,ray->o.z + ray->d.z * tAperture);
            
            // For the offset aperture case
            if ((apertureIntersect.x - xApertureOffset) * (apertureIntersect.x - xApertureOffset) +  (apertureIntersect.y - yApertureOffset) * (apertureIntersect.y - yApertureOffset) > currentAperture * currentAperture * .25)
                return 0.f;
            
            normalVec = Vector(0,0,1);
            
            intersectPoint = apertureIntersect;
            startingPoint = intersectPoint;
            
        }
        else
        {
            // ----------------------------
            // --- REGULAR ELEMENT CASE ---
            // ----------------------------
            
            //---- Find intersection point ----
            
            // "Vector(0, 0, lensRadius - lensDistance)" is where the center of the sphere is with respect to the "0" point on the z-axis, (i.e. the lens element closest to the sensor)
            intersected = IntersectLensEl(*ray, tHit, lensRadius, Vector(0, 0, lensRadius - lensDistance), normalVec);
            
            if (intersected)
            {
                intersectPoint = Point(*tHit * ray->d.x + ray->o.x,*tHit * ray->d.y + ray->o.y,*tHit * ray->d.z + ray->o.z);
                
                // Check to see if the ray is outside of the lens size (i.e. aperture)
                if (intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y >= currentAperture * currentAperture / 4.f)
                    return 0.f;
                
                
                // TL: Andy had some way of visualizing the rays...
                // vdb visualization
                
                /*
                 //float currentT = 0;
                 
                 while (currentT < *tHit)
                 {
                 Point currentPoint(0,0,0);
                 currentPoint.x = currentT * ray->d.x + ray->o.x;
                 currentPoint.y = currentT * ray->d.y  + ray->o.y;
                 currentPoint.z = currentT * ray->d.z  + ray->o.z;
                 
                 // if (currentT == 0)
                 //     vdb_color(.7, 0, 0);
                 // else
                 //     vdb_color(0, 0, .7);
                 // vdb_point(currentPoint.x, currentPoint.y, currentPoint.z);
                 currentT += .5;
                 
                 vdb_color(0, .7, 0);
                 vdb_point(intersectPoint.x, intersectPoint.y, intersectPoint.z);
                 
                 }*/
                
                
                // ---- Apply Snell's Law ----
                
                float n1 = lensEls[i].n;
                float n2;
                
                // Andy:
                // default case for n2 is 1 because that is air.  n2/n1 direction is like this....
                // n2----->|lens element| ---->n1 (n2 is previous surface n.. n1 is new surface)
                
                // TL:
                // scene ..... n2 <-----> |lens element| <----> n1  ..... sensor
                
                
                // If we're not at the lens closest to the scene, we use the next IOR
                // If we're at the end, we assume it's air in the scene
                if (i-1 >= 0){
                    
                    n2 = lensEls[i-1].n;
                    
                    // Trisha: If we're entering the aperture (n2 == 0) we skip the n2 == 0 and apply the next medium.
                    // (We put this in the current if statement so we can handle the unique case when the aperture is the first element.)
                    
                    //                 |    /
                    //                 |   |
                    //                 |  |
                    // <-----            |    <------
                    //  n2, [i-2]      |  |        n1, [i]
                    //                 |   |
                    //                 |    \
                    //
                    //           aperture, [i-1]
                    
                    if(n2 == 0)
                        n2 = lensEls[i-2].n;
                    
                }
                else{
                    n2 = 1;
                }

                
                applySnellsLaw( n1,  n2,  lensRadius, normalVec, ray );
                
                // --- Update ray starting point ---
                
                startingPoint = intersectPoint;
                //std::cout << "(" << ray->o.x << "," << ray->o.y << "," << ray->o.z << ")" << std::endl;
                
            }
            else
            {
                
                // TL: More of Andy's visualization code
                /*
                 //did not intersect -> draw red rays
                 float currentT = 0;
                 while (currentT < 10)
                 {
                 Point currentPoint(0,0,0);
                 currentPoint.x = currentT * ray->d.x + ray->o.x;
                 currentPoint.y = currentT * ray->d.y  + ray->o.y;
                 currentPoint.z = currentT * ray->d.z  + ray->o.z;
                 
                 //if (currentT == 0)
                 //    vdb_color(.7, 0, 0);
                 //else
                 //    vdb_color(1, 0, 0);
                 //vdb_point(currentPoint.x, currentPoint.y, currentPoint.z);
                 currentT += .5;
                 }*/
                
                return 0.f;
            }
        }
        
        
        // --------------------add effect of diffraction----------------
        
        if (diffractionEnabled)
        {
            double initx = 0;
            double inity = 0;
            double * x = &initx;
            double * y = &inity;
            //double currentAperture = lensRadius * 2;
            
            //calculate min distance direction
            
            //calculate radius
            //Point intersectPoint(lensU, lensV, 0.f);
            double radius = sqrt(intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y );
            
            //cout << "radius: " << radius << "\n";
            //std::cout << "rayWavelength: " << ray->wavelength << "\n";
            
            //calculate direction
            Vector direction(intersectPoint.x, intersectPoint.y, 0);
            
            Vector orthoDirection(-intersectPoint.y, intersectPoint.x, 0);
            // double direction = atan (intersectPoint.y/intersectPoint.x);
            
            double a = currentAperture/2 - radius;
            double b = sqrt((currentAperture/2 * currentAperture/2) - radius * radius);
            
            a = a;
            b = b;
            double pi = 3.14159265359;
            //double lambda = .000000550;  //550 nanometers for now
            double lambda = ray->wavelength * 1e-9;
            double sigma_x = atan(1/(sqrt(2) * a *.001 * 2 * pi/lambda));
            double sigma_y = atan(1/(sqrt(2) *b * .001 * 2 * pi/lambda));
            
            //gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double * x, double * y)
            gsl_ran_bivariate_gaussian (r, sigma_x, sigma_y, 0, x, y);    //experiment for now
            
            //add r.v. in directions of direction and orthoDirection
            
            //calculate component of these vectors based on 2 random degrees
            direction = Normalize(direction);
            orthoDirection = Normalize(orthoDirection);
            
            float noiseA = (float)(*x);
            float noiseB = (float)(*y);
            
            //            std::cout << "ray->d.x = " << ray->d.x << std::endl;
            //            std::cout << "ray->d.y = " << ray->d.y << std::endl;
            //            std::cout << "ray->d.z = " << ray->d.z << std::endl;
            //            std::cout << " " << std::endl;
            
            //project the original ray onto the new bases
            double projA = (ray->d.x * direction.x + ray->d.y * direction.y)/sqrt(direction.x * direction.x + direction.y * direction.y);
            double projB = (ray->d.x * orthoDirection.x + ray->d.y * orthoDirection.y)/sqrt(orthoDirection.x * orthoDirection.x + orthoDirection.y * orthoDirection.y);
            double projC = ray->d.z;
            
            double rA = sqrt(projA * projA + projC * projC);
            double rB = sqrt(projB * projB + projC * projC);
            double thetaA = acos(projA/rA);
            double thetaB = acos(projB/rB);
            
            //add uncertainty
            thetaA = thetaA + noiseA;   //removed for now
            thetaB = thetaB + noiseB;
            
            //convert angles back into cartesian coordinates, but in a,b space
            double newProjA = cos(thetaA) * rA;
            
            ray->d.z = sin(thetaA) * rA;
            
            projC = ray->d.z;
            rB = sqrt(projB * projB + projC * projC);
            
            //rB = sqrt(ray->d.y * ray->d.y + ray->d.z * ray->d.z); // THIS LINE IS WRONG!!!!!
            
            //need to recalculate thetaB after you modify z (this is a new addition, also questionable)
            thetaB = acos(projB/rB);
            
            double newProjB = cos(thetaB) * rB;
            ray->d.z = sin(thetaB) * rB;
            
            //convert from ab space back to x,y space
            ray->d.x = direction.x * newProjA + orthoDirection.x * newProjB;
            ray->d.y = direction.y * newProjA + orthoDirection.y * newProjB;
            
             // TL: If NaN's in the ray, return 0 or else assert will catch.
             // TODO: Why does this happen?
                        if (ray->d.HasNaNs()) {
//                            std::cout << "ray has NaN's"<< std::endl;
                            ray->d = Vector(0,0,0);
                            return 0.f;
                        }
            
            ray->d = Normalize(ray->d);
        }
        //------------end diffraction-----------------
    }
    
    ray->o = startingPoint;
    
    ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d);
    return 1.f;
    
    
    
}
