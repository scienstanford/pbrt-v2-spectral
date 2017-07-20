
#include "stdafx.h"
#include "cameras/realisticEye.h"
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
#include <gsl/gsl_roots.h> // For solving biconic surface intersections.
#include <gsl/gsl_errno.h>

using namespace std;

// -----------------------------------------
// Needed for solving intersection of ray with biconic surface
// -----------------------------------------
struct biconic_params {
    float Rx, Ry, Cx, Cy;
    Ray ray;
};
double BiconicSag(double t, void *params){
    
    struct biconic_params *p;
    p = (struct biconic_params *)params;
    
    Point intersect = p->ray(t);
    float x,y,z;
    z = intersect.z;
    x = intersect.x;
    y = intersect.y;
    
    float f,g,g_term;
    f = (x*x)/p->Rx + (y*y)/p->Ry;
    g_term = 1 - (1+p->Cx)*(x*x)/(p->Rx*p->Rx) - (1+p->Cy)*(y*y)/(p->Ry*p->Ry);
    
    if(g_term < 0){
        // TODO: What to do here? Let's just make it a small number to ensure that f/g =/= z.
        g_term = 0.001;
    }
    
    g = 1 + sqrt(g_term);
    
    return z-f/g;
}

// -----------------------------------------
// -----------------------------------------


RealisticEyeCamera *CreateRealisticEyeCamera(const ParamSet &params,
                                                             const AnimatedTransform &cam2world, Film *film) {
    
    // Generate warning about lens radius convention.
    Warning("For the realistic eye code, we have gone with the Zemax convention for the lens radius. A positive lens radius means the center of the spherical lens is in the positive direction (toward the scene) relative to the location of the lenes. This is different than the other camera classes where the convention is flipped. If necessary, one can turn on the [flipLensRadius] flag.");
    
    // -----------------------------------
    // ----- Common Camera Parameters ----
    // -----------------------------------
    
    // Even though we don't need these, the parent Camera class does.
    float hither = params.FindOneFloat("hither", -1); // Do we need this?
    float yon = params.FindOneFloat("yon", -1); // Do we need this?
    float shutteropen = params.FindOneFloat("shutteropen", -1);
    float shutterclose = params.FindOneFloat("shutterclose", -1);

    
    // --------------------------
    // ----- Eye Parameters -----
    // --------------------------
    
    // The lens file should include  columns for [radiusX radiusY thickness mediumIndex semiDiameter conicConstantX conicConstantY].
    
    // The surfaces should be ordered from the world to the retina. E.g. The last surface in the text file is the retina.
    // For the radius, positi CHECK THIS!
    
    // Check for lens file
    string specfile = params.FindOneString("specfile", "");
    if (specfile == "") {
        Severe( "No lens spec file supplied!\n" );
    }
    
    // These are additional parameters we need to specify.
    float pupDiam = params.FindOneFloat("pupilDiameter", 4.0); //mm
    float lensDecenterX = params.FindOneFloat("lensDecenterX",0); //mm
    float lensDecenterY = params.FindOneFloat("lensDecenterY",0); //mm
    float lensTiltX = params.FindOneFloat("lensTiltX",0); //degrees
    float lensTiltY = params.FindOneFloat("lensTiltY",0); //degrees
    float retinaDistance = params.FindOneFloat("retinaDistance",0); //mm
    float retinaRadius = params.FindOneFloat("retinaRadius",0); //mm
    float retinaSemiDiam = params.FindOneFloat("retinaSemiDiam",0); //mm
    
    // Check for IORspectra slots
    Spectrum ior1 = params.FindOneSpectrum("ior1", 0);
    Spectrum ior2 = params.FindOneSpectrum("ior2", 0);
    Spectrum ior3 = params.FindOneSpectrum("ior3", 0);
    Spectrum ior4 = params.FindOneSpectrum("ior4", 0);
    Spectrum ior5 = params.FindOneSpectrum("ior5", 0);
    Spectrum ior6 = params.FindOneSpectrum("ior6", 0);
    
    vector<Spectrum> iorSpectra;
    iorSpectra.push_back(ior1);
    iorSpectra.push_back(ior2);
    iorSpectra.push_back(ior3);
    iorSpectra.push_back(ior4);
    iorSpectra.push_back(ior5);
    iorSpectra.push_back(ior6);

    // Specify which material index corresponds to the GRIN lens. If nothing is specified (grinSurfaceIndex == -1) we assume that a normal lens is used instead of a GRIN lens
    int grinSurfaceIndex = params.FindOneInt("grinMaterialIndex", -1);
    string grinSurfaceFile = params.FindOneString("GRINfile", "");
    
    // Flags for realism/speed
    bool chromaticFlag = params.FindOneBool("chromaticAberrationEnabled", 0.0);
    
    // Flags for convention
    bool flipRadFlag = params.FindOneBool("flipLensRadius", 0.0);
    
    return new RealisticEyeCamera(cam2world,film,hither,yon,shutteropen,shutterclose,specfile,pupDiam,lensDecenterX,lensDecenterY,lensTiltX,lensTiltY,retinaDistance,retinaRadius,retinaSemiDiam,chromaticFlag,flipRadFlag,iorSpectra,grinSurfaceIndex,grinSurfaceFile);
}




RealisticEyeCamera::RealisticEyeCamera(const AnimatedTransform &cam2world,
                                       Film *f,
                                       float hither, float yon,
                                       float sopen, float sclose,
                                       string specfile,
                                       float pupDiam,
                                       float ldX,
                                       float ldY,
                                       float ltX,
                                       float ltY,
                                       float rD,
                                       float rR,
                                       float rSD,
                                       bool chromaticFlag,
                                       bool flipRadFlag,
                                       vector<Spectrum> iorS,
                                       int gSI,
                                       string grinFile)
: Camera(cam2world, sopen, sclose, f), ShutterOpen(sopen), ShutterClose(sclose),film(f)
{
    
    // Find the complete path for the specfile
    string lensFileName = AbsolutePath(ResolveFilename(specfile));
    string grinFileName = AbsolutePath(ResolveFilename(grinFile));
    
    pupilDiameter = pupDiam;
    lensDecenterX = ldX;
    lensDecenterY = ldY;
    lensTiltX = ltX;
    lensTiltY = ltY;
    retinaDistance = rD;
    retinaRadius = rR;
    retinaSemiDiam = rSD;
    iorSpectra = iorS;
    grinSurfaceIndex = gSI;
    if(grinSurfaceIndex != -1) grinLensFlag = true;
    
    chromaticAberrationEnabled = chromaticFlag;
    flipLensRadius = flipRadFlag;
    
    // -------------------------
    // --- Read in lens file ---
    // -------------------------

    vector<float> vals;
    
    // Check to see if there is valid input in the lens file.
    if (!ReadFloatFile(lensFileName.c_str(), &vals)) {
        Warning("Unable to read lens file!");
        return;
    }
    
    // The lens file should include  columns for [radiusX radiusY thickness mediumIndex semiDiameter conicConstantX conicConstantY].
    // Let's check then that the file is a multiple of 7 (not including the effective focal length at the top)
    if ((vals.size()-1) % 7 != 0)
    {
        Warning("Wrong number of float values in lens file! Did you forget to specify the focal length? Is this a lens file with biconic surfaces? Do you have a carriage return at the end of the data?");
        return;
    }
    
    effectiveFocalLength = vals[0];   // Read the effective focal length.
    
    for (int i = 1; i < vals.size(); i+=7)
    {
        LensElementEye currentLensEl;
        currentLensEl.radiusX = vals[i];
        currentLensEl.radiusY = vals[i+1];
        currentLensEl.thickness = vals[i+2];
        currentLensEl.mediumIndex = vals[i+3];
        currentLensEl.semiDiameter = vals[i+4];
        currentLensEl.conicConstantX = vals[i+5];
        currentLensEl.conicConstantY = vals[i+6];
        
        // Note: Zemax and PBRT-spectral seem to have different conventions for what the positive and negative sign of the radius is. In Zemax, a positive radius means that the center of lens sphere is directed toward the positive Z-axis, and vice versa. In previous PBRT-spectral iterations, this was flipped. Here I've rewritten the lens tracing code to go with the Zemax convention, however to be backward compatible we might want to have this ability to flip the radii.
        if(flipLensRadius){
            currentLensEl.radiusX = -1*currentLensEl.radiusX;
            currentLensEl.radiusY = -1*currentLensEl.radiusY;
            currentLensEl.conicConstantX = -1*currentLensEl.conicConstantX;
            currentLensEl.conicConstantY = -1*currentLensEl.conicConstantY;
            Warning("Flipping lens radius & conic convention.");
        }
        
        // A radius of zero in BOTH x and y directions indicates an aperture. We should be careful of this though, since sometimes we may want to define a flat surface...
        // If the surface is an aperture, we set it's size to be equal to the pupil diameter specified.
        if (currentLensEl.radiusX == 0 && currentLensEl.radiusY == 0 ){
            currentLensEl.semiDiameter = pupilDiameter/2;
        }
        else{
            // We have to do a semi-diameter check here. As we change accommodation, we also change the radius of curvature, and we don't want the semi-diameter to be bigger than the radius. This manifests as a square root of a negative number in equation 1 on Einighammer et al. 2009.
            // TODO: This check is sort of hack-y, is there a better mathematical way to do this?
            float smallerR = min(currentLensEl.radiusX,currentLensEl .radiusY);
            float biggerK = max(currentLensEl.conicConstantX,currentLensEl.conicConstantY);
//            if((currentLensEl.semiDiameter*currentLensEl.semiDiameter) >= (smallerR*smallerR/(1+biggerK))){
//                Warning("Changing semidiameter of an element to match radius/asphericity geometry.");
//                currentLensEl.semiDiameter = 0.95 * sqrt((smallerR*smallerR/(1+biggerK))); // 0.95 is to add some buffer zone, since rays act very strangely when they get too close to the edge of the conical surface.
//            }
            if(currentLensEl.semiDiameter*currentLensEl.semiDiameter*(1+biggerK)/(smallerR*smallerR) > 1.0f ){
                currentLensEl.semiDiameter = 0.95 * sqrt((smallerR*smallerR/(1+biggerK))); // 0.95 is to add some buffer zone, since rays act very strangely when they get too close to the edge of the conical surface.
            }
        }
        
        
        lensEls.push_back(currentLensEl);
    }
    
    // Check thickness of last element. It should be zero, since we use the "retina distance" parameter for this final "thickness."
    if(lensEls[lensEls.size()-1].thickness != 0){
        Error("Thickness of lens element closest to zero must be zero. Define thickness in 'retinaDistance' parameter instead.");
    }
    // ---------------------------
    // --- Read in GRIN file ---
    // ---------------------------
    
    // TODO: DOUBLE CHECK THAT THIS IS WORKING
    
    /* The GRIN file described the GRIN lens used. It follows the same format as Zemax, specifically:
     
     MIN_WAVELENGTH MAX_WAVELENGTH
     
     REF_WAVELENGTH
     
     K_MAX L_MAX
     
     K11 K12 K13 ... K1K_MAX
     
     K21 K22 K23 ... K2K_MAX
     
     K31 K32 K33 ... K3K_MAX
     
     L11 L12 L13 ... L1L_MAX
     
     L21 L22 L23 ... L2L_MAX
     
     L31 L32 L33 ... L3L_MAX
     
     Wavelength should be in um. Total number of values should be K_MAX*3 + L_MAX*3 + 3
     
     */
    
    
    if(grinSurfaceIndex != -1){
        
        vals.clear();
        
        // Check to see if there is valid input in the lens file.
        if (!ReadFloatFile(grinFileName.c_str(), &vals)) {
            Warning("Unable to read GRIN file!");
            return;
        }
        
        GRINelement grin;
        float minWavelength = vals[0]*1000; // We operate in nm
        float maxWavelength = vals[1]*1000;
        float refWavelength = vals[2]*1000;
        float K_max = (int)vals[3];
        float L_max = (int)vals[4];
        
        if(vals.size() != (K_max*3+L_max*3+3)){
            Warning("Wrong number of float values in lens file! Maybe check the format again?");
        }
        grin.K1 = 0; grin.K2 = 0; grin.K3 = 0;
        grin.L1 = 0; grin.L2 = 0; grin.L3 = 0;
        for(int i = 0; i < K_max; i++){
            float wavePow = pow(refWavelength,(i+1)-1);
            grin.K1 += vals[5+i]*wavePow;
            grin.K2 += vals[5+K_max+i]*wavePow;
            grin.K3 += vals[5+2*K_max+i]*wavePow;
        }
        

        
        
    }
    
    // To calculate the "film diagonal", we use the retina semi-diameter. The film diagonal is the diagonal of the rectangular image rendered out by PBRT, in real units. Since we restrict samples to a circular image, we can calculate the film diagonal to be the same as a square that circumscribes the circular image.
    filmDiag = retinaSemiDiam*1.4142*2; // sqrt(2)*2
    
    // We are going to use our own error handling for gsl to prevent it from crashing the moment a ray doesn't intersect.
    gsl_set_error_handler_off();
    
    
}


RealisticEyeCamera::~RealisticEyeCamera()
{
}



float RealisticEyeCamera::getFocalLength()
{
    return effectiveFocalLength;
}

void RealisticEyeCamera::applySnellsLaw(float n1, float n2, float lensRadius, Vector &normalVec, Ray * ray ) const
{
    
    //                  |
    // scene ..... n2   |   n1 ..... sensor
    //                  |
    //                  |
    //                 lens
    
    
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

bool RealisticEyeCamera::IntersectLensElAspheric(const Ray &r, float *tHit, LensElementEye currElement, float zShift, Vector *n) const{
    
    // This code lets us find the intersection with an aspheric surface. At tHit, the ray will intersect the surface. Therefore:
    // If
    // (x,y,z) = ray_origin + thit * ray_direction
    // then:
    // z - u(x,y) = 0
    // where u(x,y) is the SAG of the surface, defined in Eq.1 of Einighammer et al. 2009 and in the Zemax help page under biconic surfaces.
    // We can use this fact to solve for thit. If the surface is not a sphere, this is a messy polynomial. So instead we use a numeric root-finding method (Van Wijingaarden-Dekker-Brent's Method.) This method is available in the GSL library.
    
    // DEBUG
    /*
    std::cout << "r.o = " << r.o.x << "," << r.o.y << "," << r.o.z << std::endl;
    std::cout << "r.d = " << r.d.x << "," << r.d.y << "," << r.d.z << std::endl;
    */
    
    // Move ray to object(lens) space.
    Ray objSpaceRay = r;
    objSpaceRay.o = objSpaceRay.o + Vector(0,0,zShift);
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double root = 0;
    double x_lo = 0.0;
    double x_hi;
    if(currElement.thickness == 0){
        // Probably the surface closest to the retina.
        x_hi = retinaDistance*2;
    }else{
        x_hi = currElement.thickness*1.5; // thit will probably be less than this
    }
    gsl_function F;
    
    // DEBUG
    /*
    std::cout << "radiusX = " << currElement.radiusX << std::endl;
    std::cout << "radiusY = " << currElement.radiusY << std::endl;
    std::cout << "conicConstantX = " << currElement.conicConstantX << std::endl;
    std::cout << "conicConstantY = " << currElement.conicConstantY << std::endl;
    std::cout << "objSpaceRay.o = " << objSpaceRay.o.x << "," << objSpaceRay.o.y << "," << objSpaceRay.o.z << std::endl;
    std::cout << "objSpaceRay.d = " << objSpaceRay.d.x << "," << objSpaceRay.d.y << "," << objSpaceRay.d.z << std::endl;
    */
    
    struct biconic_params params = {currElement.radiusX,currElement.radiusY,currElement.conicConstantX,currElement.conicConstantY,objSpaceRay};
    F.function = &BiconicSag;
    F.params = &params;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);

    status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    if(status != 0){
        // Ray probably does not intersect. This might depend on the x_hi set above, i.e. if it's too small OR too large. TODO: Can we check this?
        return false;
    }
    
    // DEBUG
    /*
    printf ("using %s method\n",
            gsl_root_fsolver_name (s));
    
    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
            "iter", "lower", "upper", "root",
            "err", "err(est)");
    */
    
    do
    {
        iter++;
        gsl_root_fsolver_iterate (s);
        root = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.0001);

        if (status == GSL_SUCCESS){
 
            // DEBUG
            /*
            printf ("Converged:\n");
            printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                    iter, x_lo, x_hi,
                    root, root - root,
                    x_hi - x_lo);
            */
            
            gsl_root_fsolver_free (s);
            
            // DEBUG
            /*
            std::cout << "root = " << root << std::endl;
            std::cout << "r.o = " << r.o.x << "," << r.o.y << "," << r.o.z << std::endl;
            std::cout << "r.d = " << r.d.x << "," << r.d.y << "," << r.d.z << std::endl;
            */
            
            *tHit = root;
            Point intersect = r(*tHit);
            
            // Check if intersection is within the semi-diameter of the lens, if not we return false (no intersect)
            // (If we don't do this here, we might get a complex normal which would crash the rendering.)
            if(intersect.x * intersect.x + intersect.y * intersect.y > (currElement.semiDiameter * currElement.semiDiameter)){
                return false;
            }
            
            // Calculate normal at intersection
            // These equations are from Eq 2, 3, and 4 in Einighammer 2009
            float term1 = ((1+currElement.conicConstantX)*intersect.x*intersect.x)/(currElement.radiusX*currElement.radiusX);
            float term2 = ((1+currElement.conicConstantY)*intersect.y*intersect.y)/(currElement.radiusY*currElement.radiusY);
            
            float fprime_x = 2*intersect.x/currElement.radiusX;
            float gprime_x = (-1*(1+currElement.conicConstantX)*intersect.x)/(currElement.radiusX*currElement.radiusX*sqrt(1-term1-term2));
            
            float fprime_y = 2*intersect.y/currElement.radiusY;
            float gprime_y = (-1*(1+currElement.conicConstantY)*intersect.y)/(currElement.radiusY*currElement.radiusY*sqrt(1-term1-term2));
            
            float f = (intersect.x*intersect.x)/currElement.radiusX + (intersect.y*intersect.y)/currElement.radiusY;
            float g = 1+sqrt(1-term1-term2);
            
            float zprime_y = (fprime_y*g-gprime_y*f)/(g*g);
            float zprime_x = (fprime_x*g-gprime_x*f)/(g*g);
            
            Vector v_x = Vector(1,0,zprime_x);
            Vector v_y = Vector(0,1,zprime_y);
            
            *n = Normalize(Cross(v_x,v_y));
            *n = Faceforward(*n, -r.d);
            
            return true;
            
        }
        
        // For debugging
        /*
        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                iter, x_lo, x_hi,
                root, root - root,
                x_hi - x_lo);
        */
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free (s);
    
    return false;
}

bool RealisticEyeCamera:: IntersectLensEl(const Ray &r, float *tHit, float radius, Vector dist, Vector & normalVec){
    
    // TL: This is the same code that is in sphere.cpp
    // We have to be careful about where the zero axis is in all of this.
    // TODO: Double check Andy's calculations. For now I'm assuming its correct.
    
    float phi;
    Point phit;
    
    // Transform _Ray_ to object space
    Transform shiftZ =  Translate(dist);
    Ray ray = r;
    (shiftZ)(r, &ray);
    
    
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

float RealisticEyeCamera::getSensorWidth()
{
    
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    return width;
}

float RealisticEyeCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
    // DEBUG
//    std::cout << "" << std::endl;
    
    // TL:                     z=0         z= -filmDistance
    //                  |  |    ||           |
    //                  |  |    ||           |
    //  Scene <---------|--|----||---------> |
    //            +z    |  |    ||    -z     |
    //                  |  |    ||           |
    //                Lens Elements        Sensor
    //
    
    
    // GenerateRay() should return the weight of the generated ray
    
    // Determine the size of the sensor in real world units (i.e. convert from pixels to millimeters).
    float aspectRatio = (float)film->xResolution/(float)film->yResolution;
    float width = filmDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
    float height = width/aspectRatio;
    
    Point startingPoint;
    startingPoint.x = -((sample.imageX) - film->xResolution/2.f - .25)/(film->xResolution/2.f);
    startingPoint.y = ((sample.imageY) - film->yResolution/2.f - .25)/(film->yResolution/2.f);
    
    // Convert starting point units to millimeters
    startingPoint.x = startingPoint.x * width/2.f;
    startingPoint.y = startingPoint.y * height/2.f;
    startingPoint.z = -retinaDistance;
    
    if (retinaRadius != 0)
    {
        // Right now the code only lets you curve the sensor toward the scene and not the other way around. See diagram:
        /*    
        
         The distance between the zero point on the z-axis (i.e. the lens element closest to the sensor) and the dotted line will be equal to the "retinaDistance." The retina curvature is defined by the "retinaRadius" and it's height in the y and x direction is defined by the "retinaSemiDiam."
         
         
                                        :
                                     |  :
                                      | :
                         | | |         |:
           scene <------ | | | <----   |:
                         | | |         |:
                     Lens System      | :
                                     |  :
                                        :
                                  retina
            <---- +z
         
         */

        // Limit sample points to a circle within the retina semi-diameter
        if((startingPoint.x*startingPoint.x + startingPoint.y*startingPoint.y) > (retinaSemiDiam*retinaSemiDiam)){
            return 0.f;
        }
        
        // Calculate the distance of a disc that fits inside the curvature of the retina.
        float zDiscDistance = -1*sqrt(retinaRadius*retinaRadius-retinaSemiDiam*retinaSemiDiam);
        
        // If we are within this radius, project each point out onto a sphere. There may be some issues here with even sampling, since this is a direct projection...
        float el = atan(startingPoint.x/zDiscDistance);
        float az = atan(startingPoint.y/zDiscDistance);
        
        // Convert spherical coordinates to cartesian coordinates (note: we switch up the x,y,z axis to match our conventions)
        float xc,yc,zc, rcoselev;
        xc = -1*retinaRadius*sin(el); // TODO: Confirm this flip?
        rcoselev = retinaRadius*cos(el);
        zc = -1*(rcoselev*cos(az)); // The -1 is to account for the curavure described above in the diagram
        yc = -1*rcoselev*sin(az); // TODO: Confirm this flip?
        
        zc = zc + -1*retinaDistance + retinaRadius; // Move the z coordinate out to correct retina distance
        
        startingPoint = Point(xc,yc,zc);
        
    }
    
    float lensU, lensV;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    
    //We need to shoot rays toward the disc that fits inside the curvature of first lens surface. Since we no longer have spherical elements, we have to calculate it as follows:
    // TODO: This is just a guess. Is there a correct way to do this? Should we look at the exitPupil code in PBRTv3?
    
    float lensP_semiDiam = lensEls[lensEls.size()-1].semiDiameter;
    float lensP_radius = lensEls[lensEls.size()-1].radiusX;
    
    // sgn(lensP_radius)
    // It's very rare for the first lens element to have a negative radius (spherical center toward sensor)...but just in case:
    float sgn_radius = (lensP_radius > 0) - (lensP_radius < 0);
    float discDistance = sgn_radius * BiconicZ(lensP_semiDiam, 0, lensEls[lensEls.size()-1]);
    
    // Scale the normalized lens coordinates by the size of the first lens element
    lensU *= lensP_semiDiam;
    lensV *= lensP_semiDiam;
    
    Point pointOnLens = Point(lensU, lensV, discDistance);   // We aim the ray at a flat disk, will that cause problems later?
    
    float tempWavelength = ray->wavelength;
    *ray = Ray(Point(0,0,0), Normalize(Vector(startingPoint)), 0.f, INFINITY);
    ray->o = startingPoint;    //initialize ray origin
    ray->d = Normalize(pointOnLens - ray->o);
    ray->wavelength = tempWavelength;  //so that wavelength information is retained.
    
    
    // DEBUG
    /*
    // Scene to retina
    startingPoint = Point(0,0.00081193431816,-16.319999973);
    pointOnLens = Point(0,1.544122632,0.21476753379);
    ray->o = startingPoint;
    ray->d = Normalize(pointOnLens - ray->o);
    ray->wavelength = 550;
    
    // Retina to scene
    startingPoint = Point(0,0,-16.3200);
    pointOnLens = Point(0,1.5330684279,0.2229);
    ray->o = startingPoint;
    ray->d = Normalize(pointOnLens - ray->o);
    ray->wavelength = 550;
    */
    
    // --------------------------------------------------------
    // --- Trace through the lens elements of the main lens ---
    // --------------------------------------------------------
    
    float lensDistance = 0; // How far we are from the "0" on the z-axis (see diagram above)
    
    for (int i = lensEls.size()-1; i>=0 ; i--)
    {
        
        ray->o = startingPoint;
        lensDistance += lensEls[i].thickness;
        
        // If the ray direction is zero, there is probably internal reflection going on somewhere. We will just terminate the ray here to avoid assert errors.
        if(ray->d == Vector(0,0,0)){
            return 0.f;
        }
        
        // DEBUG
        // ----
        /*
        std::cout << "\n" << std::endl;
        std::cout << "i = " << i << std::endl;
        std::cout << "start " << ray->o.x << " " << ray->o.y << " " << ray->o.z << std::endl;
        std::cout << "dir " << ray->d.x << " " << ray->d.y << " " << ray->d.z << std::endl;
         */
        // ----
        
        float tHit = 0;
        bool intersected = false;
        Vector normalVec(0,0,1);
        Point intersectPoint(0,0,0);
        
        if (lensEls[i].radiusX == 0 && lensEls[i].radiusY == 0)
        {
            // ---------------------
            // --- APERTURE CASE ---
            // ---------------------
            
            float tAperture = 0;
            if (i == lensEls.size()-1)
                tAperture = retinaDistance/ray->d.z;   //special case for when aperture is the 1st element
            else
                tAperture = (lensDistance - ray->o.z)/(ray->d.z);
            
            // Point where the ray intersects the aperture plane
            Point intersectPoint = (*ray)(tAperture);
            normalVec = Vector(0,0,1);

            // Move ray to start at intersection
            startingPoint = intersectPoint;
            
            // Check if ray makes it through the aperture
            if((intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y) > (lensEls[i].semiDiameter * lensEls[i].semiDiameter)){
                return 0.f;
            }
            
            
        }
        else
        {
            // ----------------------------
            // --- REGULAR ELEMENT CASE ---
            // ----------------------------
            
            //---- Find intersection point ----
            
            // Since the surface is defined in object space with the edge at zero, so we need to move the ray accordingly from world space to object space. This is simply a translation according to the location of the surface.
            
            // z(x,y) is defined like this (flipped across the z axis for +r vs -r):
            //
            //   :   /
            //   :  /
            //   : /
            //   :|
            //   :|
            //   :|
            //   : \d_zema
            //   :  \
            //   :   \
            //  z=0
            //
            
            float zShift = -lensDistance;
            
            // Does the ray intersect the lens surface? If so, where?
            intersected = IntersectLensElAspheric(*ray, &tHit, lensEls[i], zShift, &normalVec);
            
            if (intersected)
            {
                intersectPoint = (*ray)(tHit);
                
                //DEBUG
               // std::cout << "Intersect point: " << intersectPoint.x << " " << intersectPoint.y << " " << intersectPoint.z << std::endl;
                
                
                // ---- Apply Snell's Law ----
                
                // TOOD: Check this configuration!
                //                  |
                // scene ..... n2   |   n1 ..... sensor
                //           [i-1]  |  [i]
                //                  |
                //                 lens
                //
                
                // To trace wavelength by wavelength through the lens, the spectral renderer must be selected. If it's not, let's warn the user.
                if(chromaticAberrationEnabled && ray->wavelength == 0){
                    Error("Chromatic aberration enabled but the ray has no wavelength. Are you using the spectral renderer?");
                }
                
                // The user can load IOR spectra for each ocular medium into ior1, ior2, etc. In the lens file, they can then specify with medium they would like to use. The number (X) corresponds to iorX.
                
                float n1,n2;
                
                n1 = lookUpIOR(lensEls[i].mediumIndex, *ray);
                
                // If we're at the lens surface closest to the scene, n2 should be air.
                if (i-1 >= 0){
                    
                    n2 = lookUpIOR(lensEls[i-1].mediumIndex,*ray);
                    
                    // Trisha: If we're entering the aperture (sometimes we put n2 == 0 when that is the case) we skip the n2 == 0 and apply the next medium.
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
                        n2 = lookUpIOR(lensEls[i-2].mediumIndex, *ray);
                    
                }
                else{
                    n2 = 1;
                }
                
                applySnellsLaw( n1,  n2,  0, normalVec, ray );
                // --- Update ray starting point ---
                
                startingPoint = intersectPoint;

            }
            else
            {
                return 0.f;
            }
        }
        
        
    }
    
    ray->o = startingPoint;
    
    // DEBUG
    /*
    // ----
    std::cout << "\n" << std::endl;
    std::cout << "final" << std::endl;
    std::cout << ray->o.x << " " << ray->o.y << " " << ray->o.z << std::endl;
    std::cout << ray->d.x << " " << ray->d.y << " " << ray->d.z << std::endl;
    // ----
    */
    
    ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d);
    return 1.f;
    
    
    
}

// Handy method to explicity solve for the z(x,y) at a given point (x,y),for the biconic SAG.
float RealisticEyeCamera::BiconicZ(float x, float y, LensElementEye currElement) const{
    
    float f,g,g_term;
    float Rx,Ry,Cx,Cy;
    Rx = currElement.radiusX; Ry = currElement.radiusY;
    Cx = currElement.conicConstantX; Cy = currElement.conicConstantY;
    
    f = (x*x)/Rx + (y*y)/Ry;
    g_term = 1 - (1+Cx)*(x*x)/(Rx*Rx) - (1+Cy)*(y*y)/(Ry*Ry);
    
    if(g_term < 0){
        g_term = 0.001;
        Warning("Encountered a complex value when solving for the z(x,y) of the biconic.");
    }
    
    g = 1 + sqrtf(g_term);
    
    return f/g;
    
}

// Given the mediumIndex, load up the right spectra from the ones read in through ior1, ior2, etc. Then find the corresponding IOR for the given ray wavelength.
float RealisticEyeCamera::lookUpIOR(int mediumIndex, const Ray &ray)const{
    
    float n;
    
    if(mediumIndex == grinSurfaceIndex){
        // A GRIN lens
        
    }else{
        // Standard media
        if(chromaticAberrationEnabled){
            iorSpectra[mediumIndex-1].GetValueAtWavelength(ray.wavelength,&n);
        }else{
            iorSpectra[mediumIndex-1].GetValueAtWavelength(550,&n);
        }
    }
    return n;
}


