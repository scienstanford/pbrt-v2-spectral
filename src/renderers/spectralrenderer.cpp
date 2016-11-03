
/*
 pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.
 
 This file is part of pbrt.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 
 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */


// renderers/spectralrenderer.cpp*
#include "stdafx.h"
#include "renderers/spectralrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t hash = 0, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

// SpectralRendererTask Definitions
void SpectralRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SpectralRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }
    
    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);
    
    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples]; // Total returned spectrum - all wavelengths
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];
    
    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            
            // For each sample, we loop through all the wavelength bands and trace a new ray per wavelength. We then put all the returned values in a spectrum for the original sample. (TL)
            
            // This method is slighty different from what Andy did. The advantage of this method is that when a user specifies a certain number of ray samples, that's the exact number of main rays (in a wavelength "bundle") that will be generated. In other words, there won't be the issue of having a blue image because you did not specify enough rays. The disadvantage is that in general we need more rays, since a bundle of rays with different wavelengths are all sent in the same direction and don't split until they hit the lens. In Andy's old method, each ray with its own wavelength will be shot at a different direction. In the future let's think about which way is better.
            
            // TODO: For speedup purposes, we don't actually have to trace a new ray for every single spectral sample, we can also pick a few wavelength bands and trace rays for those, and then return the values associated with the band. For example, trace a new ray for 400-450 nm, and then assign the returned spectrum from 400 to 450 nm to the final radiance at the pixel.
            
            for(int s = 0; s < nSpectralSamples; s++){
                
                Spectrum Ls_thisRay; // Returned radiance for an individual wavelength-ray
                
                // Assign a wavelength for this ray
                rays[i].wavelength = sampledLambdaStart + (sampledLambdaEnd-sampledLambdaStart)/nSpectralSamples * s;
                
                // Find camera ray for _sample[i]_
                PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
                float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
                rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
                PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);
                
                // Evaluate radiance along camera ray
                PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
                if (visualizeObjectIds) {
                    if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                        // random shading based on shape id...
                        uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                        uint32_t h = hash((char *)ids, sizeof(ids));
                        float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff),
                            float((h >> 16) & 0xff) };
                        Ls_thisRay = Spectrum::FromRGB(rgb);
                        Ls_thisRay /= 255.f;
                    }
                    else
                        Ls_thisRay= 0.f;
                }
                else {
                    if (rayWeight > 0.f)
                        // This specific ray (with an assigned wavelength band) will go through the rest of the rendering pipeline here. The included going out through the lens (where it will be refracted according to its wavelength), reflecting off objects, and finally hitting a light source. The radiance is returned here. The radiance is returned as a full spectrum, but we only care about the value associated with the ray's assigned wavelength. This is because the direction the ray exited the lens is dependent on the wavelength.
                        Ls_thisRay = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                         arena, &isects[i], &Ts[i]);
                    else {
                        Ls_thisRay = 0.f;
                        Ts[i] = 1.f;
                    }
                    
                    // Issue warning if unexpected radiance value returned
                    if (Ls_thisRay.HasNaNs()) {
                        Error("Not-a-number radiance value returned "
                              "for image sample.  Setting to black.");
                        Ls_thisRay = Spectrum(0.f);
                    }
                    else if (Ls[i].y() < -1e-5) {
                        Error("Negative luminance value, %f, returned "
                              "for image sample.  Setting to black.", Ls[i].y());
                        Ls_thisRay = Spectrum(0.f);
                    }
                    else if (isinf(Ls[i].y())) {
                        Error("Infinite luminance value returned "
                              "for image sample.  Setting to black.");
                        Ls_thisRay = Spectrum(0.f);
                    }
                }
                
                // Place the returned radiance for this wavelength into the total spectrum
                // Note, having more nSpectralSamples doesn't really affect the speed of the rendering in the rest of the pipeline, but it does affect it here because we're rendering nSpectralSamples rays per wavelength. We don't really have to worry about the speed cost of calculating the entire spectrum for each of these wavelength-dependent rays and extracting the right value (as we do below) - this calculation should be neglible speed-wise.
                float Ls_lambda;
                Ls_thisRay.GetValueAtWavelength(rays[i].wavelength, &Ls_lambda);
                Ls[i].AssignValueAtIndex(s, Ls_lambda);
                
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }
        
        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i], rays[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }
        
        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }
    
    // Clean up after _SpectralRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
                                sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// SpectralRenderer Method Definitions
SpectralRenderer::SpectralRenderer(Sampler *s, Camera *c,
                                   SurfaceIntegrator *si, VolumeIntegrator *vi,
                                   bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


SpectralRenderer::~SpectralRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void SpectralRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);
    
    // Create and launch _SpectralRendererTask_s for rendering image
    
    // Compute number of _SpectralRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new SpectralRendererTask(scene, this, camera,
                                                       reporter, sampler, sample,
                                                       visualizeObjectIds,
                                                       nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}


Spectrum SpectralRenderer::Li(const Scene *scene,
                              const RayDifferential &ray, const Sample *sample, RNG &rng,
                              MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
            Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum SpectralRenderer::Transmittance(const Scene *scene,
                                         const RayDifferential &ray, const Sample *sample, RNG &rng,
                                         MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


