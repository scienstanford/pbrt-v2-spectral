
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_SPECTRALRENDERER_H
#define PBRT_RENDERERS_SPECTRALRENDERER_H

#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"

// This renderer is virtually identical to samplerrenderer, the only difference is that is traces one ray per wavelength. We do this in order to create wavelength-dependent lens artifacts such as chromatic aberration. (In the realistic diffraction lenses, the wavelengths are bent according to their wavelength and Snell's law.) Each pixel sample is therefore multiplied by the number of wavelength bands. The speed of this renderer will be many times slower than that of the sampler renderer.
// Trisha Lian 10/11/2016

// SpectralRenderer Declarations
class SpectralRenderer : public Renderer {
public:
    // SpectralRenderer Public Methods
    SpectralRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, bool visIds, int numWaves);
    ~SpectralRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // SpectralRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
    int nWaveBands;
};



// SpectralRendererTask Declarations
class SpectralRendererTask : public Task {
public:
    // SpectralRendererTask Public Methods
    SpectralRendererTask(const Scene *sc, Renderer *ren, Camera *c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        bool visIds, int tn, int tc, int numWave)
      : reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
        nWaveBands = numWave;
    }
    void Run();
private:
    // SpectralRendererTask Private Data
    const Scene *scene;
    const Renderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
    int nWaveBands;
};



#endif // PBRT_RENDERERS_SPECTRALRENDERER_H
