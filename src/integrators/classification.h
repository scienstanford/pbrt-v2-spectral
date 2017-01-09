
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

#ifndef PBRT_INTEGRATORS_CLASSIFICATION_H
#define PBRT_INTEGRATORS_CLASSIFICATION_H

// integrators/classification.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"

// ClassificationIntegrator Declarations
class ClassificationIntegrator : public SurfaceIntegrator {
public:
    // ClassificationIntegrator Public Methods
    ClassificationIntegrator();
    ~ClassificationIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
};


ClassificationIntegrator *CreateClassificationIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_CLASSIFICATION_H
