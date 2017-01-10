
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


// integrators/classification.cpp*
#include "stdafx.h"
#include "integrators/classification.h"
#include "intersection.h"
#include "paramset.h"

// ClassificationIntegrator Method Definitions
ClassificationIntegrator::ClassificationIntegrator(ClassificationStrategy st) {
    strategy = st;
}


ClassificationIntegrator::~ClassificationIntegrator() {
}



Spectrum ClassificationIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Intersection &isect, const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L;
    switch (strategy) {
        case CLASSIFY_BY_MESH:
            L = Spectrum(isect.primitiveId);
            break;
            
        case CLASSIFY_BY_MATERIAL:
            L =  Spectrum(isect.materialId);
            break;
    }
    
    return L;
    
}


ClassificationIntegrator *CreateClassificationIntegrator(const ParamSet &params) {
    ClassificationStrategy strategy;
    string st = params.FindOneString("strategy", "mesh"); // Can be Mesh or Material
    
    if (st == "mesh") strategy = CLASSIFY_BY_MESH;
    else if (st == "material") strategy = CLASSIFY_BY_MATERIAL;
    else {
        Warning("Strategy \"%s\" for classification unknown. "
                "Using \"mesh\".", st.c_str());
        strategy = CLASSIFY_BY_MESH;
    }
    
    return new ClassificationIntegrator(strategy);
}


