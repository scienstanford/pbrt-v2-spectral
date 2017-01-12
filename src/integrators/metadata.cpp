
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


// integrators/metadata.cpp*
#include "stdafx.h"
#include "integrators/metadata.h"
#include "intersection.h"
#include "paramset.h"

// MetadataIntegrator Method Definitions
MetadataIntegrator::MetadataIntegrator(MetadataStrategy st) {
    strategy = st;
}


MetadataIntegrator::~MetadataIntegrator() {
}



Spectrum MetadataIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Intersection &isect, const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L;
    switch (strategy) {
        case MESH_MASK:
            L = Spectrum(isect.primitiveId);
            break;
            
        case MATERIAL_MASK:
            L =  Spectrum(isect.materialId);
            break;
        case DEPTH_MAP:
            L = Spectrum(ray.maxt);
            // The above works because we set:
            // r.maxt = thit in "GeometricPrimitive::Intersect"
            break;
        default:
            Error("No metadata integrator strategy specified!");
            break;
    }

    return L;
    
}


MetadataIntegrator *CreateMetadataIntegrator(const ParamSet &params) {
    MetadataStrategy strategy;
    string st = params.FindOneString("strategy", "depth"); // Can be "mesh," "material," or "depth."
    
    if (st == "mesh") strategy = MESH_MASK;
    else if (st == "material") strategy = MATERIAL_MASK;
    else if (st == "depth") strategy = DEPTH_MAP;
    else {
        Warning("Strategy \"%s\" for metadata unknown. "
                "Using \"depth\".", st.c_str());
        strategy = DEPTH_MAP;
    }
    
    return new MetadataIntegrator(strategy);
}


