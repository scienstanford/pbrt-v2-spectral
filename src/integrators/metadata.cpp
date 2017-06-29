
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
        {
            
            Point intersectPoint = isect.dg.p;
            Point rayOrigin = Point(0,0,0);
            Vector toIntersect = intersectPoint - rayOrigin;
            float distance = sqrt(toIntersect[0]*toIntersect[0] + toIntersect[1]*toIntersect[1] + toIntersect[2]*toIntersect[2]);
            L = Spectrum(distance);
            
            // DEBUG
            /*
            std::cout << "Intersect point was (" << intersectPoint[0] << "," << intersectPoint[1] << "," << intersectPoint[2] << ")" << std::endl;
            std::cout << "Ray origin was (" << ray.o[0] << "," << ray.o[1]  << "," << ray.o[2] <<  ")" << std::endl;
            std::cout << "Distance was " << distance << std::endl;
             */
            
            break;
        }
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


