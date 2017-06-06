
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

#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

// core/material.h*
#include "pbrt.h"
#include "memory.h"

// Material Declarations
class Material : public ReferenceCounted {
public:
    // Material Interface
    Material() : materialId(nextmaterialId++) {} // Added by Trisha
    virtual BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                          const DifferentialGeometry &dgShading,
                          MemoryArena &arena) const = 0;
    virtual BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                              const DifferentialGeometry &dgShading,
                              MemoryArena &arena) const {
        return NULL;
    }
    virtual ~Material();
    static void Bump(const Reference<Texture<float> > &d, const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, DifferentialGeometry *dgBump);
    
    static void NormalMap(const Reference<Texture<Spectrum> > &d, const DifferentialGeometry &dgGeom,
                       const DifferentialGeometry &dgShading, DifferentialGeometry *dgBump);
    // Material Public Data
    const uint32_t materialId; // Added by Trisha
protected:
    // Material Protected Data
    static uint32_t nextmaterialId; // Added by Trisha
};



#endif // PBRT_CORE_MATERIAL_H
