
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

#ifndef PBRT_MATERIALS_MIRROR_H
#define PBRT_MATERIALS_MIRROR_H

// materials/mirror.h*
#include "pbrt.h"
#include "material.h"

// MirrorMaterial Declarations
class MirrorMaterial : public Material {
public:
    // MirrorMaterial Public Methods
    MirrorMaterial(Reference<Texture<Spectrum> > r, Reference<Texture<float> > bump, Reference<Texture<Spectrum> > normal) {
        Kr = r;
        bumpMap = bump;
        normalMap = normal;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // MirrorMaterial Private Data
    Reference<Texture<Spectrum> > Kr, normalMap;
    Reference<Texture<float> > bumpMap;
};


MirrorMaterial *CreateMirrorMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_MIRROR_H
