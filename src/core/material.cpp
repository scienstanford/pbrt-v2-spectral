
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


// core/material.cpp*
#include "stdafx.h"
#include "material.h"
#include "primitive.h"
#include "texture.h"
#include "spectrum.h"
#include "reflection.h"

// Material Method Definitions
uint32_t Material::nextmaterialId = 1; //Added by Trisha
Material::~Material() {
}


void Material::Bump(const Reference<Texture<float> > &d,
                    const DifferentialGeometry &dgGeom,
                    const DifferentialGeometry &dgs,
                    DifferentialGeometry *dgBump) {
    // Compute offset positions and evaluate displacement texture
    DifferentialGeometry dgEval = dgs;
    
    // Shift _dgEval_ _du_ in the $u$ direction
    float du = .5f * (fabsf(dgs.dudx) + fabsf(dgs.dudy));
    if (du == 0.f) du = .01f;
    dgEval.p = dgs.p + du * dgs.dpdu;
    dgEval.u = dgs.u + du;
    dgEval.nn = Normalize((Normal)Cross(dgs.dpdu, dgs.dpdv) +
                          du * dgs.dndu);
    float uDisplace = d->Evaluate(dgEval);

    // Shift _dgEval_ _dv_ in the $v$ direction
    float dv = .5f * (fabsf(dgs.dvdx) + fabsf(dgs.dvdy));
    if (dv == 0.f) dv = .01f;
    dgEval.p = dgs.p + dv * dgs.dpdv;
    dgEval.u = dgs.u;
    dgEval.v = dgs.v + dv;
    dgEval.nn = Normalize((Normal)Cross(dgs.dpdu, dgs.dpdv) +
                          dv * dgs.dndv);
    float vDisplace = d->Evaluate(dgEval);
    float displace = d->Evaluate(dgs);

    // Compute bump-mapped differential geometry
    *dgBump = dgs;
    dgBump->dpdu = dgs.dpdu + (uDisplace - displace) / du * Vector(dgs.nn) +
                   displace * Vector(dgs.dndu);
    dgBump->dpdv = dgs.dpdv + (vDisplace - displace) / dv * Vector(dgs.nn) +
                   displace * Vector(dgs.dndv);
    dgBump->nn = Normal(Normalize(Cross(dgBump->dpdu, dgBump->dpdv)));
    
    if (dgs.shape->ReverseOrientation ^ dgs.shape->TransformSwapsHandedness)
        dgBump->nn *= -1.f;

    // Orient shading normal to match geometric normal
    dgBump->nn = Faceforward(dgBump->nn, dgGeom.nn);
    
    
    
}

void Material::NormalMap(const Reference<Texture<Spectrum> > &d,
                    const DifferentialGeometry &dgGeom,
                    const DifferentialGeometry &dgs,
                    DifferentialGeometry *dgBump) {
    
    DifferentialGeometry dgEval = dgs;
    
    
    RGBSpectrum normalSpectrum = d->EvaluateMemory(dgs); // This read the RGB values directly from the texture, without first converting to spectra.
    float normalRGB[3];
    normalSpectrum.ToRGB(normalRGB);

    // Remap normals
    // X: -1 to +1 :  Red: 0 to 1
    // Y: -1 to +1 :  Green: 0 to 1
    // Z: 0 to -1 :  Blue: 0.5 to 1
    normalRGB[0] = (normalRGB[0])*2 -1;
    normalRGB[1] = (normalRGB[1])*2 -1;
    normalRGB[2] = (normalRGB[2])*2 -1;
    
    Normal n = Normal(normalRGB[0],normalRGB[1],normalRGB[2]); // Normal from the normal map
    n = Normalize(n); // Make sure it is a unit vector, or else our rotation equations are not valid.
    
    // These are tangent-space normals. We want object-space normals.
    // To solve this, we first find the rotation between (0,0,1) and the normal in the normal map. We  then apply that rotation to the normal from geometry.
    
    // To find the rotation, we find the orthogonal axis to (0,0,1) and the normal. We then solve for the rotation around this orthogonal axis.
    Vector orthogAxis = Cross(Vector(0,0,1),n);
    float rotationAngle = acosf(Dot(Vector(0,0,1),n));

    Transform rotationTransform = Rotate(Degrees(rotationAngle), orthogAxis); // Rotate takes degrees!
    
    *dgBump = dgs;
    dgBump->nn = rotationTransform(dgs.nn);
    
    if (dgs.shape->ReverseOrientation ^ dgs.shape->TransformSwapsHandedness)
        dgBump->nn *= -1.f;

    // Orient shading normal to match geometric normal
    dgBump->nn = Faceforward(dgBump->nn, dgGeom.nn);
    
}


