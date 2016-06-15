//
//  water.cpp
//  pbrt
//
//  Created by Trisha Lian on 5/26/16.
//
//

#include "stdafx.h"
#include "water.h"
#include "paramset.h"




// WaterVolumeDensity Method Definitions
WaterVolumeDensity *CreateWaterVolumeDensity(const Transform &volume2world,
                                                               const ParamSet &params) {

    float spc = params.FindOneFloat("smallPartConc", 0.);
    float lpc = params.FindOneFloat("largePartConc", 0.);
    float cc = params.FindOneFloat("chlorophyllConc", 0.);
    float domc = params.FindOneFloat("domConc", 0.);
    
    Spectrum plankton = params.FindOneSpectrum("phytoplanktonCurve", 0.);
    
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
    
    return new WaterVolumeDensity(spc,lpc,cc,domc,plankton,BBox(p0, p1),
                                        volume2world);
}