
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


// film/image.cpp*
#include <sstream>
#include <iostream>
#include <fstream>
#include "stdafx.h"
#include "film/spectralImage.h"
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"
#include "floatfile.h"
#include "cameras/realisticDiffraction.h"
#include <typeinfo>
#include <math.h>

// SpectralImageFilm Method Definitions
SpectralImageFilm::SpectralImageFilm(int xres, int yres, Filter *filt, const float crop[4],
                                     const string &fn, bool openWindow)
: Film(xres, yres, fn) { // fn added by Trisha
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);
    
    // Allocate film image storage
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
    
    // Precompute filter weight table
#define FILTER_TABLE_SIZE 16
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
        float fy = ((float)y + .5f) *
        filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
            float fx = ((float)x + .5f) *
            filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }
    
    // Possibly open window for image display
    if (openWindow || PbrtOptions.openWindow) {
        Warning("Support for opening image display window not available in this build.");
    }
    //Andy added
    nCMRows = nSpectralSamples;
    nCMCols = nSpectralSamples;
}

void SpectralImageFilm::AddSample(const CameraSample &sample,
                                  const Spectrum &L, const Ray &currentRay) {
    
    
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
    int x0 = Ceil2Int (dimageX - filter->xWidth);
    int x1 = Floor2Int(dimageX + filter->xWidth);
    int y0 = Ceil2Int (dimageY - filter->yWidth);
    int y1 = Floor2Int(dimageY + filter->yWidth);
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
        PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(con st_cast<CameraSample *>(&sample));
        return;
    }
    
    //Andy: put in spectrum data here
    float origC[nSpectralSamples];
    L.GetOrigC(origC);

    // Precompute $x$ and $y$ filter table offsets
    int *ifx = ALLOCA(int, x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
        float fx = fabsf((x - dimageX) *
                         filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
    }
    int *ify = ALLOCA(int, y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
        float fy = fabsf((y - dimageY) *
                         filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
    }
    bool syncNeeded = (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
    
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];
            
            // Update pixel values with filtered sample contribution
            Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
            if (!syncNeeded) {
                
                // Loop through all wavelengths
                for (int i = 0; i < nSpectralSamples; i++)
                {
                    pixel.c[i] += filterWt * origC[i]; // Multiply by filter weight
                }
                pixel.weightSum += filterWt;
            }
            else {
                // Safely update _Lxyz_ and _weightSum_ even with concurrency
                // (These atomic add operations have something to do with safe memory allocation when using multiple threads. More info on p.1036 in the PBRT textbook.)
                for (int i = 0; i < nSpectralSamples; i++)
                {
                    AtomicAdd(&pixel.c[i], filterWt * origC[i]);
                }
                AtomicAdd(&pixel.weightSum, filterWt);
            }
            
            // Thi is Andy's old way of doing depth.
            // Why multiply by the filterWt?
            pixel.Z += currentRay.maxt * filterWt;    //add the depth map data member to a pixel - depth is inherently stored in Ray
            
        }
    }

}



//Andy: this needs to be modified to allow for more than 3 channel splatting.  Try to understand this more.  Do we need depth map processing here?
// Trisha: I think this has already been figured out below. The above comment is probably dated.
void SpectralImageFilm::Splat(const CameraSample &sample, const Spectrum &L) {
    if (L.HasNaNs()) {
        Warning("SpectralImageFilm ignoring splatted spectrum with NaN values");
        return;
    }
    
    int x = Floor2Int(sample.imageX), y = Floor2Int(sample.imageY);
    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) return;
    Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
    
    for (int i = 0; i < nSpectralSamples; i++)
    {
        AtomicAdd(&pixel.splatC[i], pixel.splatC[i]);
    }
}


void SpectralImageFilm::GetSampleExtent(int *xstart, int *xend,
                                        int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Floor2Int(xPixelStart + 0.5f + xPixelCount  +
                        filter->xWidth);
    
    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Floor2Int(yPixelStart + 0.5f + yPixelCount +
                        filter->yWidth);
}


void SpectralImageFilm::GetPixelExtent(int *xstart, int *xend,
                                       int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}

//Andy added
// Trisha Note: Andy says this was an early attempt to have matrix multiplies (e.g. for color conversion) directly in PBRT. He says "in retrospect it's kind of useless and tough to use."
// TODO: Clear out conversion matrix stuff.

void SpectralImageFilm::ParseConversionMatrix(string filename){
    string fn = AbsolutePath(ResolveFilename(filename));
    nCMRows = nSpectralSamples;
    nCMCols = nSpectralSamples;
    
    //default wavespecify
    waveSpecify = new float[nCMRows];
    for (int i = 0; i < nCMCols; i++)
        waveSpecify[i] = sampledLambdaStart + (((sampledLambdaEnd - sampledLambdaStart)/nSpectralSamples) * (i - .5));  //TODO: double check int and float issues
    //populate the identity matrix
    float * identity = new float[nCMRows * nCMCols];
    for (int i =0; i < nCMRows * nCMCols; i++)
        identity[i] = 0.f;
    for (int i = 0; i < nCMRows; i++)
        identity[i * nCMCols + i] = 1.f;
    
    
    // Trisha: Get rid of annoying warnings for now. In the future we will remove all this conversion matrix stuff.
    conversionMatrix = identity;
    /*
     vector<float> vals;
     if (!ReadFloatFile(fn.c_str(), &vals)) {
     
     Warning("Unable to read conversion matrix file \"%s\".  Using identity matrix.",
     fn.c_str());
     
     conversionMatrix = identity;
     return;
     }
     // TODO: consider list of wavelengths
     if (vals.size() < 2 || ((vals.size() - 2 - nCMRows) % nSpectralSamples != 0))  //changed to consider the CMRows amount of wave specifications
     {
     Warning("Incorrect conversion matrix file format (wrong number of matrix elements! Using identity matrix.");
     conversionMatrix = identity;
     return;
     }
     nCMRows = vals[0];
     nCMCols = vals[1];
     
     waveSpecify = new float[nCMRows];
     //read wavelength representation data
     for (int i = 2; i < 2  + nCMRows; i++)
     {
     waveSpecify[i-2] = vals[i];
     }
     
     if (nCMRows * nCMCols != vals.size() -2 - nCMRows || nCMCols != nSpectralSamples)
     {
     Warning("Incorrect conversion matrix file format (wrong number of matrix elements! Using identity matrix.");
     conversionMatrix = identity;
     return;
     }
     
     conversionMatrix = new float[nCMCols * nCMRows];
     
     for (int i = 2 + nSpectralSamples; i< vals.size(); i++)
     {
     conversionMatrix[i-2-nSpectralSamples] = vals[i];
     if (debugMode)
     std::cout << "conversionMatrix[" << i-2 << "]=" << conversionMatrix[i-2-nCMRows] << "\t";
     }
     */
    return;
    
}

//Andy changed
void SpectralImageFilm::WriteImage(float splatScale) {
    
    // Convert image to RGB and compute final pixel values
    int nPix = xPixelCount * yPixelCount;
    float *finalC = new float[nSpectralSamples * nPix];
    
    int offset = 0;
    for (int x = 0; x < xPixelCount; ++x) {
        for (int y = 0; y < yPixelCount; ++y) {
            
            for (int ind = 0; ind < nSpectralSamples; ind++)
            {
                finalC[(y * xPixelCount + x)*nSpectralSamples + ind] = (*pixels)(x,y).c[ind];
            }
            
            // TODO: we may need to eliminate this "filtering" later
            float weightSum = (*pixels)(x, y).weightSum;  // Andy: will still need this, but need to make a for loop for all channels
            if (weightSum != 0.f) {
                float invWt = 1.f / weightSum;
                
                for (int i = 0; i < nSpectralSamples; i++)
                {
                    finalC[nSpectralSamples * offset + i] =  max(0.f, finalC[nSpectralSamples*offset + i  ]);
                    
                    // DEBUG
                    //std::cout << "finalC*invWt = " << finalC[nSpectralSamples * offset + i]*invWt << std::endl;
                    
                    // No invWt here.
                    //  Trisha Note: For our simulation, we simply add up all the "rays" that hit the sensor pixel - we don't do any normalization at all. In other words, the higher the number of pixel samples, the larger the pixel values will be. The spectrums given to the lights are all in relative values. In ISET, we scale the output pixel values to a mean luminance to get an actual physical value.
                }
                
                //                // Depth map weighting
                //                finalZ[3 * offset ] = max(0.f, finalZ[3*offset] * invWt);
                //                finalZ[3 * offset + 1] = max(0.f, finalZ[3*offset + 1] * invWt);
                //                finalZ[3 * offset + 2] = max(0.f, finalZ[3*offset + 2] * invWt);
                
            }
            
            
            //Add splat value at pixel
            float * splatC = (*pixels)(x, y).splatC;
            
            for (int i = 0; i < nSpectralSamples; i++)
            {
                
                finalC[nSpectralSamples * offset + i] += splatScale * splatC[nSpectralSamples];
            }
            ++offset;
        }
    }
    
    
    float *finalCMultiplied = new float[nCMRows * nPix];
    //Andy added: multiply by conversion matrix to produce proper output
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            //for each pixel, multiply conversion matrix by existing data
            //these nested for loops are for matrix multiplication
            for (int row = 0; row < nCMRows; row++)
            {
                float tempSum = 0;
                for (int iter = 0; iter < nCMCols; iter++)
                {
                    tempSum += conversionMatrix[row * nCMCols + iter] * finalC[nSpectralSamples * (y * xPixelCount + x) + iter];    //matrix multiplication by column vector
                }
                finalCMultiplied[nCMRows * (x *  yPixelCount + y) + row] = tempSum;  //check this later
            }
            
        }
    }
    
    //Write to file
    std::ofstream myfile;
    int lastPos = imageOutputName.find_last_of(".");
    string newFileName = imageOutputName.substr(0, lastPos) + ".dat";
    myfile.open (newFileName.c_str());
    
    //print out the dimensions of the image
    myfile << xPixelCount << " " << yPixelCount << " " << nCMRows << "\n";
    
    //print out field of view information
    float d = sensorWidth;
    float fieldOfView = 2 * atan(d/(2 * focalLength)) / 3.1415926539 * 180;
    
    myfile << focalLength << " " << fStop << " " << fieldOfView << "\n";
    
    myfile.close();
    
    //open file for binary writing purposes
    FILE * spectralImageBin;
    spectralImageBin = fopen(newFileName.c_str(), "a");
    
    //Write Binary image
    for (int i = 0; i < nCMRows; i++)
    {
        for (int j = 0; j < nPix; j++)
        {
            double r = (double)finalCMultiplied[nCMRows * j + i];
            fwrite((void*)(&r), sizeof(r), 1, spectralImageBin);
        }
    }
    
    fclose(spectralImageBin);
    
    // Release temporary image memory
    delete[] finalC;
    delete[] finalCMultiplied;
    
    // Clear pixels, in case we're doing another render with the same film (e.g. camerasrenderer)
    delete pixels;
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
}


void SpectralImageFilm::UpdateDisplay(int x0, int y0, int x1, int y1,
                                      float splatScale) {
}

//Andy: added this additional construction function to allow for camera pointer - for FOV calculations
SpectralImageFilm *CreateSpectralImageFilm(const ParamSet &params, Filter *filter, Camera * baseCamera)
{
    
    
    // ANDY: THIS WILL REQUIRE SOME RETHINKING FOR MAXIMUM COMPATIBILITY!!!    ... but works for now.  but what if baseCamera is not a compatible type?!
    
    //baseCamera contains information about the camera: information such as focal length is very important to calculate field of view.
    //float focalLength = baseCamera->focalLength
    //also need film dimensions
    //float horizontalFOV = 2 * arctan(d/2f)
    
    //RealisticDiffractionCamera*
    //dynamic_cast<RealisticDiffractionCamera*> (baseCamera);
    
    string myType = typeid(*baseCamera).name();
    std::cout << "typeOfCamera: " << myType << "\n\n";
    
    
    string comparison("26RealisticDiffractionCamera"); //("17PerspectiveCamera"); // ("26RealisticDiffractionCamera");	//**IF COMPARISON IS 0, then strings are THE SAME!!!
    std::cout << "comparison: " << myType.compare(comparison)  << "\n\n";
    
    std::cout << "\n\nin constructor!!\n\n";
    float tfocalLength = 0;
    float tfStop = 0;
    float tSensWidth = 0;
    
    if (myType.compare(comparison) == 0)
    {
        tfocalLength = ((RealisticDiffractionCamera*)baseCamera)->getFocalLength();  // NEED TO ADD THIS STILL
        tfStop = ((RealisticDiffractionCamera*)baseCamera)->getFStop();  //NEED TO ADD THIS STILL
        tSensWidth = ((RealisticDiffractionCamera*)baseCamera)->getSensorWidth();
        
        std::cout << "focalLength: " << tfocalLength;
        std::cout << "\nfStop: " << tfStop << "\n";
    }
    else
    {
        std::cout << "no focal length or fStop information for this camera!\n";
        std::cout << "myType: " << myType << "\n";
    }
    
    
    SpectralImageFilm * newImageFilm = CreateSpectralImageFilm(params, filter);
    newImageFilm->SetFStop(tfStop);
    newImageFilm->SetFocalLength(tfocalLength);
    newImageFilm->SetSensorWidth(tSensWidth);
    return newImageFilm;
}


//Andy: added these for lens information and FOV output to ISET
void SpectralImageFilm::SetFStop(float inputFStop)
{
    fStop = inputFStop;
}
void SpectralImageFilm::SetFocalLength(float inputFocalLength)
{
    focalLength = inputFocalLength;
}
void SpectralImageFilm::SetSensorWidth(float sensWidth)
{
    sensorWidth = sensWidth;
}


//Andy: added conversionmatrix file here
SpectralImageFilm *CreateSpectralImageFilm(const ParamSet &params, Filter *filter) {
    string filename = params.FindOneString("filename", PbrtOptions.imageFile);
    if (filename == "")
#ifdef PBRT_HAS_OPENEXR
        filename = "pbrt.exr";
#else
    filename = "pbrt.tga";
#endif
    
    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);
    
    if (PbrtOptions.quickRender) xres = max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    bool openwin = params.FindOneBool("display", false);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }
    
    
    //Andy changed to allow for parsing of conversionmatrixfile name
    SpectralImageFilm * newFilm = new SpectralImageFilm(xres, yres, filter, crop, filename, openwin);
    
    // Trisha: Comment this out for now since it generates an annoying warning everytime we render. I don't think we use conversion matrix at all right now, so I guess this is a TODO.
    
    string conversionMatrixFilename = params.FindOneString("conversionmatrixfile", "");   //params contains all the parsed attributes of film... we now look for the "conversionmatrixfile" entry, which is a string
    // std::cout<< "\nconversionMatrixFilename: " << conversionMatrixFilename << "\n\n";
    newFilm->ParseConversionMatrix(conversionMatrixFilename);
    
    
    return newFilm;
}


