/***************************************************************************
**                                                                        **
**  Blaze - a volume rendering and analytics program                      **
**  Copyright (C) 2016-2018 Graphics Research Group, IIIT Delhi           **
**                                                                        **
**  This program is free software: you can redistribute it and/or modify  **
**  it under the terms of the GNU General Public License as published by  **
**  the Free Software Foundation, either version 3 of the License, or     **
**  (at your option) any later version.                                   **
**                                                                        **
**  This program is distributed in the hope that it will be useful,       **
**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
**  GNU General Public License for more details.                          **
**                                                                        **
**  You should have received a copy of the GNU General Public License     **
**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
**                                                                        **
****************************************************************************
**           Author: Ojaswa Sharma                                        **
**           E-mail: ojaswa@iiitd.ac.in                                   **
**           Date  : 14.12.2016                                           **
****************************************************************************/

#include "occlusionspectrum.h"

#include <vtkImageData.h>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <itkImage.h>
#include <itkFFTConvolutionImageFilter.h>
#include <itkVTKImageToImageFilter.h>
#include <itkImportImageFilter.h>

OcclusionSpectrum::OcclusionSpectrum(vtkImageData *_img)
{
    m_image = _img;
    m_occlusion = NULL;

    int dims[3];
    m_image->GetDimensions(dims);
    m_radius = 0.1 * (dims[0] + dims[1] + dims[2]) / 3;//Radius should capture the object sizes in the volume.
}

OcclusionSpectrum::~OcclusionSpectrum()
{
    if(m_occlusion) delete []m_occlusion;
}

double* OcclusionSpectrum::getOcclusionSpectrum()
{
    return m_occlusion;
}

void OcclusionSpectrum::compute()
{
    int* psum = prefixSum();
    unsigned int r2 = m_radius*m_radius;
    int dims[3];
    m_image->GetDimensions(dims);
    long nelems = dims[0]*dims[1]*dims[2];
    m_occlusion = new double[nelems];

    //Loop through all points and create an occlusion volume

    for(int z = 0; z<dims[2]; z++)
        for(int y = 0; y<dims[1]; y++)
            for(int x = 0; x<dims[0]; x++) {
                //Calculate neighborhood extent
                int nextent[6] =
                {
                    std::max(x - m_radius, 0), std::min(x+m_radius, dims[0]-1),
                    std::max(y - m_radius, 0), std::min(y+m_radius, dims[1]-1),
                    std::max(z - m_radius, 0), std::min(z+m_radius, dims[2]-1)
                };

                long sum = 0;
                int num = 0;

                //Loop through each grid point in the projected disk of neighbour sphere
                for(int k=nextent[4]; k<=nextent[5]; k++) //Along z
                    for(int j=nextent[2]; j<=nextent[3]; j++) { // Along y
                        int i = r2 - (y-j)*(y-j) - (z-k)*(z-k);
                        if (i<0) continue; //Outside sphere
                        int t = (k*dims[1] + j)*dims[0]; // Row offset in volume
                        i = (int)sqrt((float)i);//Radius of intersection in X direction
                        int const i0 = std::max(x-i, 0);
                        int const i1 = std::min(x+i, dims[0]-1);

                        //Sum up everything in [i0, i1]
                        sum += psum[t+i1] - (0==i0 ? 0:psum[t+i0-1]);
                        num += i1-i0+1;
                    }
                m_occlusion[x + dims[0]*(y + dims[1]*z)] = num ? ((double)sum/num) : 0;
            }
    delete []psum;
}

//Compute prefix sum along x dimension
int* OcclusionSpectrum::prefixSum()
{
    int dims[3];
    m_image->GetDimensions(dims);
    long nelems = dims[0]*dims[1]*dims[2];

    int* psum = new int[nelems];
    memset(psum, 0, nelems*sizeof(int));

    int* ptr_image = static_cast<int*>(m_image->GetScalarPointer());
    int startIdx, endIdx;
    for(int z=0; z<dims[2]; z++)
        for(int y=0; y<dims[1]; y++)
        {
            startIdx = (y + z*dims[1])*dims[0];
            endIdx = startIdx + dims[0];
            std::partial_sum(ptr_image + startIdx, ptr_image + endIdx, psum + startIdx);
        }
    return psum;
}

// Use Fast Fourier Transform to compute convolution
// This is *much* faster than the other method compute().
void OcclusionSpectrum::computeWithFFT(OcclusionMapType type)
{
    // Create a spherical kernel
    int kernelSize = m_radius*2 + 1;
    typedef itk::Image<double, 3> ImageType;
    typedef itk::FFTConvolutionImageFilter<ImageType> FilterType;
    typedef itk::VTKImageToImageFilter<ImageType> ConvertType;
    typedef itk::ImportImageFilter<double, 3> ImportFilterType;

    // Create kernel data
    ImageType::Pointer kernel = ImageType::New();
    ImageType::SizeType imsize;
    imsize[0] = kernelSize; imsize[1] = kernelSize; imsize[2] = kernelSize;
    ImportFilterType::IndexType start;
    start.Fill(0);
    ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(imsize);
    kernel->SetRegions(region);
    kernel->Allocate();
    kernel->GetPixelContainer()->SetContainerManageMemory(true);
    double *ptr_kernel = reinterpret_cast<double*>(kernel->GetPixelContainer()->GetImportPointer());
    createSphericalKernel(ptr_kernel, kernelSize, type);

    // Create convolution kernel
    FilterType::Pointer convolutionFilter = FilterType::New();
    // Get VTK image into ITK
    ConvertType::Pointer convert = ConvertType::New();
    convert->SetInput(m_image);
    convolutionFilter->SetInput(convert->GetOutput());
    // Set kernel
    convolutionFilter->SetKernelImage(kernel);
    convolutionFilter->Update();

    // Get data from convolution filter
    convolutionFilter->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);
    m_occlusion = reinterpret_cast<double*>(convolutionFilter->GetOutput()->GetPixelContainer()->GetImportPointer());
    convolutionFilter->Update();      
}

void OcclusionSpectrum::createSphericalKernel(double *ptr_kernel, int kernelSize, OcclusionMapType mapType)
{
    int kernelElems = 0;
    int r2 = m_radius*m_radius;
    int ii, jj, kk;
    int dist2, pred;
    for(int k=0; k<kernelSize; k++)
        for(int j=0; j<kernelSize; j++)
            for(int i=0; i<kernelSize; i++) {
                ii = (i - m_radius);
                jj = (j - m_radius);
                kk = (k - m_radius);
                dist2 = ii*ii + jj*jj  + kk*kk;
                pred = r2 - dist2;
                if (pred < 0)
                    ptr_kernel[i + kernelSize*(j + k*kernelSize)] = 0.0; // Outside sphere
                else {
                    ptr_kernel[i + kernelSize*(j + k*kernelSize)] = (mapType == OcclusionMapExponential)? exp(-dist2) : 1.0; //Inside Sphere
                    kernelElems++;
                }
            }
    int kernelSize3 = kernelSize*kernelSize*kernelSize;
    for(int i=0; i<kernelSize3; i++) ptr_kernel[i] /= (double)kernelElems;
}
