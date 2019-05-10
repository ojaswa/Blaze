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

#ifndef LHHISTOGRAM_H
#define LHHISTOGRAM_H

#include <vtkSmartPointer.h>
#include <QString>

#define TINY 1e-12
#define LH_CUM_HISTOGRAM_FRACTION 0.9
#define LH_BASIC_THRESHOLDING 0
#define LH_INTERPOLATE_IMAGE_WHILE_RK_TACKING 1
#define LH_USE_CANNY 0
#define LH_RK_TRACKING_DELTA 1
#define LH_GAUSSIAN_SMOOTH_SIGMA 2.0

class vtkImageGradientMagnitude;
class vtkImageData;
class vtkImageInterpolator;
class VolumeManager;

class LHHistogram
{
public:
    LHHistogram(VolumeManager *vm, const char* volumefile);
    ~LHHistogram();
    double* getLHDensity();

    template<class T> void writeLH();
    void compute();
    vtkImageData* getVTKImageData();
private:
    //Data
    char *m_fileName;
    VolumeManager *m_parent;
    vtkSmartPointer<vtkImageData> m_imageData;
    vtkSmartPointer<vtkImageData> m_imageLH; // 32-bit LH for 16-bit volume, 16-bit LH for 8-bit volume. Stores LH pair per voxel: LS bytes store L value, MS byts store H value
    vtkSmartPointer<vtkImageData> m_imageLHDensityPlot;

    //Functions
    void readVolume(const char* filename);
    template<class T> void computeLH();
    double computeEpsilon(vtkImageGradientMagnitude *gmag);
    template<class T> vtkImageData* computeGradientLabelImageLH(vtkImageData *smoothedImage);
    double gradientTrack(vtkImageInterpolator *image, vtkImageInterpolator *gradient, double delta_t, int x, int y, int z);
    void rungeKuttaEstimate(vtkImageInterpolator *gradient, double delta_t, double *pos);
};

// Extern template declarations
extern template void LHHistogram::computeLH<unsigned char>();
extern template void LHHistogram::computeLH<unsigned short>();

extern template void LHHistogram::writeLH<unsigned char>();
extern template void LHHistogram::writeLH<unsigned short>();

#endif // LHHISTOGRAM_H
