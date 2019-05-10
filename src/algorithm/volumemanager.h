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

#ifndef VOLUMEMANAGER_H
#define VOLUMEMANAGER_H

#include <QObject> //Need this to use Signal-Slot mechanism

#define OCCLUSION_DOWNSIZE_VOLUME 1 // Downsize volume to compute occlusion spectrum
#define OCCLUSION_FAST_COMPUTE 1 // Use FFT based fast computation of convolution for occlusion
#define TINY 1e-12
#define FILTER_CLASSIFIED_VOLUME 1

struct Histogram {
    int m_nbins;
    float* m_logFreq;
    float *m_freq;
};

#include <itkImage.h>
#include <itkCovariantVector.h>
#include <itkInterpolateImageFilter.h>
#include "lhhistogram.h"
#include "material.h"

typedef itk::Image<float, 3> FloatImageType;
typedef itk::CovariantVector<float, 3> GradientType;
typedef itk::Image<GradientType, 3> GradientImageType;
typedef itk::InterpolateImageFunction<FloatImageType, float> ImageInterpolatorType;
typedef itk::InterpolateImageFunction<GradientImageType, float> GradientInterpolatorType;

class VolumeManager : public QObject
{
    Q_OBJECT

public:
    VolumeManager();
    ~VolumeManager();
    void readNHDR(const char *filename);
    void readVTK(const char* filename);
    int const & width() const { return m_width;}
    int const & height() const { return m_height;}
    int const & depth() const { return m_depth;}
    float const & spacingX() const { return m_spacingX;}
    float const & spacingY() const { return m_spacingY;}
    float const & spacingZ() const { return m_spacingZ;}
    float* data() { return m_data;}
    unsigned char* getCannyEdges() { return m_cannyEdges;}
    float* gradient() { return m_gradient;}
    double* getLHDensityImage() { return m_LHhistogram->getLHDensity(); }
    double* getOcclusionDensityImage();
    itk::Image<float, 3>::Pointer getITKImage();
    Histogram const & histogram() const { return m_histogram; }
    void preprocess(); //Perform preprocessing and data preparation
    void computeMaterialGraph(int LHMinClustSz = -1, int LHMinSamples = -1, bool restorePrevious = false, bool runOnlyAlgo=false);
    vector<Material*>& getMaterials () {return m_materials; }
    vector<QPoint*>& getGraphEdges() {return m_graphEdges; }
    int getBackgroundMaterial() {return m_backgroundMaterial; }
    double* getIntensityRange() {return m_intensityRange; }
    double* getOcclusionRange() {return m_occlusionRange; }
    void saveRGBA(QString filename);
    void saveRGB(QString filename);
    void saveClassifiedVolume(QString filename);
    bool readClassifiedVolume(QString filename);
    void computeRGBAVolume(unsigned char *rgbaVolume, bool storePremultipliedColors = true, int nChannels = 4);
    void computeRGBAVolumeExternal(unsigned char *rgbaVolume);//To be invoked after reading an external material volume

signals:
    void volumeDataCreated(VolumeManager *vm);
    void volumeLHComputed(VolumeManager *vm);
    void volumeEdgesComputed(VolumeManager *vm);
    void volumeGradientComputed(VolumeManager *vm);
    void volumePreprocessCompleted(VolumeManager *vm);
    void materialGraphComputed(VolumeManager *vm);
    void volumeOcclusionComputed(VolumeManager *vm);
    void classifiedVolumeComputed(VolumeManager *vm);
    void classifiedVolumeRead(VolumeManager *vm);

protected slots:
    void computeClassifiedVolume();

private:
    int m_width, m_height, m_depth;
    float m_spacingX, m_spacingY, m_spacingZ;
    float m_min, m_max; // Min, Max of the original data.
    char* m_volumeName;
    char* filePathName;
    unsigned char* m_classifiedVolume;

    Histogram m_histogram;
    LHHistogram *m_LHhistogram;
    vector<Material*> m_materials;
    vector<QPoint*> m_graphEdges; //Store material graph edges in QPoint as (material idx1, material idx 2)
    int m_backgroundMaterial; //Index to background material detected.

    //Derived data
    float *m_data; // Normalized voxel values in the range [0, 1]
    unsigned char *m_cannyEdges;// Edge voxels marked as 255
    float *m_gradient; // Stored as gx, gy, gz, gx, gy, gz, ...
    vtkImageData* m_vtkImageData;
    vtkSmartPointer<vtkImageData> m_imageOcclusionDensityPlot;
    double m_intensityRange[2], m_occlusionRange[2];

    //Private functions
    void computeHistogram();
    void computeCannyEdges(); // Canny edge detection on volume
    void computeGradient();
    void computeLH();
    void computeOcclusionSpectrum();
    void readMaterials();
    void readMaterialEdges();
};

#endif // VOLUMEMANAGER_H
