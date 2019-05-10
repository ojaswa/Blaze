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

#include "volumemanager.h"
#include "defines.h"

#include <fstream>
#include <string>
#include <QFileInfo>
#include <QDir>
#include <QFuture>
#include <QFutureWatcher>
#include <QtConcurrent>
#include <QColor>
#include <cstdlib>
#include <math.h>
#include <float.h>

//ITK includes
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkImportImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImageToHistogramFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkContinuousIndex.h>
#include <vtkImageCast.h>
#include <itkMedianImageFilter.h>
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "lhhistogram.h"
#include "occlusionspectrum.h"
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkImageResize.h>
#include <vtkImageCast.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkNrrdReader.h>

using namespace std;

#define EXTERNAL_RGBA_RENDER_ALPHA 0.1

VolumeManager::VolumeManager()
{
    m_data = NULL;
    m_width = m_height = m_depth = 0;
    m_min = m_max = 0.0;
    m_volumeName = new char[256];
    m_histogram.m_nbins = 256;
    m_histogram.m_logFreq = new float[m_histogram.m_nbins];
    m_histogram.m_freq = new float[m_histogram.m_nbins];
    m_LHhistogram = NULL;
    m_intensityRange[0] = -1; m_intensityRange[1] = -1;
    m_occlusionRange[0] = -1; m_occlusionRange[1] = -1;

    if (m_materials.size() > 0) {
        for (int i=0; i<m_materials.size(); i++)
            delete m_materials.at(i);
        m_materials.clear();
    }

    if(m_graphEdges.size() > 0) {
        for (int i=0; i<m_graphEdges.size(); i++)
            delete m_graphEdges.at(i);
        m_graphEdges.clear();
    }

    m_cannyEdges = NULL;
    m_gradient = NULL;
    m_classifiedVolume = NULL;
    //m_LH = NULL;
}

VolumeManager::~VolumeManager()
{
    delete []m_volumeName;
    if(m_data) delete []m_data;
    if(m_histogram.m_freq) delete []m_histogram.m_freq;
    if(m_histogram.m_logFreq) delete []m_histogram.m_logFreq;
    m_histogram.m_nbins = 0;
    if(m_cannyEdges) delete []m_cannyEdges;
    if(m_gradient) delete []m_gradient;
    if(m_LHhistogram) delete m_LHhistogram;
    if(m_classifiedVolume) delete []m_classifiedVolume;
}

double* VolumeManager::getOcclusionDensityImage()
{
    double *densityPtr = static_cast<double*>(m_imageOcclusionDensityPlot->GetScalarPointer());
    return densityPtr;
}

void VolumeManager::readNHDR(const char *filename)
{
    m_LHhistogram = new LHHistogram(this, filename); // Initialize instance of LH Histogram
    m_vtkImageData = m_LHhistogram->getVTKImageData(); //Get pointer to VTK data for other (vtk-based) functions to make use of it.

    ofstream pathNameFile;
    pathNameFile.open("./process/path.txt");

    //Read the .nhdr first in text mode
    string line;
    int linecount = 0;
    ifstream fid(filename);
    int vol_typeSize;
    char datafilename[256];
    if(fid) {
        pathNameFile<<filename<<"\n";
        while (getline(fid, line)) {
            if (strncmp(line.c_str(), "encoding", 8) == 0) {
                //Has to be raw
                if(line.find("raw") == string::npos) {
                    fprintf(stderr, "NHDR file not in RAW format!");
                    fid.close();
                    return;
                }
            } else if (strncmp(line.c_str(), "content", 7) == 0) {
                //Name of the volume
                sscanf(line.c_str(), "content: %s", m_volumeName);
            } else if (strncmp(line.c_str(), "type", 4) == 0) {
                char type_str[256];
                sscanf(line.c_str(), "type: %[^\n\r]\n", type_str);
                if (strncmp(type_str, "unsigned char", 13) == 0){
                    vol_typeSize = 1;
                    pathNameFile<<vol_typeSize;
                }
                else if (strncmp(type_str, "unsigned short", 14) == 0){
                    vol_typeSize = 2;
                    pathNameFile<<vol_typeSize;
                }
                else {
                    fprintf(stderr, "Unknown data type: %s\n", type_str);
                    pathNameFile.close();
                    fid.close();
                    return;
                }
            } else if (strncmp(line.c_str(), "sizes", 5) == 0) {
                sscanf(line.c_str(), "sizes: %d %d %d\n", &m_width, &m_height, &m_depth);
            } else if (strncmp(line.c_str(), "spacings", 8) == 0) {
                sscanf(line.c_str(), "spacings: %f %f %f\n", &m_spacingX, &m_spacingY, &m_spacingZ);
            } else if (strncmp(line.c_str(), "dimension", 9) == 0) {
                int dim;
                sscanf(line.c_str(), "dimension: %d\n", &dim);
                if(dim != 3) {
                    fprintf(stderr, "Not a 3D data!");
                    fid.close();
                    return;
                }
            } else if (strncmp(line.c_str(), "data file", 9) == 0) {
                sscanf(line.c_str(), "data file: %[^\n\r]\n", datafilename);
            }
            linecount++;
        }
    }
    pathNameFile.close();
    fid.close();

    //Load data from binary raw file
    int nelements = m_width*m_height*m_depth;
    m_data  = new float[nelements];

    QFileInfo nhdr_file(filename);
    QFileInfo raw_file(datafilename);
    QString qdatafile("");
    qdatafile.append(nhdr_file.path()).append("/").append(raw_file.fileName());
    FILE *data_fid = fopen(qdatafile.toStdString().c_str(), "rb");
    if(data_fid) {
#define IO_BLOCK_SIZE 4096
        char *block4k = new char[IO_BLOCK_SIZE*vol_typeSize];
        int elements_read;
        int k=0;
        do {
            elements_read = fread((void*)block4k, vol_typeSize, IO_BLOCK_SIZE, data_fid);
            if(elements_read < 0) break;
            //Convert from NHRD data type to float and store in volume
            if(vol_typeSize == 1)
                for(int i=0; i<elements_read; i++) *(m_data + k++) = (float)*((unsigned char*)block4k + i);
            else if(vol_typeSize == 2)
                for(int i=0; i<elements_read; i++) *(m_data + k++) = (float)*((unsigned short*)block4k + i);
        } while(elements_read == IO_BLOCK_SIZE);

        delete []block4k;
        fclose(data_fid);
    }

    //Find min-max
    m_min = m_data[0];
    m_max = m_data[0];
    for(int i=1; i<nelements; i++) {
        if (m_data[i] < m_min) m_min = m_data[i];
        if (m_data[i] > m_max) m_max = m_data[i];
    }

    //Rescale data to [0, 1]
    for(int i=0; i<nelements; i++) {
        m_data[i] = (m_data[i] - m_min) / (m_max - m_min);
    }

    fprintf(stderr, "Read volume: %s\n", filename);
    fprintf(stderr, "\tName: %s\n", m_volumeName);
    fprintf(stderr, "\tType: %s\n", (vol_typeSize==1)?"unsigned char":(vol_typeSize == 2)?"unsigned short":"unknown");
    fprintf(stderr, "\tSize: %d x %d x %d\n", m_width, m_height, m_depth);
    fprintf(stderr, "\tSpacing: %f x %f x %f\n", m_spacingX, m_spacingY, m_spacingZ);
    fprintf(stderr, "\tData range: [%f, %f] normalized to [0, 1]\n", m_min, m_max);

    computeHistogram();

    emit volumeDataCreated(this);
}

void VolumeManager::preprocess()
{
    //Add any volume preprocessing code here.
    fprintf(stderr, "Processing volume:\n");

#if TIME_PROCESSES
    computeCannyEdges();// Compute Canny before LH, since LH may use Canny edges.
    computeGradient();
    computeLH();
    computeOcclusionSpectrum();
#else
    //Run independent tasks in seperate threads with QFuture class
    //N.B.: Any task that depends on another must be started after the previous one has finished (use waitForFinished())

    QFuture<void> futureComputeCannyEdges = QtConcurrent::run(this, &VolumeManager::computeCannyEdges);
    QFuture<void> futureComputeGradient = QtConcurrent::run(this, &VolumeManager::computeGradient);

    // Finish up Canny computation before LH, since LH may use Canny edges.
    futureComputeCannyEdges.waitForFinished();
    // Finish up Gradient computation before LH, since LH may use Gradient edges.
    futureComputeGradient.waitForFinished();

    QFuture<void> futureComputeLH = QtConcurrent::run(this, &VolumeManager::computeLH);
    QFuture<void> futureComputeOcclusion = QtConcurrent::run(this, &VolumeManager::computeOcclusionSpectrum);

    // Sync-up tasks here.
    futureComputeLH.waitForFinished();
    futureComputeOcclusion.waitForFinished();
#endif

    emit volumePreprocessCompleted(this);
    fprintf(stderr, "Done.\n");
}

void VolumeManager::readVTK(const char *filename)
{
    //TODO:
}

void VolumeManager::computeLH() {
    fprintf(stderr, "\tComputing LH histogram... \n");

    m_LHhistogram->compute();

    //Signal task completion to the application
    emit volumeLHComputed(this);
}

void VolumeManager::computeMaterialGraph(int LHMinClustSz, int LHMinSamples, bool restorePrevious, bool runOnlyAlgo){
    fprintf(stderr, "Computing material graph from edge information...\n");

    char cmd[512];
    sprintf(cmd,"python3 ./scripts/LHDataAnalysis.py %d %d %d %d", LHMinClustSz, LHMinSamples, restorePrevious?1:0, runOnlyAlgo?1:0);
    fprintf(stderr, "Command: %s\n", cmd);
    int signal=system(cmd);
    fprintf(stderr, "Material computation exited with signal %d.\n",signal);

    readMaterials();
    computeClassifiedVolume();
    readMaterialEdges();

    //Signal task completion to the application
    emit materialGraphComputed(this);
}

void VolumeManager::readMaterials()
{
    //Clear existing materials, if any.
    if(m_materials.size() > 0) {
        for (int i=0; i<m_materials.size(); i++)
            delete m_materials.at(i);
    }
    m_materials.clear();

    // Read material_reduce file
    ifstream mat_file;
    char mat_filename[512] = "./process/materials_reduced.txt";
    mat_file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        mat_file.open(mat_filename);
    }
    catch (std::ifstream::failure e) {
       fprintf(stderr, "Exception openinig file: %s\n\t%s\n", mat_filename, e.what());
       return;
    }

    int num_mat;

    mat_file>>num_mat;
    unsigned char rgba[4];
    double min_occlussion;
    double max_occlussion;
    m_backgroundMaterial = -1;
    double min_muval=numeric_limits<double>::infinity();
    double mu, sigma,bkg_val;
    vector<string> tokens;
    string token;

    for(int i=0;i<num_mat;i++) {
        string temp;
        mat_file>>temp;        
        stringstream ss(temp);
        tokens.clear();

        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        mu = stod(tokens.at(0));
        sigma = stod(tokens.at(1));

        rgba[0]=(unsigned char)stoi(tokens.at(2));
        rgba[1]=(unsigned char)stoi(tokens.at(3));
        rgba[2]=(unsigned char)stoi(tokens.at(4));
        rgba[3]=(unsigned char)stoi(tokens.at(5));
        bkg_val=stod(tokens.at(6));
        if(bkg_val<min_muval){
            min_muval=bkg_val;
            m_backgroundMaterial=i;
        }
        Material *m = new Material();
        m->m_mean = mu;
        m->m_sigma = sigma;
        m->m_color.r = rgba[0]; m->m_color.g = rgba[1]; m->m_color.b = rgba[2]; m->m_color.a = rgba[3];
        m->m_visible = true;
        m->m_isBackgroud = false;
        m_materials.push_back(m);
    }
    fprintf(stderr,"Background material idx %d at val %lf\n",m_backgroundMaterial,min_muval);
    (m_materials[m_backgroundMaterial])->m_color.a = 0; //Explicitely set background alpha to zero.
    (m_materials[m_backgroundMaterial])->m_visible = false; //Set it  to invisible
    (m_materials[m_backgroundMaterial])->m_isBackgroud = true; //Set corresponding flag

    mat_file>>min_occlussion;
    mat_file>>max_occlussion;
    mat_file.close();
}

void VolumeManager::readMaterialEdges()
{
    //Read material graph file
    m_graphEdges.clear();
    char graph_filename[512] = "./process/edges_reduced.txt";
    ifstream graph_file;
    graph_file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        graph_file.open(graph_filename);
    }
    catch (std::ifstream::failure e) {
       fprintf(stderr, "Exception openinig file: %s\n\t%s\n", graph_filename, e.what());
       return;
    }

    vector<string> tokens;
    string token;
    int num_mat = m_materials.size();
    int num_nodes, num_edges;
    {
        string temp;
        graph_file>>temp;
        stringstream ss(temp);
        tokens.clear();

        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
    }
    num_nodes = stoi(tokens.at(0));
    num_edges = stoi(tokens.at(1));
    if(num_nodes != num_mat) {
        fprintf(stderr, "\nMaterial mismatch detected after graph computation!\n\n");
    }

    for(int i=0; i<num_edges; i++) {
        string temp;
        graph_file>>temp;
        stringstream ss(temp);
        tokens.clear();

        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        QPoint *e = new QPoint(stoi(tokens.at(0)), stoi(tokens.at(1)));
        m_graphEdges.push_back(e);
    }
    graph_file.close();
}

void VolumeManager::computeClassifiedVolume()
{
    fprintf(stderr, "Computing classified volume... \n");

#if TIME_PROCESSES
    QElapsedTimer timer;
    timer.start();
#endif
    // Interpolating all occlusion values to alpha values
    // create rgba volume
    int dim[3];
    m_vtkImageData->GetDimensions(dim);
    if(m_classifiedVolume == NULL)
        m_classifiedVolume= new unsigned char[dim[0]*dim[1]*dim[2]];
    vtkSmartPointer<vtkImageCast> ic = vtkSmartPointer<vtkImageCast>::New();
    ic->SetInputData(m_vtkImageData);
    ic->SetOutputScalarTypeToDouble();
    ic->Update();
    double *intensityPtr = static_cast<double*>(ic->GetOutput()->GetScalarPointer());

    // Adding statistical ML Classifier for material assignment
    int points_allocated[m_materials.size()];
    for(int i=0;i<m_materials.size();i++){points_allocated[i]=0;}
    int num_mat = m_materials.size();
    int flag, curridx;
    double mu, sigma;
    bool all_zero;
    double ival;
    double max_likelihood;
    double likelihood;
    int ival_idx;
    for(int z=0;z<dim[2];z++){
        for(int y=0;y<dim[1];y++){
            for(int x=0;x<dim[0];x++){
                flag=0;
                curridx=x+(y+z*dim[1])*dim[0];
                ival = intensityPtr[curridx];
                max_likelihood=-1;
                likelihood = 0.0;
                ival_idx=-1;
                all_zero=true;
                for(int i=0;i<num_mat;i++){
                    mu = m_materials[i]->m_mean;
                    sigma = m_materials[i]->m_sigma;
                    if(sigma==0.0){
                        sigma=1e-6;
                    }
                    likelihood= (exp(-1*(pow(ival-mu,2))/(2*pow(sigma,2))))*(1/(sqrt((2*M_PI))*sigma));			
                    if(likelihood>max_likelihood){
                        max_likelihood=likelihood;
                        ival_idx=i;
                        flag=1;
                    }
                    if(max_likelihood>1e-8){
                        all_zero=false;
                    }
                    // flag=1;
                }
                if(all_zero){
                    max_likelihood=numeric_limits<double>::max();
                    for(int i=0;i<num_mat;i++){
                        mu=m_materials[i]->m_mean;
                        if(fabs(mu-ival)<max_likelihood){
                            max_likelihood=fabs(mu-ival);
                            ival_idx=i;
                        }
                    }
                }
                // if(max_likelihood<=numeric_limits<double>::min()){flag=0;}
                m_classifiedVolume[curridx] = (flag==0)?0:(ival_idx+1);// 0: invalid/unclassified class ID, valid IDs start from 1
                points_allocated[ival_idx]+=1;
            }
        }
    }

#if FILTER_CLASSIFIED_VOLUME
    typedef itk::Image<unsigned char, 3> ImageType;
    //Create ITK image around raw pointer
    typedef itk::ImportImageFilter<unsigned char, 3>  ImportFilterType;
    ImportFilterType::Pointer importFilter = ImportFilterType::New();

    ImageType::SizeType imsize;
    imsize[0] = m_width;
    imsize[1] = m_height;
    imsize[2] = m_depth;

    ImportFilterType::IndexType start;
    start.Fill(0);
    ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(imsize);
    importFilter->SetRegion(region);

    itk::SpacePrecisionType origin[3] = {0.0, 0.0, 0.0};
    importFilter->SetOrigin(origin);

    itk::SpacePrecisionType spacing[3];
    spacing[0] = m_spacingX;
    spacing[1] = m_spacingY;
    spacing[2] = m_spacingZ;
    importFilter->SetSpacing(spacing);

    const unsigned int numberOfVoxels = imsize[0] * imsize[1] * imsize[2];
    importFilter->SetImportPointer(m_classifiedVolume, numberOfVoxels, false);
    itk::Image<unsigned char, 3>::Pointer classifiedVolume = importFilter->GetOutput();
    classifiedVolume->Update();

    //Perform median filtering to produce spatially coherent classification
    typedef itk::MedianImageFilter<ImageType, ImageType> FilterType;
    FilterType::InputSizeType radius;
    radius.Fill(1);

    FilterType::Pointer medianFilter = FilterType::New();
    medianFilter->SetRadius(radius);
    medianFilter->SetInput(classifiedVolume);
    medianFilter->Update();

    unsigned char *medianVol = static_cast<unsigned char*>(medianFilter->GetOutput()->GetPixelContainer()->GetImportPointer());
    memcpy(m_classifiedVolume,medianVol, numberOfVoxels*sizeof(unsigned char)); //Copy back. Required? Perhaps yes.

    //Update statistics
    for(int i=0;i<m_materials.size();i++){points_allocated[i]=0;}
    for(long i=0; i<numberOfVoxels; i++)
        points_allocated[m_classifiedVolume[i]-1]++;
#endif

#if TIME_PROCESSES
    int nMilliseconds = timer.elapsed();
    fprintf(stderr, "Volume classification took %d ms.\n", nMilliseconds);
#endif

    fprintf(stderr, "Points allocation distribution in classified volume: \n");
    for(int i=0;i<m_materials.size();i++){
        fprintf(stderr, "Material %d: mean = %04.4lf, composition = %02.2f %%\n",i,m_materials[i]->m_mean,(float)points_allocated[i]*100.0/(dim[0]*dim[1]*dim[2]));
    }
    fprintf(stderr, "done.\n");
    emit classifiedVolumeComputed(this);
    
}

void VolumeManager::computeRGBAVolume(unsigned char *rgbaVolume, bool storePremultipliedColors, int nChannels)
{
    //Create a table of (premultiplied by default) colors for RGBA assignment
    unsigned char *colors = new unsigned char[4*(1+m_materials.size())];
    float alpha = 0.0;
    colors[0] = 0;
    colors[1] = 0;
    colors[2] = 0;
    colors[3] = 0;
    if(storePremultipliedColors) {
        for(int i=1; i<=m_materials.size(); i++) {    
            if(!m_materials[i-1]->m_visible) {
                alpha = 0.0;
                colors[i*4+3] = 0;
            } else {
                alpha = float(m_materials[i-1]->m_color.a)/255.0;
                colors[i*4+3] = m_materials[i-1]->m_color.a;
            }
            colors[i*4] = (unsigned char)(m_materials[i-1]->m_color.r*alpha);
            colors[i*4+1] = (unsigned char)(m_materials[i-1]->m_color.g*alpha);
            colors[i*4+2] = (unsigned char)(m_materials[i-1]->m_color.b*alpha);
        }
    } else {
        for(int i=1; i<=m_materials.size(); i++) {
            colors[i*4] = (unsigned char)(m_materials[i-1]->m_color.r);
            colors[i*4+1] = (unsigned char)(m_materials[i-1]->m_color.g);
            colors[i*4+2] = (unsigned char)(m_materials[i-1]->m_color.b);
            colors[i*4+3] = m_materials[i-1]->m_color.a;
        }
    }

    //Assign colors to the entire volume
    if (nChannels != 3) nChannels = 4; // RGB or RGBA volume? RGB only may be needed for save to disk.
    if(storePremultipliedColors) nChannels = 4; // Always use RGBA for OpenGL display.
    int flag, curridx, i;
    int classID = 0;
    for(int z=0;z<m_depth;z++){
        for(int y=0;y<m_height;y++){
            for(int x=0;x<m_width;x++){
                flag=0;
                curridx=x+(y+z*m_height)*m_width;
                classID = (m_classifiedVolume[curridx])<<2;
                if (nChannels == 3) curridx *=3;
                else curridx <<=2;
                rgbaVolume[curridx] = colors[classID];
                rgbaVolume[curridx+1] = colors[classID+1];
                rgbaVolume[curridx+2] = colors[classID+2];
                if(nChannels == 4) rgbaVolume[curridx+3] = colors[classID+3];
			}
		}
	}
    delete []colors;
}	

void VolumeManager::computeHistogram()
{
    int count = m_width*m_height*m_depth;
    for(int i=0; i<m_histogram.m_nbins; i++)
        m_histogram.m_freq[i] = 0.0;

    for(int i=0; i<count; i++)
        m_histogram.m_freq[(int)ceil(m_data[i]*m_histogram.m_nbins)]++;

    for(int i=0; i<m_histogram.m_nbins; i++)
        m_histogram.m_logFreq[i] = log(1.0 + m_histogram.m_freq[i]);
}

void VolumeManager::computeCannyEdges()
{
    fprintf(stderr, "\tDetecting Canny edges... \n");

    typedef itk::Image<float, 3> InputImageType;
    typedef itk::Image<unsigned char, 3> OutputImageType;

    // Just some good parameters for Canny. Change if required.
    float variance = 1.0;
    float lowerThreshold = 0.05;
    float upperThreshold = 0.1;

    typedef itk::CannyEdgeDetectionImageFilter<InputImageType, InputImageType> FilterType;
    FilterType::Pointer cannyFilter = FilterType::New();
    cannyFilter->SetInput(getITKImage());
    cannyFilter->SetVariance(variance);
    cannyFilter->SetLowerThreshold(lowerThreshold);
    cannyFilter->SetUpperThreshold(upperThreshold);

    typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> RescaleType;
    RescaleType::Pointer rescalar = RescaleType::New();
    rescalar->SetInput(cannyFilter->GetOutput());
    rescalar->SetOutputMinimum(0);
    rescalar->SetOutputMaximum(255);

    //Erode Canny output
    typedef itk::BinaryBallStructuringElement<OutputImageType::PixelType,3> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(1);
    structuringElement.CreateStructuringElement();
    typedef itk::BinaryDilateImageFilter <OutputImageType, OutputImageType, StructuringElementType> BinaryDilateImageFilterType;
    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(rescalar->GetOutput());
    dilateFilter->SetKernel(structuringElement);

    OutputImageType *output = dilateFilter->GetOutput();
    output->Update();
    OutputImageType::PixelContainer *container;
    container = output->GetPixelContainer();
    container->SetContainerManageMemory(false);
    m_cannyEdges = (unsigned char*) container->GetImportPointer();

    //Signal task completion to the application
    emit volumeEdgesComputed(this);
}

void  VolumeManager::computeGradient() {
    fprintf(stderr, "\tComputing gradient... \n");

    typedef itk::CovariantVector<float, 3> GradientType;
    typedef itk::Image<float, 3> InputImageType;
    typedef itk::Image<GradientType, 3> OutputImageType;

    using FilterType = itk::SmoothingRecursiveGaussianImageFilter< InputImageType, InputImageType >;
    FilterType::Pointer smoothFilter = FilterType::New();
    smoothFilter->SetSigma(2.0);
    smoothFilter->SetInput(getITKImage());

    typedef itk::GradientImageFilter<InputImageType, float, float, OutputImageType> GradientFilterType;
    GradientFilterType::Pointer gradientFilter  = GradientFilterType::New();
    gradientFilter->SetInput(smoothFilter->GetOutput());
    gradientFilter->Update();

    gradientFilter->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);
    m_gradient = reinterpret_cast<float*>(gradientFilter->GetOutput()->GetPixelContainer()->GetImportPointer());
    gradientFilter->Update();

    //Signal task completion to the application
    emit volumeGradientComputed(this);
}

itk::Image<float, 3>::Pointer VolumeManager::getITKImage()
{
    typedef itk::Image<float, 3> ImageType;
    typedef itk::ImportImageFilter<float, 3>  ImportFilterType;
    ImportFilterType::Pointer importFilter = ImportFilterType::New();

    ImageType::SizeType imsize;
    imsize[0] = m_width;
    imsize[1] = m_height;
    imsize[2] = m_depth;

    ImportFilterType::IndexType start;
    start.Fill(0);
    ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(imsize);
    importFilter->SetRegion(region);

    itk::SpacePrecisionType origin[3] = {0.0, 0.0, 0.0};
    importFilter->SetOrigin(origin);

    itk::SpacePrecisionType spacing[3];
    spacing[0] = m_spacingX;
    spacing[1] = m_spacingY;
    spacing[2] = m_spacingZ;

    importFilter->SetSpacing(spacing);

    const unsigned int numberOfVoxels = imsize[0] * imsize[1] * imsize[2];
    importFilter->SetImportPointer(m_data, numberOfVoxels, false);

    itk::Image<float, 3>::Pointer retImg = importFilter->GetOutput();
    retImg->Update();
    return retImg;
}

void VolumeManager::computeOcclusionSpectrum()
{
    fprintf(stderr, "\tComputing occlusion spectrum... \n");

#if TIME_PROCESSES
    QElapsedTimer timer;
    timer.start();
#endif

    int dim[3];
    //Downscale intensity volume. Otherwise computations take a long time.
    //Maximum dimension is set to 100.
#if OCCLUSION_DOWNSIZE_VOLUME
    m_vtkImageData->GetDimensions(dim);
    int maxdim = max(dim[0], max(dim[1], dim[2]));
    double scale = 100.0/double(maxdim);
    // double scale = double(1);
    vtkSmartPointer<vtkImageResize> imresize = vtkSmartPointer<vtkImageResize>::New();
    imresize->SetInputData(m_vtkImageData);
    int newdim[3];
    newdim[0] = ceil(scale*double(dim[0]));
    newdim[1] = ceil(scale*double(dim[1]));
    newdim[2] = ceil(scale*double(dim[2]));
    imresize->SetOutputDimensions(newdim[0], newdim[1], newdim[2]);
    imresize->Update();
#endif

    vtkImageData* intensityImage;
#if OCCLUSION_DOWNSIZE_VOLUME
    intensityImage = imresize->GetOutput();
#else
    intensityImage = m_vtkImageData;
#endif

    //Compute occlusion spectrum. Typecase original image to Int
    //NOTE: This works only with input images of type: CHAR, SHORT, and INT
    vtkSmartPointer<vtkImageCast> ic = vtkSmartPointer<vtkImageCast>::New();
    ic->SetInputData(intensityImage);
#if OCCLUSION_FAST_COMPUTE
    ic->SetOutputScalarTypeToDouble();
#else
    ic->SetOutputScalarTypeToInt();
#endif
    ic->Update();

    OcclusionSpectrum *os = new OcclusionSpectrum(ic->GetOutput());
#if OCCLUSION_FAST_COMPUTE
    os->computeWithFFT(OcclusionMapLinear);
#else
    os->compute();
#endif

#if TIME_PROCESSES
    int nMilliseconds = timer.elapsed();
    fprintf(stderr, "\tOcclusion computation took %d ms.\n", nMilliseconds);
#endif

    // Export occlusion and intensity data for analysis
    int i,j,k;
    double* occlusionData = os->getOcclusionSpectrum();
    ic->SetOutputScalarTypeToDouble();
    ic->Update();
    ic->GetOutput()->GetDimensions(dim);
    double* intensityData = reinterpret_cast<double*>(ic->GetOutput()->GetScalarPointer());
    //Write to disk
    ofstream tempout;
    tempout.open("./process/intensity_data.txt");
    tempout<<dim[0]<<" "<<dim[1]<<" "<<dim[2]<<"\n";
    for(i=0;i<dim[2];i++)
        for(j=0;j<dim[1];j++)
            for(k=0;k<dim[0];k++){
                tempout<<intensityData[k+j*dim[0]+i*dim[1]*dim[0]]<<"\n";
            }
    ofstream occ;
    occ.open("./process/occlusion_data.txt");
    for(i=0;i<dim[2];i++)
        for(j=0;j<dim[1];j++)
            for(k=0;k<dim[0];k++){
                occ<<occlusionData[k+j*dim[0]+i*dim[1]*dim[0]]<<"\n";
            }
    occ.close();
    tempout.close();

    //Create a density image of bivariate histogram (for display purpose)
    //Scale input image to 8-bit
    unsigned char* scaledImageData;
    double valuesRange[2];
    intensityImage->GetScalarRange(valuesRange);
    m_intensityRange[0] = valuesRange[0];
    m_intensityRange[1] = valuesRange[1];
    if(intensityImage->GetScalarType() != VTK_UNSIGNED_CHAR)
    {
        intensityImage->GetDimensions(dim);
        scale = 255.0/(valuesRange[1] - valuesRange[0]);
        vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
        scaler->SetOutputScalarTypeToUnsignedChar();
        scaler->SetScale(scale);
        scaler->SetShift(-valuesRange[0]);
        scaler->SetInputData(intensityImage);
        scaler->Update();
        vtkSmartPointer<vtkImageCast> castFilter = vtkSmartPointer<vtkImageCast>::New();
        castFilter->SetInputConnection(scaler->GetOutputPort());
        castFilter->SetOutputScalarTypeToUnsignedChar();
        castFilter->Update();
        scaledImageData = reinterpret_cast<unsigned char*> (castFilter->GetOutput()->GetScalarPointer());
    } else
        scaledImageData = reinterpret_cast<unsigned char*> (intensityImage->GetScalarPointer());

    //Create density plot
    m_imageOcclusionDensityPlot = vtkSmartPointer<vtkImageData>::New();
    m_imageOcclusionDensityPlot->SetDimensions(256, 256, 1); // Is a 2D image of counts
    m_imageOcclusionDensityPlot->SetSpacing(1, 1, 1);
    m_imageOcclusionDensityPlot->SetOrigin(0, 0, 0);
    m_imageOcclusionDensityPlot->AllocateScalars(VTK_DOUBLE, 1);

    valuesRange[0] = DBL_MAX; valuesRange[1] = DBL_MIN;
    intensityImage->GetDimensions(dim);
    long nelems = dim[0]*dim[1]*dim[2];
    for(int i=0; i<nelems; i++) {
        if(occlusionData[i] < valuesRange[0]) valuesRange[0] = occlusionData[i];
        if(occlusionData[i] > valuesRange[1]) valuesRange[1] = occlusionData[i];
    }
    //fprintf(stderr,"Occlusion range: %f, %f\n", valuesRange[0], valuesRange[1]);
    m_occlusionRange[0] = valuesRange[0];
    m_occlusionRange[1] = valuesRange[1];

    scale = 255.0/(valuesRange[1] - valuesRange[0]);
    double *densityPtr = static_cast<double*>(m_imageOcclusionDensityPlot->GetScalarPointer());
    int I, O;
    long ndens = 256*256;
    double freq;
    for(long i=0; i<ndens; i++)
        *(densityPtr + i) = 0.0; // Initialize with 1 so that taking log is not problematic.
    for(long i=0; i<nelems; i++) {
        I = *(scaledImageData + i);
        O = (int)((occlusionData[i] - valuesRange[0])*scale);
        freq = *(densityPtr + I + O*256);
        *(densityPtr + I + O*256) = freq + 1.0;
    }
    for(long i=0; i<ndens; i++) //Take log for better color contrast
        *(densityPtr + i) = logf(1.0 + *(densityPtr + i));
    
    delete os;

    //Signal task completion to the application    
    emit volumeOcclusionComputed(this);
}

void VolumeManager::saveRGB(QString filename)
{
    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    vtkSmartPointer<vtkImageData> rgbVolume = vtkSmartPointer<vtkImageData>::New();
    rgbVolume->SetDimensions(m_width, m_height, m_depth);
    rgbVolume->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    unsigned char* rgbptr = static_cast<unsigned char*> (rgbVolume->GetScalarPointer());
    unsigned char* voxels = new unsigned char[3*m_width*m_height*m_depth];
    computeRGBAVolume(voxels, false, 3);
    memcpy(rgbptr, voxels, m_width*m_height*m_depth*3);
    delete []voxels;

    //Save RGB volume
    vtkSmartPointer<vtkImageToStructuredPoints> image2points = vtkSmartPointer<vtkImageToStructuredPoints>::New();
    image2points->SetInputData(rgbVolume);
    image2points->Update();

    writer->SetInputData(image2points->GetOutput());
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename.toStdString().c_str());
    writer->Update();
}

void VolumeManager::saveRGBA(QString filename)
{
    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    vtkSmartPointer<vtkImageData> rgbaVolume = vtkSmartPointer<vtkImageData>::New();
    rgbaVolume->SetDimensions(m_width, m_height, m_depth);
    rgbaVolume->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
    unsigned char* rgbaptr = static_cast<unsigned char*> (rgbaVolume->GetScalarPointer());
    unsigned char* voxels = new unsigned char[4*m_width*m_height*m_depth];
    computeRGBAVolume(voxels, false);
    memcpy(rgbaptr, voxels, m_width*m_height*m_depth*4);
    delete []voxels;

    //Save
    vtkSmartPointer<vtkImageToStructuredPoints> image2points = vtkSmartPointer<vtkImageToStructuredPoints>::New();
    image2points->SetInputData(rgbaVolume);
    image2points->Update();

    writer->SetInputData(image2points->GetOutput());
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename.toStdString().c_str());
    writer->Update();
}

void VolumeManager::saveClassifiedVolume(QString filename)
{
    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    vtkSmartPointer<vtkImageData> classVolume = vtkSmartPointer<vtkImageData>::New();
    classVolume->SetDimensions(m_width, m_height, m_depth);
    classVolume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
    unsigned char* classptr = static_cast<unsigned char*> (classVolume->GetScalarPointer());
    memcpy(classptr, m_classifiedVolume, m_width*m_height*m_depth);
    long nelem = m_width*m_height*m_depth;
    for(long i=0; i<nelem; i++) //Merge the background class with invalid class ID 0. I.e. set background ID to 0
        if(classptr[i] == (m_backgroundMaterial+1)) classptr[i] = 0;

    //Save RGB volume
    vtkSmartPointer<vtkImageToStructuredPoints> image2points = vtkSmartPointer<vtkImageToStructuredPoints>::New();
    image2points->SetInputData(classVolume);
    image2points->Update();

    writer->SetInputData(image2points->GetOutput());
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename.toStdString().c_str());
    writer->Update();
}

bool VolumeManager::readClassifiedVolume(QString filename)
{
    QFileInfo fi(filename);
    QString filetype = fi.suffix();
    vtkSmartPointer<vtkImageData> classifiedVolume;

    if (filetype.compare("nhdr")==0 | filetype.compare("nrrd")==0) //Load NRRD file
    {
        vtkSmartPointer<vtkNrrdReader> reader = vtkSmartPointer<vtkNrrdReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->Update();
        classifiedVolume = reader->GetOutput();
    }
    else if (filetype.compare("vtk") == 0) // Load VTK file
    {
        vtkSmartPointer<vtkStructuredPointsReader> reader  = vtkSmartPointer<vtkStructuredPointsReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->Update();
        classifiedVolume = reader->GetOutput();
    }

    //Check if the input file dimensions match the in-memory volume data
    int dim_vol[3];
    m_vtkImageData->GetDimensions(dim_vol);
    int dim_cls[3];
    classifiedVolume->GetDimensions(dim_cls);
    if((dim_vol[0] != dim_cls[0]) ||(dim_vol[1] != dim_cls[1]) || (dim_vol[2] != dim_cls[2])) {
        fprintf(stderr, "Voxel size of classified volume %s does not match with that of in-memory grayscale volume. Aborting classified volume read...\n", filename.toStdString().c_str());
        return false;
    }

    //Copy
    bool overwritten = true;
    long nelem = dim_vol[0]*dim_vol[1]*dim_vol[2];
    if(m_classifiedVolume == NULL) {
        m_classifiedVolume= new unsigned char[nelem];
        overwritten = false;
    }
    unsigned char* classptr = static_cast<unsigned char*> (classifiedVolume->GetScalarPointer());
    memcpy(m_classifiedVolume, classptr, sizeof(unsigned char)*nelem);
    fprintf(stderr, "Input classified volume %s from file %s.\n", overwritten?"overwritten":"read", filename.toStdString().c_str());

    emit classifiedVolumeRead(this);
    return true;
}

void VolumeManager::computeRGBAVolumeExternal(unsigned char *rgbaVolume)
{
    //Calculate well seperated colors for labels
    int dim[3];
    m_vtkImageData->GetDimensions(dim);
    long nelem = dim[0]*dim[1]*dim[2];
    unsigned char hist[256];
    for(int i=0; i<256; i++) hist[i] = 0;
    for(long i=0; i<nelem; i++) //Count number of unique materials in the volume
        hist[m_classifiedVolume[i]]++;

    //Count total number of materials
    unsigned char n_materials = 0;
    for(int i=0; i<256; i++)
        if(hist[i] > 0) n_materials++;
    fprintf(stderr, "A total of %d class labels exist in the material volume.\n", n_materials);

    // Create well seperated colors and saved as premultiplied
    unsigned char n_colors = n_materials;
    unsigned char *colors = new unsigned char[4*256];
    QColor color;
    float h, s, v, a;
    int start=0;
    if (hist[0] > 0) {//There are voxels corresponding to background. Assign white with full transparency
        fprintf(stderr, "Label 0 is considered as transparent background.\n");
        n_colors = n_materials - 1;
        start = 1;
        colors[0] = 0;
        colors[1] = 0;
        colors[2] = 0;
        colors[3] = 0;
    }
    int j=0;
    for(int i=start; i<256; i++) {//Assign colors for existing labels
        if(hist[i] == 0) continue;
        h = (float)(j++)/(float)(n_colors-1); // Uniformly sample Hue space
        s = 1.0;
        v = 1.0;
        a = EXTERNAL_RGBA_RENDER_ALPHA;
        color.setHsvF(h, s, v, a);//Save HSV and convert to RGB below
        colors[4*i+0] = (unsigned char)(color.red()*a);
        colors[4*i+1] = (unsigned char)(color.green()*a);
        colors[4*i+2] = (unsigned char)(color.blue()*a);
        colors[4*i+3] = (unsigned char)color.alpha();
    }

    //Assign colors to each voxel
    int mat;
    for (long i=0; i<nelem; i++) {
        mat = m_classifiedVolume[i];
        rgbaVolume[i*4+0] = colors[4*mat+0];
        rgbaVolume[i*4+1] = colors[4*mat+1];
        rgbaVolume[i*4+2] = colors[4*mat+2];
        rgbaVolume[i*4+3] = colors[4*mat+3];
    }

    delete []colors;
}
