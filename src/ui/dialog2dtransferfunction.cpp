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

#include "dialog2dtransferfunction.h"
#include "ui_dialog2dtransferfunction.h"

Dialog2DTransferFunction::Dialog2DTransferFunction(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog2DTransferFunction)
{
    ui->setupUi(this);
    setWindowFlags(Qt::Tool);
    setFixedSize(this->width(),this->height());
    show(); m_windowPosition = pos(); hide();
    ui->customPlot->sizePolicy().setHeightForWidth(true);
    ui->customPlot->xAxis->setRange(0.0, 1.0);
    ui->customPlot->yAxis->setRange(0.0, 1.0);


    //Connections
    connect(parent, SIGNAL(volumeLHComputed(VolumeManager*)), this, SLOT(on_volumeLHComputed(VolumeManager*)));
    connect(parent, SIGNAL(volumeOcclusionComputed(VolumeManager*)), this, SLOT(on_volumeOcclusionComputed(VolumeManager*)));
    connect(parent, SIGNAL(materialGraphComputed(VolumeManager*)), this, SLOT(on_materialGraphComputed(VolumeManager*)));
    connect(ui->graphWidget, SIGNAL(updateRGBAInteractionWidgets()), this, SLOT(updateRGBAInteractionWidgets()));

    //Initialization
    ui->radioButton_LH->setEnabled(false);
    ui->radioButton_Occlusion->setEnabled(false);
    ui->radioButton_Graph->setEnabled(false);
    enableRGBAInteractionWidgets(false);

    m_materials = NULL;
    m_intensityRange = NULL;
    m_occlusionRange = NULL;
}

Dialog2DTransferFunction::~Dialog2DTransferFunction()
{
    delete ui;
}

void Dialog2DTransferFunction::on_volumeLHComputed(VolumeManager *vm)
{
    m_densityLH = vm->getLHDensityImage();
    ui->radioButton_LH->setEnabled(true);
    ui->radioButton_LH->setChecked(true);
    if(!ui->radioButton_Occlusion->isChecked())
        setupDisplay(m_densityLH, "L", "H");
}

void Dialog2DTransferFunction::on_volumeOcclusionComputed(VolumeManager *vm)
{
    m_densityOcclusion = vm->getOcclusionDensityImage();
    ui->radioButton_Occlusion->setEnabled(true);
    ui->radioButton_Occlusion->setChecked(true);
    if(!ui->radioButton_LH->isChecked()) {
        setupDisplay(m_densityOcclusion, "Intensity", "Occlusion");
    }
}

void Dialog2DTransferFunction::on_materialGraphComputed(VolumeManager *vm)
{
    ui->radioButton_Graph->setEnabled(true);
    m_materials = &(vm->getMaterials());
    m_graphEdges = &(vm->getGraphEdges());

    m_intensityRange = vm->getIntensityRange();
    m_occlusionRange = vm->getOcclusionRange();
    if(!ui->radioButton_LH->isChecked()) {
        ui->customPlot->clearPlottables();
        setupDisplay(m_densityOcclusion, "Intensity", "Occlusion");
    }

    ui->graphWidget->setupGraph(m_materials, m_graphEdges);
}

void Dialog2DTransferFunction::setupDisplay(double* density, QString xlabel, QString ylabel)
{
    // Initialize customPlot
    // Configure axis rect:
    //ui->customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    ui->customPlot->axisRect()->setupFullAxesBox(true);
    ui->customPlot->xAxis->setLabel(xlabel);
    ui->customPlot->yAxis->setLabel(ylabel);
    ui->customPlot->yAxis->setScaleRatio(ui->customPlot->xAxis, 1.0);

    // Set up the QCPColorMap
    QCPColorMap *colorMap = new QCPColorMap(ui->customPlot->xAxis, ui->customPlot->yAxis);
    int nx = 256;
    int ny = 256;
    colorMap->data()->setSize(nx, ny); // We want the color map to have nx * ny data points
    colorMap->data()->setRange(QCPRange(0.0, 1.0), QCPRange(0.0, 1.0)); // And span the coordinate range 0..1 in both key (x) and value (y) dimensions
    // Now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double z;
    unsigned char alpha;
    for (int xIndex=0; xIndex<nx; ++xIndex)
      for (int yIndex=0; yIndex<ny; ++yIndex)
      {
        //z = log(double(m_densityLH[xIndex + yIndex*nx]) + 1.0);
        z = density[xIndex + yIndex*nx];
        colorMap->data()->setCell(xIndex, yIndex, z);
        alpha = (unsigned char)(255.0*(2.0/(1.0 + exp(-exp(z) + 1.0)) - 1.0));
        colorMap->data()->setAlpha(xIndex, yIndex, alpha);
      }

    // Set the color gradient of the color map to one of the presets
    colorMap->setGradient(QCPColorGradient(QCPColorGradient::gpSpectrum).inverted());

    // Rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient
    colorMap->rescaleDataRange();

    // Make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->customPlot);
    ui->customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    // Rescale the key (x) and value (y) axes so the whole color map is visible:
    ui->customPlot->rescaleAxes();

    ui->customPlot->replot();
}

void Dialog2DTransferFunction::on_radioButton_LH_toggled(bool checked)
{
    ui->stackedWidget->setCurrentIndex(0); //Display the customPlot widget
    enableRGBAInteractionWidgets(false);

    if(checked){
        ui->customPlot->clearPlottables();
        setupDisplay(m_densityLH, "L", "H");
    }
}

void Dialog2DTransferFunction::on_radioButton_Occlusion_toggled(bool checked)
{
    ui->stackedWidget->setCurrentIndex(0); //Display the customPlot widget
    enableRGBAInteractionWidgets(false);
    if(checked) {
        ui->customPlot->clearPlottables();
        setupDisplay(m_densityOcclusion, "Intensity", "Occlusion");
    }
}

void Dialog2DTransferFunction::on_radioButton_Graph_toggled(bool checked)
{
    ui->stackedWidget->setCurrentIndex(1); //Display the graphWidget widget
    enableRGBAInteractionWidgets(true);
}

void Dialog2DTransferFunction::on_screenshotToolButton_clicked()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                            tr("Save Plot"), QCoreApplication::applicationDirPath(),
                                            tr("Image Files (*.pdf *.png *.jpg *.bmp);;Adobe PDF(*.pdf);; PNG (*.png);;JPEG (*.jpg);; Bitmap (*.bmp)"));
    if(filename.isEmpty() || filename.isNull())
        return;
    QFileInfo fi(filename);

    if (ui->stackedWidget->currentIndex() == 0) { // Save customPlot screenshot
        if (fi.suffix() == "pdf")
            ui->customPlot->savePdf(filename);
        else if(fi.suffix() == "png")
            ui->customPlot->savePng(filename, 0, 0, 2.0);
        else if(fi.suffix() == "jpg")
            ui->customPlot->saveJpg(filename, 0, 0, 2.0);
        else if(fi.suffix() == "bmp")
            ui->customPlot->saveBmp(filename, 0, 0, 2.0);
    } else if (ui->stackedWidget->currentIndex() == 1) { // Save graphWidget screenshot
        QPixmap pix = ui->graphWidget->grab();
        pix.save(filename);
    }
}

void Dialog2DTransferFunction::on_renderUpdateToolButton_clicked()
{
    emit updateRGBAVolume();
}

void Dialog2DTransferFunction::enableRGBAInteractionWidgets(bool visible)
{
    int matIdx = ui->graphWidget->selectedMaterial;

    ui->opacityLabel->setEnabled(visible);
    ui->renderUpdateToolButton->setEnabled(visible);
    ui->opacityValueLabel->setEnabled(visible);

    if(!visible) {
        ui->opacityCheckBox->setEnabled(false);
        ui->opacityHorizontalSlider->setEnabled(false);
    } else
        updateRGBAInteractionWidgets();
}

void Dialog2DTransferFunction::updateRGBAInteractionWidgets()
{
    int matIdx = ui->graphWidget->selectedMaterial;
    if(matIdx == -1) {
        ui->opacityCheckBox->setEnabled(false);
        ui->opacityHorizontalSlider->setEnabled(false);

        ui->opacityCheckBox->setChecked(false);
        ui->opacityHorizontalSlider->setValue(0);
        ui->opacityValueLabel->setText(QString::number(0));
    } else {
        Material *mat = m_materials->at(matIdx);
        ui->opacityCheckBox->setEnabled(true);
        ui->opacityHorizontalSlider->setEnabled(mat->m_visible);

        ui->opacityCheckBox->setChecked(mat->m_visible);
        ui->opacityHorizontalSlider->setValue(mat->m_color.a);
        ui->opacityValueLabel->setText(QString::number(mat->m_color.a));
    }
}

void Dialog2DTransferFunction::on_opacityHorizontalSlider_valueChanged(int value)
{
    int matIdx = ui->graphWidget->selectedMaterial;
    if(matIdx < 0) return;

    Material *mat = m_materials->at(matIdx);
    mat->m_color.a = value;
    ui->opacityValueLabel->setText(QString::number(value));
}

void Dialog2DTransferFunction::on_opacityCheckBox_toggled(bool checked)
{
    int matIdx = ui->graphWidget->selectedMaterial;
    if(matIdx < 0) return;

    Material *mat = m_materials->at(matIdx);
    mat->m_visible = checked;
    ui->opacityHorizontalSlider->setEnabled(checked);
    ui->graphWidget->viewport()->repaint();
}
