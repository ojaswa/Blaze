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

#include "mainwindow.h"
#include "ui_mainwindow.h"

//Qt includes
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //this->setUnifiedTitleAndToolBarOnMac(true);
    m_volumeManager = new VolumeManager();
    m_1DTFDialog = new Dialog1DTransferFunction(this);
    m_2DTFDialog = new Dialog2DTransferFunction(this);
    m_raycastingSettingsDialog = new DialogRaycastingSettings(this);
    m_materialSettingsDialog = new DialogMaterialSettings(this);

    //Enable/disable actions
    ui->action_Read->setEnabled(true);
    ui->action1D_TF->setEnabled(false);
    ui->action2D_TF->setEnabled(false);
    ui->actionRaycasting_settings->setEnabled(false);
    ui->actionMaterial_settings->setEnabled(false);
    ui->actionCompute_material_graph->setEnabled(false);
    ui->actionSave_RGBA->setEnabled(false);
    ui->actionSave_material->setEnabled(false);
    ui->actionRead_material->setEnabled(false);

    //Connections
    connect(m_1DTFDialog, SIGNAL(finished(int)), ui->action1D_TF, SLOT(toggle()));
    connect(m_1DTFDialog, SIGNAL(TFChanged(unsigned char*)), ui->centralWidget, SLOT(TF1DChanged(unsigned char*)));
    connect(m_2DTFDialog, SIGNAL(finished(int)), ui->action2D_TF, SLOT(toggle()));
    connect(m_2DTFDialog, SIGNAL(updateRGBAVolume()), ui->centralWidget, SLOT(on_classifiedVolumeComputed()));
    connect(m_raycastingSettingsDialog, SIGNAL(finished(int)), ui->actionRaycasting_settings, SLOT(toggle()));
    connect(m_raycastingSettingsDialog, SIGNAL(stepSizeChanged(float)), ui->centralWidget, SLOT(raycasterStepSizeChanged(float)));
    connect(m_raycastingSettingsDialog, SIGNAL(interpolationTypeChanged(RaycastingInterpolationType)), ui->centralWidget, SLOT(raycasterInterpolationTypeChanged(RaycastingInterpolationType)));
    connect(m_raycastingSettingsDialog, SIGNAL(enableJitteredSampling(bool)), ui->centralWidget, SLOT(enableJitteredSampling(bool)));
    connect(m_raycastingSettingsDialog, SIGNAL(toggleDirectRGBARendering(bool)), ui->centralWidget, SLOT(toggleDirectRGBARendering(bool)));
    connect(m_raycastingSettingsDialog, SIGNAL(toggleDirectRGBARendering(bool)), this, SLOT(toggleDirectRGBARendering(bool)));
    connect(m_raycastingSettingsDialog, SIGNAL(togglePhongShading(bool)), ui->centralWidget, SLOT(togglePhongShading(bool)));
    connect(m_volumeManager, SIGNAL(volumeDataCreated(VolumeManager *)), this, SLOT(on_volumeReadFinished()));
    connect(m_volumeManager, SIGNAL(volumeLHComputed(VolumeManager *)), this, SLOT(on_volumeLHComputed()));
    connect(m_volumeManager, SIGNAL(volumeOcclusionComputed(VolumeManager *)), this, SLOT(on_volumeOcclusionComputed()));
    connect(m_volumeManager, SIGNAL(volumeEdgesComputed(VolumeManager*)), this, SLOT(on_volumeEdgesComputed()));
    connect(m_volumeManager, SIGNAL(volumePreprocessCompleted(VolumeManager*)), this, SLOT(on_volumePreprocessCompleted()));
    connect(m_volumeManager, SIGNAL(classifiedVolumeComputed(VolumeManager*)), this, SLOT(on_classifiedVolumeComputed()));
    connect(m_volumeManager, SIGNAL(classifiedVolumeComputed(VolumeManager*)), ui->centralWidget, SLOT(on_classifiedVolumeComputed()));
    connect(m_volumeManager, SIGNAL(materialGraphComputed(VolumeManager*)), this, SLOT(on_materialGraphComputed()));
    connect(m_volumeManager, SIGNAL(volumeGradientComputed(VolumeManager*)), ui->centralWidget, SLOT(on_volumeGradientComputed()));
    connect(m_materialSettingsDialog, SIGNAL(finished(int)), ui->actionMaterial_settings, SLOT(toggle()));
    connect(m_volumeManager, SIGNAL(classifiedVolumeRead(VolumeManager*)), ui->centralWidget, SLOT(on_classifiedVolumeRead()));
}

MainWindow::~MainWindow()
{
    delete m_raycastingSettingsDialog;
    delete m_materialSettingsDialog;
    delete m_1DTFDialog;
    delete m_2DTFDialog;
    delete m_volumeManager;
    delete ui;
}

void MainWindow::on_volumeReadFinished()
{
    //Begin preprocessing asynchronously while allowing user to explore the volume render.
    QtConcurrent::run(m_volumeManager, &VolumeManager::preprocess); //Preprocess volume

    //Resume UI thread jobs
    emit volumeDataCreated(m_volumeManager);
    ui->action1D_TF->setEnabled(true);
    ui->actionRaycasting_settings->setEnabled(true);
}

void MainWindow::on_volumeLHComputed()
{
    emit volumeLHComputed(m_volumeManager);;
    ui->action2D_TF->setEnabled(true);
}

void MainWindow::on_volumeOcclusionComputed()
{
    emit volumeOcclusionComputed(m_volumeManager);;
    ui->action2D_TF->setEnabled(true);
}

void MainWindow::on_volumePreprocessCompleted()
{
    emit volumePreprocessCompleted(m_volumeManager);
    ui->actionCompute_material_graph->setEnabled(true);
    ui->actionMaterial_settings->setEnabled(true);
}

void MainWindow::on_materialGraphComputed()
{
    emit materialGraphComputed(m_volumeManager);
    //TODO: Enable disable UI elements that depend on completion of material graph computation.
}

void MainWindow::on_volumeEdgesComputed()
{
    emit volumeEdgesComputed(m_volumeManager);
}

void MainWindow::on_action_Read_triggered()
{
    QString selfilter = tr("NRRD (*.nhdr *.nrrd)");
    QString filename = QFileDialog::getOpenFileName(this, "Open a volume", QCoreApplication::applicationDirPath(),
                                              tr("Supported formats (*.nhdr *.nrrd *.vtk);;NRRD (*.nhdr *.nrrd);;VTK (*.vtk)"),
                                              &selfilter);
    if(filename.isEmpty() || filename.isNull())
        return;

    readVolume(filename);
    ui->action_Read->setEnabled(false);
    ui->actionRead_material->setEnabled(true);
}

bool MainWindow::readVolume(QString filename)
{
    //Before  reading a new volume file, clear old contents of the process dir.
    QDir dir( "./process");
    if(dir.exists()) {
        dir.setFilter( QDir::NoDotAndDotDot | QDir::Files );
        foreach( QString dirItem, dir.entryList() )
            dir.remove( dirItem );

        dir.setFilter( QDir::NoDotAndDotDot | QDir::Dirs );
        foreach( QString dirItem, dir.entryList() )
        {
            QDir subDir( dir.absoluteFilePath( dirItem ) );
            subDir.removeRecursively();
        }
    }
    else
        QDir().mkdir(dir.path());

    //Return true of false appropriatey if the volume was successfully read or not.
    QFileInfo fi(filename);
    QString filetype = fi.suffix();

    fprintf(stderr, "Reading %s from disk...\n", filename.toStdString().c_str());
    if (filetype.compare("nhdr")==0 | filetype.compare("nrrd")==0) //Load NRRD file
    {
        m_volumeManager->readNHDR(filename.toStdString().c_str());
    }
    else if (filetype.compare("vtk") == 0) // Load VTK file
    {
        //TODO:
    }

    return true;
}


void MainWindow::on_action_Quit_triggered()
{
    QApplication::quit();
}

void MainWindow::on_action1D_TF_toggled(bool arg1)
{
    if(arg1) {
        m_1DTFDialog->move(m_1DTFDialog->m_windowPosition);
        m_1DTFDialog->show();
    }
    else {
        m_1DTFDialog->m_windowPosition = m_1DTFDialog->pos();
        m_1DTFDialog->hide();
    }
}

void MainWindow::on_action2D_TF_toggled(bool arg1)
{
    if(arg1) {
        m_2DTFDialog->move(m_2DTFDialog->m_windowPosition);
        m_2DTFDialog->show();
    }
    else {
        m_2DTFDialog->m_windowPosition = m_2DTFDialog->pos();
        m_2DTFDialog->hide();
    }
}

void MainWindow::on_actionRaycasting_settings_toggled(bool arg1)
{
    if(arg1) {
        m_raycastingSettingsDialog->move(m_raycastingSettingsDialog->m_windowPosition);
        m_raycastingSettingsDialog->show();
    }
    else {
        m_raycastingSettingsDialog->m_windowPosition = m_raycastingSettingsDialog->pos();
        m_raycastingSettingsDialog->hide();
    }
}

void MainWindow::on_actionMaterial_settings_toggled(bool arg1)
{
    if(arg1) {

        m_materialSettingsDialog->move(m_materialSettingsDialog->m_windowPosition);
        m_materialSettingsDialog->show();
    }
    else {
        m_materialSettingsDialog->m_windowPosition = m_materialSettingsDialog->pos();
        m_materialSettingsDialog->hide();
    }
}

void MainWindow::on_actionCompute_material_graph_triggered()
{
    ui->actionCompute_material_graph->setEnabled(false);
    //Run material graph computation in a seperate thread so that user can still interact with UI
    int LH_mcs = m_materialSettingsDialog->getLHMinClusterSize();
    int LH_ms = m_materialSettingsDialog->getLHMinSamples();
    bool restore = m_materialSettingsDialog->getRestoreFlag();
    bool onlyAlgo = false;
    QtConcurrent::run(m_volumeManager, &VolumeManager::computeMaterialGraph, LH_mcs, LH_ms, restore,onlyAlgo);
}

void MainWindow::on_classifiedVolumeComputed()
{
    ui->action1D_TF->setEnabled(false); //Rendering should use RGBA volume instead
    ui->actionSave_RGBA->setEnabled(true);
    ui->actionSave_material->setEnabled(true);
    m_1DTFDialog->m_windowPosition = m_1DTFDialog->pos();
    m_1DTFDialog->hide();
    m_raycastingSettingsDialog->enableDirectRenderingToggle();
}

void MainWindow::toggleDirectRGBARendering(bool arg1)
{
    if (arg1) {
        ui->action1D_TF->setEnabled(false);
        if(!m_1DTFDialog->isHidden()) {
            m_1DTFDialog->m_windowPosition = m_1DTFDialog->pos();
            m_1DTFDialog->hide();
        }
    }
    else
        ui->action1D_TF->setEnabled(true);
}


void MainWindow::on_actionSave_screenshot_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                            tr("Save volume render"), QCoreApplication::applicationDirPath(),
                                            tr("PNG (*.png)"));
    if(filename.isEmpty() || filename.isNull())
        return;
    QFileInfo fi(filename);

    QImage img = ui->centralWidget->grabFramebuffer();
    img.save(filename);
}

void MainWindow::on_actionSave_RGBA_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                            tr("Save RGB volume"), QCoreApplication::applicationDirPath(),
                                            tr("VTK (*.vtk)"));
    if(filename.isEmpty() || filename.isNull())
        return;
    //m_volumeManager->saveRGBA(filename); //Save RGB instead
    m_volumeManager->saveRGB(filename);
}

void MainWindow::on_actionSave_material_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                            tr("Save material volume"), QCoreApplication::applicationDirPath(),
                                            tr("VTK (*.vtk)"));
    if(filename.isEmpty() || filename.isNull())
        return;
    m_volumeManager->saveClassifiedVolume(filename);
}

void MainWindow::on_actionRead_material_triggered()
{
    QString selfilter = tr("NRRD (*.nhdr *.nrrd)");
    QString filename = QFileDialog::getOpenFileName(this,
                                            tr("Open material volume"), QCoreApplication::applicationDirPath(),
                                            tr("Supported formats (*.nhdr *.nrrd *.vtk);;NRRD (*.nhdr *.nrrd);;VTK (*.vtk)"));
                                            //&selfilter);
    if(filename.isEmpty() || filename.isNull())
        return;
    if(m_volumeManager->readClassifiedVolume(filename)) {
        ui->action1D_TF->setEnabled(false);
        ui->action2D_TF->setEnabled(false);

        ui->actionSave_RGBA->setEnabled(false); // Disallow saving the usual route --FIX it later (requires updating materials in VolumeManager)
        ui->actionSave_material->setEnabled(false);

        m_1DTFDialog->m_windowPosition = m_1DTFDialog->pos();
        m_1DTFDialog->hide();
        m_2DTFDialog->m_windowPosition = m_2DTFDialog->pos();
        m_2DTFDialog->hide();

        m_raycastingSettingsDialog->enableDirectRenderingToggle(); //Rendering should use RGBA volume instead
        ui->actionCompute_material_graph->setEnabled(true); //Enable (re)computing the material graph
    }
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("About Blaze"),
             tr(   "<p style=\"text-align: center;font-size: 20px\"><strong>Blaze</strong></p>"
                   "<span style=\"font-weight:normal; text-align: center\">A volume rendering and analytics program</span><br><br>"
                   "<span style=\"font-weight:normal; text-align: center\">Version 2.0</span><br><br>"
                   "<span style=\"font-weight:normal; text-align: center\">Copyright \u00A9 2016-2019 Graphics Research Group</span><br>"
                   "<span style=\"font-weight:normal; text-align: center\">IIIT Delhi</span>"
               ));
}
