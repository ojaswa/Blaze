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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QToolButton>
#include <QFuture>
#include <QFutureWatcher>
#include <QtConcurrent>
#include "glwidget.h"
#include "algorithm/volumemanager.h"
#include "ui/dialog1dtransferfunction.h"
#include "ui/dialog2dtransferfunction.h"
#include "ui/dialograycastingsettings.h"
#include "ui/dialogmaterialsettings.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    bool readVolume(QString filename);

signals:
    void volumeDataCreated(VolumeManager *vm);
    void volumeLHComputed(VolumeManager *vm);
    void volumeOcclusionComputed(VolumeManager *vm);
    void volumeEdgesComputed(VolumeManager *vm);
    void materialGraphComputed(VolumeManager *vm);
    void volumePreprocessCompleted(VolumeManager *vm);

private slots:
    void on_action_Read_triggered();
    void on_action_Quit_triggered();
    void on_action1D_TF_toggled(bool arg1);
    void on_action2D_TF_toggled(bool arg1);
    void on_actionRaycasting_settings_toggled(bool arg1);
    void on_actionAbout_triggered();
    void toggleDirectRGBARendering(bool arg1);

    //Asynchronous job callbacks
    void on_volumeReadFinished();
    void on_volumeLHComputed();
    void on_volumeOcclusionComputed();
    void on_volumeEdgesComputed();
    void on_volumePreprocessCompleted();
    void on_classifiedVolumeComputed();
    void on_materialGraphComputed();
    void on_actionCompute_material_graph_triggered();
    void on_actionSave_screenshot_triggered();
    void on_actionSave_RGBA_triggered();
    void on_actionMaterial_settings_toggled(bool arg1);
    void on_actionSave_material_triggered();
    void on_actionRead_material_triggered();

private:
    Ui::MainWindow *ui;
    VolumeManager *m_volumeManager;
    Dialog1DTransferFunction *m_1DTFDialog;
    Dialog2DTransferFunction *m_2DTFDialog;
    DialogRaycastingSettings *m_raycastingSettingsDialog;  
    DialogMaterialSettings *m_materialSettingsDialog;
};

#endif // MAINWINDOW_H
