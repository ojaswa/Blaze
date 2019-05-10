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

#ifndef DIALOG2DTRANSFERFUNCTION_H
#define DIALOG2DTRANSFERFUNCTION_H

#include <QDialog>
#include "qcustomplot.h"
#include "algorithm/volumemanager.h"
#include "algorithm/material.h"
#include <map>

using namespace std;

namespace Ui {
class Dialog2DTransferFunction;
}

class Dialog2DTransferFunction : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog2DTransferFunction(QWidget *parent = 0);
    ~Dialog2DTransferFunction();
    QPoint m_windowPosition;

protected slots:
    void on_volumeLHComputed(VolumeManager *vm);
    void on_volumeOcclusionComputed(VolumeManager *vm);
    void on_materialGraphComputed(VolumeManager *vm);
    void setupDisplay(double* density, QString xlabel, QString ylabel);

private slots:
    void on_radioButton_LH_toggled(bool checked);
    void on_radioButton_Occlusion_toggled(bool checked);
    void on_screenshotToolButton_clicked();
    void on_renderUpdateToolButton_clicked();
    void on_radioButton_Graph_toggled(bool checked);
    void updateRGBAInteractionWidgets();
    void on_opacityHorizontalSlider_valueChanged(int value);
    void on_opacityCheckBox_toggled(bool checked);

private:
    void enableRGBAInteractionWidgets(bool visible);

    //Member variables
    Ui::Dialog2DTransferFunction *ui;
    double *m_densityLH;
    double *m_densityOcclusion;
    vector<Material*> *m_materials;
    vector<QPoint*> *m_graphEdges;
    double *m_intensityRange, *m_occlusionRange;
    map<void*, int> m_polyIDMap;

signals:
    void updateRGBAVolume();
};

#endif // DIALOG2DTRANSFERFUNCTION_H
