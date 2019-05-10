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

#include "dialogmaterialsettings.h"
#include "ui_dialogmaterialsettings.h"

DialogMaterialSettings::DialogMaterialSettings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogMaterialSettings)
{
    ui->setupUi(this);
    setWindowFlags(Qt::Tool);
    setFixedSize(this->width(),this->height());
    show(); m_windowPosition = pos(); hide();

    ui->checkBoxOverrideLH->setChecked(false);
    ui->checkBoxRestorePrevious->setChecked(false);
    ui->spinBoxLHMCS->setEnabled(false);
    ui->spinBoxLHMS->setEnabled(false);
}

DialogMaterialSettings::~DialogMaterialSettings()
{
    delete ui;
}

void DialogMaterialSettings::on_checkBoxOverrideLH_toggled(bool checked)
{
    ui->spinBoxLHMCS->setEnabled(checked);
    ui->spinBoxLHMS->setEnabled(checked);
}

int DialogMaterialSettings::getLHMinClusterSize()
{
    if(ui->checkBoxOverrideLH->isChecked())
        return ui->spinBoxLHMCS->value();
    else
        return -1;
}

int DialogMaterialSettings::getLHMinSamples()
{
    if(ui->checkBoxOverrideLH->isChecked())
        return ui->spinBoxLHMS->value();
    else
        return -1;
}

bool DialogMaterialSettings::getRestoreFlag()
{
    return ui->checkBoxRestorePrevious->isChecked();
}
