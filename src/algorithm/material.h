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

#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

using namespace std;

struct Color {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
};

struct Material
{
    double m_mean;
    double m_sigma;
    Color m_color;
    bool m_visible;
    bool m_isBackgroud;

    //vector<vector<double> > m_occlusionSpectrumBoundaryI;//Multiple contours of a material - Intensity values
    //vector<vector<double> > m_occlusionSpectrumBoundaryO;//Multiple contours of a material - Occlusion values
    vector<double> m_occlusionSpectrumModeI;
    vector<double> m_occlusionSpectrumModeO;
};

#endif // MATERIAL_H
