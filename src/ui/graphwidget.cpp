/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "graphwidget.h"
#include "edge.h"
#include "node.h"

#include <math.h>
#include <ui/dialog2dtransferfunction.h>
#include <QWheelEvent>

#define SCENE_SIZE 100

GraphWidget::GraphWidget(QWidget *parent)
    : QGraphicsView(parent), timerId(0), numIterations(0), selectedMaterial(-1)
{
    QGraphicsScene *scene = new QGraphicsScene(this);
    scene->setItemIndexMethod(QGraphicsScene::NoIndex);
    //scene->setSceneRect(-SCENE_SIZE, -SCENE_SIZE, SCENE_SIZE, SCENE_SIZE);
    setScene(scene);
    setCacheMode(CacheBackground);
    setViewportUpdateMode(BoundingRectViewportUpdate);
    setRenderHint(QPainter::Antialiasing);
    setTransformationAnchor(AnchorUnderMouse);
    setFrameStyle(QFrame::NoFrame);
    //scale(qreal(0.8), qreal(0.8));
    //setMinimumSize(SCENE_SIZE, SCENE_SIZE);
    setWindowTitle(tr("Elastic Nodes"));
}

void GraphWidget::setupGraph(vector<Material *> *nodes, vector<QPoint *> *edges)
{
    QGraphicsScene *scene = this->scene();
    m_materials = nodes;

    //Add all material nodes
    int posx, posy;
    for(int i=0; i<nodes->size(); i++) {
        Node *node = new Node(this, i);

        //Generate random position
        Material *material = nodes->at(i);
        posx = SCENE_SIZE*float(qrand())/float(RAND_MAX);
        posy = SCENE_SIZE*float(qrand())/float(RAND_MAX);
        node->setPos(posx, posy);
        node->setColor((QColor(material->m_color.r, material->m_color.g, material->m_color.b)));

        m_nodes.push_back(node);
        scene->addItem(node);
    }

    //Add all edges
    QPoint *edge;
    for(int i=0; i<edges->size(); i++) {
        edge = edges->at(i);
        scene->addItem(new Edge(m_nodes[edge->x()], m_nodes[edge->y()]));
    }
}

void GraphWidget::itemMoved()
{
    if (!timerId) {
        timerId = startTimer(1000 / 25);
        numIterations = 0;
    }
}

void GraphWidget::timerEvent(QTimerEvent *event)
{
    numIterations++;
    Q_UNUSED(event);

    QList<Node *> nodes;
    foreach (QGraphicsItem *item, scene()->items()) {
        if (Node *node = qgraphicsitem_cast<Node *>(item))
            nodes << node;
    }

    foreach (Node *node, nodes)
        node->calculateForces();

    bool itemsMoved = false;
    foreach (Node *node, nodes) {
        if (node->advance())
            itemsMoved = true;
    }

    if ((!itemsMoved) || (numIterations > MAX_ITERATIONS)) {
        killTimer(timerId);
        timerId = 0;
    }
}

#ifndef QT_NO_WHEELEVENT

void GraphWidget::wheelEvent(QWheelEvent *event)
{
    scaleView(pow((double)2, -event->delta() / 240.0));
}

#endif

void GraphWidget::mouseReleaseEvent(QMouseEvent *event)
{
    //Note: Having this callback here disables graph autoshape functionality. FIXME.
    selectedMaterial =-1; //Deselct any selected material
    viewport()->repaint(); //Force repaint of all nodes
    emit updateRGBAInteractionWidgets();

    QGraphicsView::mouseReleaseEvent(event);
}

void GraphWidget::scaleView(qreal scaleFactor)
{
    qreal factor = transform().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    if (factor < 0.07 || factor > 100)
        return;

    scale(scaleFactor, scaleFactor);
}

void GraphWidget::shuffle()
{
    foreach (QGraphicsItem *item, scene()->items()) {
        if (qgraphicsitem_cast<Node *>(item))
            item->setPos(-150 + qrand() % 300, -150 + qrand() % 300);
    }
}

void GraphWidget::zoomIn()
{
    scaleView(qreal(1.2));
}

void GraphWidget::zoomOut()
{
    scaleView(1 / qreal(1.2));
}
