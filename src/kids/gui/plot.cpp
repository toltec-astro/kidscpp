#include "plot.h"

SweepDataView::SweepDataView(QWidget *parent): QCustomPlot(parent) {
    this->setInteraction(QCP::iRangeDrag, true) ;
    this->setInteraction(QCP::iRangeZoom, true) ;
    this->setInteraction(QCP::iSelectItems, true);

    this->xAxis->setLabel("Frequency (MHz)");
    this->xAxis->grid()->setVisible(false);

    this->yAxis->setLabel("dB");
    this->yAxis->grid()->setVisible(false);
    this->yAxis->setRange(-20, 5);

    this->yAxis2->setVisible(true);
    this->yAxis2->setLabel("Log Amplitude (arbi. unit)");
    this->yAxis2->setScaleType(QCPAxis::stLogarithmic);
    this->yAxis2->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));
    syncYAxis2(this->yAxis->range());
    connect(this->yAxis, SIGNAL(rangeChanged(QCPRange)),
            this, SLOT(syncXAxis2(QCPRange)));

    title = new QCPItemText(this);
    title->position->setType(QCPItemPosition::ptAxisRectRatio);
    title->position->setCoords(0.03, 0.04);
    title->setSelectable(false);
    title->setPositionAlignment(Qt::AlignLeft |Qt::AlignTop);
    title->setText("N/A");

    stats = new QCPItemText(this);
    stats->position->setParentAnchor(title->position);
    stats->position->setCoords(0, title->bottom->pixelPosition().y() - title->top->pixelPosition().y());
    stats->setSelectable(false);
    stats->setPositionAlignment(Qt::AlignLeft |Qt::AlignTop);
    stats->setTextAlignment(Qt::AlignLeft);
    stats->setText("N/A");

    mouseTracer = new QCPItemTracer(this);
    mouseTracer->setStyle(QCPItemTracer::tsPlus);
    mouseTracer->setPen(QPen(solarized::cyan));
    mouseTracer->setBrush(QBrush(solarized::cyan));
    mouseTracer->setSize(8);
    mouseTracer->setSelectable(false);
    mouseTracer->position->setType(QCPItemPosition::ptPlotCoords);

    mouseSelector = new QCPItemTracer(this);
    mouseSelector->setStyle(QCPItemTracer::tsCircle);
    mouseSelector->setPen(QPen(solarized::base02));
    mouseSelector->setBrush(QBrush(solarized::base02));
    mouseSelector->setSize(8);
    mouseSelector->setSelectable(false);
    mouseSelector->position->setType(QCPItemPosition::ptPlotCoords);

    connect(this, &SweepDataView::mouseMove,
            this, &SweepDataView::mouseTrace);

    // graphs
    /*
    graphList.setItemPrototype(new GraphItem());
    vGraph vg(
        this->addGraph(),
        "sweepdata",
        true, 1.,
        solarized::blue
        );
    graphList.add(std::move(vg));
    */
    // create sweepdata layer
    if (this->addLayer("sweepdata", this->layer("main"))) {
        this->setCurrentLayer("sweepdata");
        this->addGraph();
    }
    qDebug() << "sweep data view initialized";
}

SweepDataView::~SweepDataView() {
    qDebug() << "sweep data view destroyed";
}

void SweepDataView::mouseZoom(QWheelEvent *) {
    if (isZoomX && isZoomY) {
        this->axisRect()->setRangeZoom(Qt::Horizontal|Qt::Vertical);
    }
    else if (isZoomX) {
        this->axisRect()->setRangeZoom(this->xAxis->orientation());
    }
    else if (isZoomY) {
        this->axisRect()->setRangeZoom(this->yAxis->orientation());
    }
    else {
        this->axisRect()->setRangeZoom(0);
    }
}

void SweepDataView::mouseTrace(QMouseEvent *event) {
    double x = this->xAxis->pixelToCoord(event->x());
    double y = this->yAxis->pixelToCoord(event->y());
    QString msg("mouse: %1, %2; nearest: %3, %4");
    msg = msg.arg(QString::number(x, 'f', 3),
                  QString::number(y, 'f', 3));
    mouseTracer->setGraphKey(x);
    msg = msg.arg(QString::number(mouseTracer->position->key(), 'f', 3),
                  QString::number(mouseTracer->position->value(), 'f', 3)
                  );
    this->replot(QCustomPlot::rpQueuedRefresh);
    emit mouseMessage(msg);
}

void SweepDataView::mouseSelect(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton)
    {
        mouseSelector->position->setCoords(mouseTracer->position->coords());
        this->replot(QCustomPlot::rpImmediateRefresh);
        emit selected(mouseTracer->position->key());
    }
}

void SweepDataView::refresh() {
    qDebug() << "refresh canvas";
    this->replot(QCustomPlot::rpImmediateRefresh);
}

void SweepDataView::syncXAxis2(QCPRange range) {}

void SweepDataView::save() {
    QString fileName;
    QString defaultFileName;
    defaultFileName = QString("fig_%1.png").arg(title->text());
    fileName = QFileDialog::getSaveFileName(this, tr("save figure to"),
                                            defaultFileName, tr("Images (*.png *.jpeg *.jpg *.pdf)"));
    if (fileName.size() > 0)
    {
        //turn off the mouse trace items
        mouseTracer->setVisible(false);
        mouseSelector->setVisible(false);
        bool success;
        QFileInfo fi(fileName);
        if (fi.suffix().toLower() == "png")
        {
            success = this->savePng(fileName);
        }
        else if (fi.suffix().toLower() == "jpeg" || fi.suffix().toLower() == "jpg")
        {
            success = this->saveJpg(fileName);
        }
        else if (fi.suffix().toLower() == "pdf")
        {
            success = this->savePdf(fileName);
        }
        else
        {
            fileName += ".png";
            success = this->savePng(fileName);
        }
        if (success)
        {
            emit savedMessage(QString("Figure saved: %1").arg(fileName), 3000);
        }
        // restore
        mouseTracer->setVisible(true);
        mouseSelector->setVisible(true);
    }
}

void SweepDataView::setZoomX(bool zoomX) { isZoomX = zoomX; }

void SweepDataView::setZoomY(bool zoomY) { isZoomY = zoomY; }

void SweepDataView::setShowCutStats(bool showCutStats) {
    isShowCutStats = showCutStats;
}
