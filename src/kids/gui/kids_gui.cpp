#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QShowEvent>

#include "../sweepdata.h"
#include "kids_gui.h"
#include "plot.h"
#include <gitversion.h>
#include <ui_kids_gui.h>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    setWindowTitle(QString("KIDs GUI (commit@%1 %2)")
                       .arg(GIT_REVISION)
                       .arg(BUILD_TIMESTAMP));
    // uri root with history
    setUpDataRootUrlModelview();
    // uri root -> data file listing
    setUpDataFileModelView();
    // data file listing -> loaded data
    setUpDataModelView();
    // plot
    connect(ui->leftPlotCanvas, &SweepDataView::mouseMessage,
            [ui=ui](const QString& msg){
                ui->statusbar->showMessage(msg, 0);
            });
    ui->leftPlotControlView->setModel(ui->leftPlotCanvas->graphListModel());
    ui->leftPlotControlView->horizontalHeader()->setSectionResizeMode(
        QHeaderView::ResizeToContents);

    // done
    qDebug() << "mainwindow initialized";
}

MainWindow::~MainWindow() {
    delete ui;
    qDebug() << "mainwindow destroyed";
}

void MainWindow::showEvent(QShowEvent *event) {
    QMainWindow::showEvent(event);
    QApplication::processEvents();
    if (!startUp || event->spontaneous())
        return;
    // run statup setting
    startUp = false;
    // setup initial directory
    QStringList args = QCoreApplication::arguments();
    if (args.count() > 2)
        qDebug() << "ignore extra commandline arguments:" << args.mid(2);
    if (args.count() > 1) {
        QString dir = args.at(1);
        QMetaObject::invokeMethod(this,
                                  [this, dir] { this->setDataRootUrl(dir); },
                                  Qt::ConnectionType::QueuedConnection);
    }
    return;
}

void MainWindow::setDataRootUrl() {
    QUrl url = QFileDialog::getExistingDirectoryUrl(
        this, tr("Open Directory"), QString(),
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks |
            QFileDialog::DontUseNativeDialog);
    QApplication::processEvents();
    if (!url.isEmpty()) {
        QMetaObject::invokeMethod(this, [this, url] { setDataRootUrl(url); },
                                  Qt::ConnectionType::QueuedConnection);
    }
}

void MainWindow::setDataRootUrl(const QUrl &url) {
    qDebug() << "set data root to " << url;
    // dataRootUrlList.add(url);
    dataRootUrlList.insertRows(0, 1);
    dataRootUrlList.setData(dataRootUrlList.index(0), url.toString());
}

void MainWindow::loadDataFile(const QFileInfo &fileInfo) {
    DataItem *dataItemRoot =
        reinterpret_cast<DataItem *>(dataList.invisibleRootItem());
    DataItem *dataItem = nullptr;
    qDebug() << "open data file " << fileInfo;
    try {
        using logging::timeit;
        auto sweepdata =
            timeit("open data file", kids::SweepData::fromPath,
                   fileInfo.filePath().toStdString(), kids::SweepKind::Any);
        QVariant v;
        v.setValue(vSweepData(std::move(sweepdata), fileInfo));
        dataItem = new DataItem(std::move(v));
    } catch (const std::runtime_error &e) {
        QMessageBox::critical(this, "Exception", e.what());
    }
    if (dataItem != nullptr) {
        dataItemRoot->insertRow(0, dataItem);
    }
}

void MainWindow::processDataItem(DataItem* dataItem) {
    auto data = dataItem->data(Qt::UserRole + DataItem::Role::Data);
    if (data.canConvert<vSweepData>()) {
        vSweepData v = data.value<vSweepData>();
        qDebug() << "get sweep from variant" << v.text();
        processData(v);
    } else if (data.canConvert<vSweepTargetData>()) {
        vSweepTargetData v = data.value<vSweepTargetData>();
        qDebug() << "get sweep target from variant" << v.text();
        processData(v);
    } else {
        qDebug() << "process not implemented for" << data.type();
    }
}

void MainWindow::processData(const vSweepData &data) {
    qDebug() << "process sweepdata" << data.text();
    ui->leftPlotCanvas->plot(data.data());
}

void MainWindow::processData(const vSweepTargetData &data) {
    qDebug() << "process sweep target" << data.text();
    using logging::timeit;
    constexpr auto model = kids::SweepModel::S21;
    auto fitter = kids::rsnfitter<model, true>();
    auto results = timeit("fit target", fitter,
                         data.data().freqpack,
                         data.data().datapack);
    SPDLOG_DEBUG(
                "fit result{}", results
                );
}

void MainWindow::setUpDataRootUrlModelview()
{
    // connect data root history to the combo box
    ui->dataRootUrlEdit->setModel(&dataRootUrlList);
    connect(&dataRootUrlList, &QStringListModel::dataChanged,
            [ui = ui]() { ui->dataRootUrlEdit->setCurrentIndex(0); });

    // data file listing -> double click dir -> set new root
    connect(ui->dataFileTreeView, &QTreeView::doubleClicked,
            [this, fileList = &dataFileList]
            (const QModelIndex &index) {
                QFileInfo fi = fileList->filePath(index);
                if (fi.isDir()) {
                    setDataRootUrl(QUrl(fi.filePath()));
                }
            });
}

void MainWindow::setUpDataFileModelView()
{
    ui->dataFileTreeView->header()->setSectionResizeMode(
        QHeaderView::ResizeToContents);
    // 3rd col is modification time
    ui->dataFileTreeView->sortByColumn(3, Qt::DescendingOrder);
    // connect data file model to the view
    ui->dataFileTreeView->setModel(&dataFileList);
    // response to root url setting
    connect(ui->dataRootUrlEdit, &QComboBox::currentTextChanged,
            [model = &dataFileList, ui=ui](const QString &url) {
                QString path = QUrl(url).path();
                model->setRootPath(path);
                if (!path.isEmpty()) {
                    qDebug() << "set root path" << path;
                    const QModelIndex rootIndex =
                        model->index(QDir::cleanPath(path));
                    if (rootIndex.isValid()) {
                        ui->dataFileTreeView->setRootIndex(rootIndex);
                        qDebug() << "view at index:" << rootIndex;
                    } else {
                        qDebug() << "invalid path:" << rootIndex;
                    }
                } else {
                    qDebug() << "empty path";
                }
            });

}

void MainWindow::setUpDataModelView()
{

    // setup the data list model
    dataList.setItemPrototype(new DataItem());
    // used to check existing loaded data
    dataFilter.setSourceModel(&dataList);
    // set model to data list view
    ui->dataTreeView->setModel(&dataList);

    // data file listing -> double click file -> load data
    connect(ui->dataFileTreeView, &QTreeView::doubleClicked,
            [this, fileList = &dataFileList,
             ui = ui](const QModelIndex &index) {
                QFileInfo fi = fileList->filePath(index);
                if (fi.isFile()) {
                    // check if already opened
                    this->dataFilter.setFilterFixedString(fi.filePath());
                    this->dataFilter.setFilterRole(Qt::UserRole +
                                                   DataItem::Role::FilePath);
                    if (this->dataFilter.hasChildren()) {
                        qDebug() << "already opened:" << fi;
                        // highlight the opened
                        ui->dataTreeView->collapseAll();
                        auto index = this->dataFilter.mapToSource(
                            this->dataFilter.index(0, 0));
                        ui->dataTreeView->setExpanded(index, true);
                        auto select = ui->dataTreeView->selectionModel();
                        select->setCurrentIndex(
                            index, QItemSelectionModel::ClearAndSelect |
                                       QItemSelectionModel::Rows);
                    } else {
                        // load data
                        this->loadDataFile(fi);
                    }
                }
            });
    // data list view track data model update by selecting the newly added item
    connect(&dataList, &QStandardItemModel::rowsInserted,
            [this, ui = ui](const QModelIndex &parent, int first, int) {
                auto index = this->dataList.index(first, 0, parent);
                ui->dataTreeView->selectionModel()->setCurrentIndex(
                    index, QItemSelectionModel::ClearAndSelect |
                               QItemSelectionModel::Rows);
                ui->dataTreeView->collapseAll();
                ui->dataTreeView->setExpanded(index, true);
            });
    // data list -> double click -> process data
    connect(ui->dataTreeView, &QTreeView::doubleClicked,
            [this](const QModelIndex &index) {
                // get item from index
                DataItem *item = reinterpret_cast<DataItem *>(
                    this->dataList.itemFromIndex(index));
                qDebug() << "process data item" << item;
                this->processDataItem(item);
            });
}

vSweepData::vSweepData(kids::SweepData data, const QFileInfo &fileInfo)
    : m_data(std::move(data)), m_fileInfo(fileInfo) {
    qDebug() << "vsweepdata from moved data" <<
        QString::fromStdString(fmt::format("{}", m_data));
}

vSweepData::vSweepData() { qDebug() << "vSweepData default constructed"; }

vSweepData::vSweepData(const vSweepData &other)
    : vSweepData(other.data().copy(), other.fileInfo()) {
    qDebug() << "vSweepData from copied data" <<
            QString::fromStdString(fmt::format("{}", m_data));
}

vSweepData::vSweepData(vSweepData &&) {
    qDebug() << "VSweepData from moved other";
}

const kids::SweepData &vSweepData::data() const {
    SPDLOG_TRACE("get data{}", m_data);
    return m_data;
}

QList<vSweepTargetData> vSweepData::targets() const {
    QList<vSweepTargetData> t;
    // for (int i = 0; i < 3; ++i) {
    for (int i = 0; i < data().freqpack.cols(); ++i ) {
        // create sweepdata
        kids::SweepData targetdata{data().kind, data().freqpack.col(i),
                                   data().datapack.col(i)};
        qDebug() << QString::fromStdString(
            fmt::format("target {}: {}", i, targetdata));
        vSweepTargetData v{std::move(targetdata), fileInfo(), i};
        t.append(std::move(v));
        qDebug() << "added target" << t[i].text();
    }
    return t;
}

QString vSweepData::text() const {
    return QString::fromStdString(fmt::format(
        "{} [{}, {}, {}]", fileInfo().baseName().toStdString(),
        // fileInfo().birthTime().toString("hh:mm:ss.zzz").toStdString(),
        m_data.kind, m_data.freqpack.rows(), m_data.freqpack.cols()));
}

const QFileInfo &vSweepData::fileInfo() const { return m_fileInfo; }

vSweepTargetData::vSweepTargetData(kids::SweepData data,
                                   const QFileInfo &fileInfo, int index)
    : m_data(std::move(data)), m_fileInfo(fileInfo), m_index(index),
      fmin(m_data.freqpack.minCoeff() / 1e6), fmax(m_data.freqpack.maxCoeff() / 1e6)
{
    qDebug() << "vsweeptargetdata from moved data";
    qDebug() << "file:" << fileInfo.fileName() << "target:" << index;
    qDebug() << this->text();
}

QString vSweepTargetData::text() const {
    return QString::fromStdString(
        fmt::format("{}: [{:5.2f}, {:5.2f}] MHz [{}, {}]", index(), fmin, fmax,
                    data().kind, data().freqpack.rows()));
}

const QFileInfo &vSweepTargetData::fileInfo() const { return m_fileInfo; }
