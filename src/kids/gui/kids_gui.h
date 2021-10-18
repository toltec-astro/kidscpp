#pragma once
#include "../rsnfitter.h"
#include "../sweepdata.h"
#include <QDebug>
#include <QFileSystemModel>
#include <QMainWindow>
#include <QSortFilterProxyModel>
#include <QStandardItem>
#include <QStringListModel>
#include <QUrl>

namespace Ui {
class MainWindow;
}

class vSweepTargetData;
class vSweepData {
public:
    vSweepData(kids::SweepData data, const QFileInfo &fileInfo);

    vSweepData();
    vSweepData(const vSweepData &other);
    vSweepData(vSweepData &&other);
    ~vSweepData() = default;

    const kids::SweepData &data() const;

    QList<vSweepTargetData> targets() const;

    QString text() const;

    const QFileInfo &fileInfo() const;

private:
    kids::SweepData m_data;
    QFileInfo m_fileInfo;
};
Q_DECLARE_METATYPE(vSweepData);

// wrapper around
class vSweepTargetData {
public:
    vSweepTargetData(kids::SweepData data, const QFileInfo &fileInfo, int index);
    vSweepTargetData() { qDebug() << "vSweepTargetData default constructed"; }
    vSweepTargetData(const vSweepTargetData &other)
     : vSweepTargetData(other.data().copy(), other.fileInfo(), other.index())
    {
        qDebug() << "vSweepTargetData copied from" << other.text();
        qDebug() << "vSweepTargetData this:" << this->text();
    };
    ~vSweepTargetData() = default;

    const kids::SweepData &data() const {return m_data;}

    int index() const { return m_index; }

    QString text() const;
    const QFileInfo &fileInfo() const;;
private:
    kids::SweepData m_data;
    QFileInfo m_fileInfo;
    int m_index = 0;
    double fmin = 0;
    double fmax = 0;
};

Q_DECLARE_METATYPE(vSweepTargetData);

namespace internal {

// registere QVariant enabled wrapper types
template <typename T> struct vtype_impl {
    using type =
        std::conditional_t<std::is_same_v<std::decay_t<T>, kids::SweepData>,
                           vSweepData, void
                           //               std::conditional_t<
                           //               std::is_same_v<std::decay_t<T>,
                           //          kids::SweepTargetData>, vSweepTargetData,
                           //               void>
                           >;
};
inline void register_vtypes() {
    qRegisterMetaType<vSweepData>();
    qRegisterMetaType<vSweepTargetData>();
}

template <typename T> using vtype = typename vtype_impl<T>::type;
} // namespace internal

class DataItem : public QStandardItem {
public:
    using QStandardItem::QStandardItem;
    DataItem() { internal::register_vtypes(); }
    explicit DataItem(QVariant &&data) : DataItem() {
        qDebug() << "create DataItem from variant:" << data;
        setData(std::move(data), Qt::UserRole + Role::Data);
        // dispatch to setup according to types
        const QVariant &d = this->data(Qt::UserRole + Role::Data);
        qDebug() << "data stored:" << d;
        if (d.canConvert<vSweepData>()) {
            setUp(d.value<vSweepData>());
        } else if (d.canConvert<vSweepTargetData>()) {
            setUp(d.value<vSweepTargetData>());
        } else if (d.canConvert<QWidget *>()) {
            // do not know why it does this
        } else {
            qDebug() << "ignore type" << data.type();
        }
    }
    enum Role { Data = 1, FilePath = 2 };

    void setUp(const vSweepData &v_sweepdata) {
        setText(v_sweepdata.text());
        setData(v_sweepdata.fileInfo().filePath(),
                Qt::UserRole + Role::FilePath);
        qDebug() << "setup vsweepdata" << v_sweepdata.text();
        auto targets = v_sweepdata.targets();
        for (auto &target : targets) {
            qDebug() << "target: " << target.text();
            QVariant v;
            v.setValue(std::move(target));
            this->add(std::move(v));
        }
    }
    void setUp(const vSweepTargetData &v_sweeptargetdata) {
        setText(v_sweeptargetdata.text());
        qDebug() << "setup vsweeptargetdata" << v_sweeptargetdata.text();
        // setData(v_sweepdata.fileInfo().filePath(), Qt::UserRole +
        // Role::Data);
    }

    virtual QStandardItem *clone() const { return new DataItem(); }

    void add(QVariant &&data) {
        qDebug() << "add" << data << "to" << this->text();
        auto item = new DataItem(std::move(data));
        qDebug() << "created data item" << item;
        this->appendRow(item);
    }

private:
    QVariant m_data;
};

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;
    void showEvent(QShowEvent *event) override;

public slots:
    void setDataRootUrl();
    void setDataRootUrl(const QUrl &url);
    void loadDataFile(const QFileInfo &fileInfo);
    void processDataItem(DataItem* dataItem);

private:
    void processData(const vSweepData &data);
    void processData(const vSweepTargetData &data);
    Ui::MainWindow *ui;
    bool startUp = true;  // used to handle startup arguments
    // data models
    QStringListModel dataRootUrlList{};
    QFileSystemModel dataFileList;
    QStandardItemModel dataList{};
    QSortFilterProxyModel dataFilter;
    void setUpDataRootUrlModelview();
    void setUpDataFileModelView();
    void setUpDataModelView();
};
