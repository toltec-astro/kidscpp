#pragma once

#include "kids_gui.h"
#include "qcustomplot/qcustomplot.h"

class vGraph {
public:
    vGraph(QCPGraph* graph, const QString& title,
           bool visible, double width,
           const QColor& color,
           Qt::PenStyle penStyle=Qt::SolidLine,
           QCPGraph::LineStyle lineStyle=QCPGraph::LineStyle::lsLine)
        :m_graph(graph),  m_title(title),
          m_visible(visible), m_width(width), m_color(color),
          m_penStyle(penStyle), m_lineStyle(lineStyle)
                              {}
    vGraph() = default;
    /*
    vGraph(const vGraph& other):
                                  vGraph(other.graph(), other.title(), other.visible(),
                                         other.width(), other.color(), other.penStyle(), other.lineStyle()) {}
*/
    vGraph(const vGraph& other) = default;
    ~vGraph() = default;

    QCPGraph* graph() const {
        return m_graph;
    }
    const QString& title() const {
        return m_title;
    }
    bool visible() const {
        return m_visible;
    }
    double width() const {
        return m_width;
    }
    const QColor& color() const {
        return m_color;
    }
    Qt::PenStyle penStyle() const {
        return m_penStyle;
    }
    QCPGraph::LineStyle lineStyle() const {
        return m_lineStyle;
    }

private:
    QCPGraph* m_graph = nullptr;
    QString m_title = "N/A";
    bool m_visible = true;
    double m_width = 1.;
    QColor m_color;
    Qt::PenStyle m_penStyle;
    QCPGraph::LineStyle m_lineStyle;
};
Q_DECLARE_METATYPE(vGraph);

class GraphItem: public QStandardItem {
public:
    using QStandardItem::QStandardItem;
    GraphItem() {qRegisterMetaType<vGraph>();}
    explicit GraphItem(vGraph&& vgraph): GraphItem() {
        qDebug() << "create GraphItem from graph: " << vgraph.title();
        QVariant v;
        v.setValue(std::move(vgraph));
        setData(std::move(v), Qt::UserRole);
        // setText(vgraph.title());
    }
    /*
    explicit GraphItem(bool visible): GraphItem() {
        setData(visible, Qt::CheckStateRole);
    }
    explicit GraphItem(int width): GraphItem() {
        setData(width, Qt::UserRole + Role::Width);
    }
    explicit GraphItem(const QColor& color): GraphItem() {
        setData(color, Qt::BackgroundRole);
    }
    explicit GraphItem(Qt::PenStyle penStyle): GraphItem() {
        QVariant v;
        v.setValue(penStyle);
        setData(v, Qt::UserRole + Role::PenStyle);
    }
    explicit GraphItem(QCPGraph::LineStyle lineStyle): GraphItem() {
        QVariant v;
        v.setValue(lineStyle);
        setData(v, Qt::UserRole + Role::LineStyle);
    }
    */

    // DisplayRole, BackgroundRole
    // enum Role { Data, Width, PenStyle, LineStyle };
    virtual QStandardItem *clone() const { return new GraphItem(); }
};

class GraphListModel: public QStandardItemModel {
    Q_OBJECT
public:
    using QStandardItemModel::QStandardItemModel;
    enum Column {
        Title,
        Visible,
        Width,
        Color,
        PenStyle,
        LineStyle,
        __COUNT__
    };
    Q_ENUM(Column);
    void add(vGraph&& vgraph) {
        auto item = new GraphItem(std::move(vgraph));
        QList<QStandardItem*> row;
        for (int i = 0; i < this->columnCount(QModelIndex()); ++i) {
            row.append(item);
        }
        this->appendRow(row);
    }

    int columnCount(const QModelIndex &/*parent*/) const override
    {
        return static_cast<int>(Column::__COUNT__);
    };
    /*
    Qt::ItemFlags flags(const QModelIndex &index) const override
    {
        if (!index.isValid())
            return 0;
        return Qt::ItemIsEditable | Q::flags(index) |Qt::ItemIsUserCheckable;
    }
    */
    QVariant data(
        const QModelIndex &index, int role = Qt::DisplayRole) const override
    {
        if (!index.isValid())
            return QVariant();
        GraphItem* item = reinterpret_cast<GraphItem*>(this->item(index.row(), 0));
        QVariant v = item->data(Qt::UserRole);
        vGraph vgraph = v.value<vGraph>();

        switch (static_cast<Column>(index.column())) {
        case Column::Title: {
            if (role != Qt::DisplayRole) {
                return QVariant();
            }
            return vgraph.title();
        }
        case Column::Visible: {
            if (role != Qt::CheckStateRole) {
                return QVariant();
            }
            return vgraph.visible()?Qt::Checked:Qt::Unchecked;
        }
        case Column::Width: {
            if (role != Qt::DisplayRole) {
                return QVariant();
            }
            return QString("%1").arg(vgraph.width());
        }
        case Column::Color: {
            if (role != Qt::BackgroundRole) {
                return QVariant();
            }
            return QBrush(vgraph.color());
        }
        case Column::PenStyle: {
            if (role != Qt::DisplayRole) {
                return QVariant();
            }
            QMetaEnum metaEnum = QMetaEnum::fromType<Qt::PenStyle>();
            return metaEnum.valueToKey(vgraph.penStyle());
        }
        case Column::LineStyle: {
            if (role != Qt::DisplayRole) {
                return QVariant();
            }
            return vgraph.lineStyle();
        }
        default: {
            qDebug() << "invalid column";
        }
        }
        qDebug() << "do not know how to handle data";
        return QVariant();
    }
    QVariant headerData(int section, Qt::Orientation orientation,
                                  int role) const override
    {
        if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
            QVariant v;
            v.setValue(static_cast<Column>(section));
            return v.toString();
        }
        return QVariant();
    }
    /*
    add(QCPGraph* graph, const QString& title,
        bool visible, int width, const QColor& color,
        Qt::PenStyle penStyle, QCPGraph::LineStyle lineStyle) {
        QVariant v;
        v.setValue(vGraph(graph, title));
        GraphItem* graphItem = new GraphItem(v);
        GraphItem* visibleItem = new GraphItem(visible);
        GraphItem* widthItem = new GraphItem(width);
        GraphItem* colorItem = new GraphItem(color);
        GraphItem* Item = new GraphItem(color);
                        GraphItem* colorItem = new GraphItem(color);
    }
    data()
    */
};

class SweepDataView: public QCustomPlot
{
    Q_OBJECT

public:
    explicit SweepDataView(QWidget *parent = 0);
    ~SweepDataView();;

public slots:
    // mouse interaction
    void mouseZoom(QWheelEvent*);
    void mouseTrace(QMouseEvent *event);
    void mouseSelect(QMouseEvent *event);

    // canvas controls
    void refresh();
    void syncXAxis2(QCPRange range);
    void save();
    void setZoomX(bool zoomX);
    void setZoomY(bool zoomY);
    void setShowCutStats(bool showCutStats);
    void syncYAxis2(QCPRange range) {
        this->xAxis2->setRange(pow(10, range.lower / 20),
                               pow(10, range.upper / 20));
    }

    // plots
    void plot(const kids::SweepData& sweepdata) {
        namespace eiu = eigen_utils;
        // get data
        SPDLOG_TRACE("freqpack{}", sweepdata.freqpack);
        SPDLOG_TRACE("datapack{}", sweepdata.datapack);

        Eigen::ArrayXd freq = sweepdata.freqpack.array() / 1e6;
        auto fmin = freq.minCoeff() / 1e6;
        auto fmax = freq.maxCoeff() / 1e6;

        Eigen::ArrayXd ampl = sweepdata.datapack.array().abs();
        auto max = ampl.maxCoeff();
        ampl = 20. * (ampl / max).log10();
        SPDLOG_TRACE("freq{}", freq);
        SPDLOG_TRACE("ampl{}", ampl);

        auto toQV = [](const auto m) {
            QVector<double> ret;
            for (int i = 0; i < m.size(); ++i) {
                ret.append(m(i));
            }
            return ret;
        };

        auto x = toQV(Eigen::Map<Eigen::VectorXd>(freq.data(), freq.size()));
        auto y = toQV(Eigen::Map<Eigen::VectorXd>(ampl.data(), ampl.size()));
        qDebug() << "plot sweepdata (" << x.size() << ", " << y.size() << ")";
        this->graph(0)->setData(x, y);
        this->xAxis->setRange(fmin, fmax);
        this->replot();
    };
    void plotFinderResult() {}
    void plotFitterResult() {}

    GraphListModel* graphListModel() {
        return &graphList;
    }

signals:
    // mouse
    void mouseMessage(const QString& msg);

    void savedMessage(QString s, int t);

    void selected(int index);

private:

    // gui states
    bool isShowCutStats;
    bool isZoomX;
    bool isZoomY;

    // plot elements
    QCPItemText *title;
    QCPItemText *stats;
    QCPItemTracer *mouseTracer;
    QCPItemTracer *mouseSelector;
    QMap<QString, QCPGraph*> graphs;

    GraphListModel graphList;
};

namespace solarized
{
inline static const QColor base03  = QColor (  0,  43,  54);
inline static const QColor base02  = QColor (  7,  54,  66);
inline static const QColor base01  = QColor ( 88, 110, 117);
inline static const QColor base00  = QColor (101, 123, 131);
inline static const QColor base0   = QColor (131, 148, 150);
inline static const QColor base1   = QColor (147, 161, 161);
inline static const QColor base2   = QColor (238, 232, 213);
inline static const QColor base3   = QColor (253, 246, 227);
inline static const QColor yellow  = QColor (181, 137,   0);
inline static const QColor orange  = QColor (203,  75,  22);
inline static const QColor red     = QColor (220,  50,  47);
inline static const QColor magenta = QColor (211,  54, 130);
inline static const QColor violet  = QColor (108, 113, 196);
inline static const QColor blue    = QColor ( 38, 139, 210);
inline static const QColor cyan    = QColor ( 42, 161, 152);
inline static const QColor green   = QColor (133, 153,   0);
}   // namespace solarized


