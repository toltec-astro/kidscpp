#pragma once

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>

class DoubleRangeSlider: public QWidget
{
    Q_OBJECT
    Q_PROPERTY(double minimum MEMBER mMinimum NOTIFY minimumChanged)
    Q_PROPERTY(double maximum MEMBER mMaximum NOTIFY maximumChanged)
    Q_PROPERTY(double lowerValue MEMBER mLowerValue NOTIFY lowerValueChanged)
    Q_PROPERTY(double upperValue MEMBER mUpperValue READ upperValue WRITE setUpperValue NOTIFY upperValueChanged USER true)

public:
    DoubleRangeSlider(QWidget* aParent = Q_NULLPTR);

    QSize minimumSizeHint() const override;

    double minimum() const;
    double maximum() const;
    double lowerValue() const;
    double upperValue() const;

protected:
    void paintEvent(QPaintEvent* aEvent) override;
    void mousePressEvent(QMouseEvent* aEvent) override;
    void mouseMoveEvent(QMouseEvent* aEvent) override;
    void mouseReleaseEvent(QMouseEvent* aEvent) override;
    void changeEvent(QEvent* aEvent) override;

    QRectF firstHandleRect() const;
    QRectF secondHandleRect() const;
    QRectF barHandleRect() const;
    QRectF handleRect(double aValue, double yoffset=0., double wscale=1., double hscale=1.) const;

signals:
    void lowerValueChanged(double aLowerValue);
    void upperValueChanged(double aUpperValue);
    void valuesChanged(double aLowerValue, double aUpperValue);

    void minimumChanged(double aMinimum);
    void maximumChanged(double aMaximum);
    void rangeChanged(double aMinimum, double aMaximum);

public slots:

    void setLowerValue(double aLowerValue);
    void setUpperValue(double aUpperValue);
    void setValues(double aLowerValue, double aUpperValue);

    void setMinimum(double aMinimum);
    void setMaximum(double aMaximum);
    void setRange(double aMinimum, double aMaximum);

    void setValuesOffset(double offset);

private:
    Q_DISABLE_COPY(DoubleRangeSlider)
    int validWidth() const;
    bool isInverted() const;
    QColor getBackgroundColorEnabled() const;

    double mMinimum;
    double mMaximum;
    double mLowerValue;
    double mUpperValue;
    bool mFirstHandlePressed;
    bool mSecondHandlePressed;
    bool mBarHandlePressed;
    double mInterval;
    double mDelta;
    QColor mBackgroundColorEnabled;
    QColor mBackgroundColorDisabled;
    QColor mBackgroundColorEnabledInverted;
    QColor mBackgroundColor;
};
