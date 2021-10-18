#include "doublerangeslider.h"
#include <QDebug>

namespace
{

const int scHandleSideLength = 21;
const int scSliderBarHeight = 2;
const int scLeftRightMargin = 11;

}


DoubleRangeSlider::DoubleRangeSlider(QWidget* aParent)
    : QWidget(aParent),
      mMinimum(0.),
      mMaximum(1.),
      mLowerValue(0.),
      mUpperValue(1.),
      mFirstHandlePressed(false),
      mSecondHandlePressed(false),
      mBarHandlePressed(false),
      mInterval(mMaximum - mMinimum),
      mBackgroundColorEnabled(QColor(0x1E, 0x90, 0xFF)),
      mBackgroundColorDisabled(Qt::darkGray),
      mBackgroundColorEnabledInverted(QColor(0xFF, 0x90, 0x1E)),
      mBackgroundColor(mBackgroundColorEnabled)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
}

void DoubleRangeSlider::paintEvent(QPaintEvent* aEvent)
{
    Q_UNUSED(aEvent);
    QPainter painter(this);

    // Background
    QRectF backgroundRect = QRectF(
                scLeftRightMargin,
                (height() - scSliderBarHeight) / 2.,
                validWidth(),
                scSliderBarHeight
                );
    QPen pen(Qt::gray, 0.8);
    painter.setPen(pen);
    painter.setRenderHint(QPainter::Qt4CompatiblePainting);
    QBrush backgroundBrush(QColor(0xD0, 0xD0, 0xD0));
    painter.setBrush(backgroundBrush);
    painter.drawRoundedRect(backgroundRect, 1, 1);

    // First value handle rect  // top
    pen.setColor(Qt::darkGray);
    pen.setWidthF(0.5);
    painter.setPen(pen);
    painter.setRenderHint(QPainter::Antialiasing);
    QBrush handleBrush(QColor(0xFA, 0xFA, 0xFA));
    painter.setBrush(handleBrush);

    QRectF leftHandleRect = firstHandleRect();
    // painter.drawRoundedRect(leftHandleRect, 2, 2);
    QPainterPath leftHandlePath;
    leftHandlePath.moveTo(leftHandleRect.left() + (leftHandleRect.width() / 2), leftHandleRect.bottom());
    leftHandlePath.lineTo(leftHandleRect.topLeft());
    leftHandlePath.lineTo(leftHandleRect.topRight());
    leftHandlePath.lineTo(leftHandleRect.left() + (leftHandleRect.width() / 2), leftHandleRect.bottom());
    painter.strokePath(leftHandlePath, pen);
    painter.fillPath(leftHandlePath, handleBrush);

    // Second value handle rect
    QRectF rightHandleRect = secondHandleRect();
    // painter.drawRoundedRect(rightHandleRect, 2, 2);
    QPainterPath rightHandlePath;
    rightHandlePath.moveTo(rightHandleRect.left() + (rightHandleRect.width() / 2), rightHandleRect.top());
    rightHandlePath.lineTo(rightHandleRect.bottomLeft());
    rightHandlePath.lineTo(rightHandleRect.bottomRight());
    rightHandlePath.lineTo(rightHandleRect.left() + (rightHandleRect.width() / 2), rightHandleRect.top());
    painter.strokePath(rightHandlePath, pen);
    painter.fillPath(rightHandlePath, handleBrush);

    // bar handle
    QRectF centerHandleRect = barHandleRect();
    painter.drawRoundedRect(centerHandleRect, 1, 1);

    // bar
    painter.setRenderHint(QPainter::Antialiasing, false);
    QRectF selectedRect(backgroundRect);
    double l = leftHandleRect.center().x();
    double r = rightHandleRect.center().x();
    if (l > r) std::swap(l, r);
    selectedRect.setLeft(l);
    selectedRect.setRight(r);
    QBrush selectedBrush(getBackgroundColorEnabled());
    // QBrush selectedBrush(QColor(0xFA, 0x00, 0x00));
    painter.setBrush(selectedBrush);
    painter.drawRect(selectedRect);

    // bar handle filled
    QRectF selectedHandleRect = barHandleRect();
    selectedHandleRect.setLeft(
                l > centerHandleRect.left()? l:centerHandleRect.left()
                );
    selectedHandleRect.setRight(
                r < centerHandleRect.right()? r:centerHandleRect.right()
                );
    painter.drawRoundedRect(selectedHandleRect, 1, 1);

}

QRectF DoubleRangeSlider::firstHandleRect() const
{
    double percentage = (mLowerValue - mMinimum) / mInterval;
    return handleRect(percentage * validWidth(), -1., 0.577, 0.5);
}

QRectF DoubleRangeSlider::secondHandleRect() const
{
    double percentage = (mUpperValue - mMinimum) / mInterval;
    return handleRect(percentage * validWidth(), 1., 0.577, 0.5);
}

QRectF DoubleRangeSlider::barHandleRect() const
{
    double center = ((mUpperValue + mLowerValue) / 2. - mMinimum) / mInterval;
    return handleRect(center * validWidth(), 0., 1., 0.333);
}

QRectF DoubleRangeSlider::handleRect(double aValue, double yoffset, double wscale, double hscale) const
{
    double wsize = scHandleSideLength * wscale;
    double hsize = scHandleSideLength * hscale;
    return QRectF(aValue + scLeftRightMargin - wsize / 2., (height()-hsize) / 2. + yoffset * hsize / 2., wsize, hsize);
}

void DoubleRangeSlider::mousePressEvent(QMouseEvent* aEvent)
{
    if(aEvent->buttons() & Qt::LeftButton)
    {
        mBarHandlePressed = barHandleRect().contains(aEvent->pos());
        mSecondHandlePressed = !mBarHandlePressed && secondHandleRect().contains(aEvent->pos());
        mFirstHandlePressed = !mBarHandlePressed && !mSecondHandlePressed && firstHandleRect().contains(aEvent->pos());
        if (mBarHandlePressed)
        {
            mDelta = aEvent->pos().x() - barHandleRect().center().x();
        }
        else if(mSecondHandlePressed)
        {
            mDelta = aEvent->pos().x() - secondHandleRect().center().x();
        }
        else if (mFirstHandlePressed)
        {
            mDelta = aEvent->pos().x() - firstHandleRect().center().x();
        }
        if(aEvent->pos().y() >= 2
           && aEvent->pos().y() <= height() - 2)
        {
            double step = mInterval / 20.;
            double l = isInverted()?secondHandleRect().left():firstHandleRect().left();
            double r = isInverted()?firstHandleRect().right():secondHandleRect().right();
            if (aEvent->pos().x() < l)
            {
                setValuesOffset(-step);
            }
            else if (aEvent->pos().x() > r)
            {
                setValuesOffset(step);
            }
        }
    }
}

void DoubleRangeSlider::mouseMoveEvent(QMouseEvent* aEvent)
{
    if(aEvent->buttons() & Qt::LeftButton)
    {
        if(mFirstHandlePressed)
        {
            setLowerValue((aEvent->pos().x() - mDelta - scLeftRightMargin) * 1.0 / validWidth() * mInterval + mMinimum);
        }
        else if(mSecondHandlePressed)
        {
            setUpperValue((aEvent->pos().x() - mDelta - scLeftRightMargin) * 1.0 / validWidth() * mInterval + mMinimum);
        }
        else if (mBarHandlePressed)
        {
            setValuesOffset((aEvent->pos().x() - mDelta - scLeftRightMargin) * 1.0 / validWidth() * mInterval + mMinimum - (mLowerValue + mUpperValue) / 2.0);
        }
    }
}

void DoubleRangeSlider::mouseReleaseEvent(QMouseEvent* aEvent)
{
    mFirstHandlePressed = false;
    mSecondHandlePressed = false;
    mBarHandlePressed = false;
    QWidget::mouseReleaseEvent(aEvent);
}

void DoubleRangeSlider::changeEvent(QEvent* aEvent)
{
    if(aEvent->type() == QEvent::EnabledChange)
    {
        if(isEnabled())
        {
            mBackgroundColor = getBackgroundColorEnabled();
        }
        else
        {
            mBackgroundColor = mBackgroundColorDisabled;
        }
        update();
    }
}

QSize DoubleRangeSlider::minimumSizeHint() const
{
    return QSize(scLeftRightMargin * 2 + 100, scHandleSideLength);
}

double DoubleRangeSlider::minimum() const
{
    return mMinimum;
}

double DoubleRangeSlider::maximum() const
{
    return mMaximum;
}

double DoubleRangeSlider::lowerValue() const
{
    return mLowerValue;
}

double DoubleRangeSlider::upperValue() const
{
    return mUpperValue;
}

void DoubleRangeSlider::setLowerValue(double aLowerValue)
{
    if(aLowerValue > mMaximum)
    {
        aLowerValue = mMaximum;
    }

    if(aLowerValue < mMinimum)
    {
        aLowerValue = mMinimum;
    }
    mLowerValue = aLowerValue;
    update();
    emit lowerValueChanged(mLowerValue);
    emit valuesChanged(mLowerValue, mUpperValue);
}

void DoubleRangeSlider::setUpperValue(double aUpperValue)
{
    if(aUpperValue > mMaximum)
    {
        aUpperValue = mMaximum;
    }

    if(aUpperValue < mMinimum)
    {
        aUpperValue = mMinimum;
    }

    mUpperValue = aUpperValue;
    update();
    emit upperValueChanged(mUpperValue);
    emit valuesChanged(mLowerValue, mUpperValue);
    // qDebug() << "upper value changed to " << mUpperValue;
}

void DoubleRangeSlider::setValues(double aLowerValue, double aUpperValue)
{
    setLowerValue(aLowerValue);
    setUpperValue(aUpperValue);
}

void DoubleRangeSlider::setMinimum(double aMinimum)
{
    if(aMinimum <= mMaximum)
    {
        mMinimum = aMinimum;
    }
    else
    {
        double oldMax = mMaximum;
        mMinimum = oldMax;
        mMaximum = aMinimum;
    }
    mInterval = mMaximum - mMinimum;
    update();
    // setLowerValue(mLowerValue);
    // setUpperValue(mUpperValue);
    emit minimumChanged(mMinimum);
    emit rangeChanged(mMinimum, mMaximum);
}

void DoubleRangeSlider::setMaximum(double aMaximum)
{
    if(aMaximum >= mMinimum)
    {
        mMaximum = aMaximum;
    }
    else
    {
        double oldMin = mMinimum;
        mMaximum = oldMin;
        mMinimum = aMaximum;
    }
    mInterval = mMaximum - mMinimum;
    update();
    // setLowerValue(mLowerValue);
    // setUpperValue(mUpperValue);
    emit rangeChanged(mMinimum, mMaximum);
    emit maximumChanged(mMaximum);
}

void DoubleRangeSlider::setValuesOffset(double offset)
{
    double l = mLowerValue + offset;
    double u = mUpperValue + offset;
    if (isInverted())
    {
        if (u < mMinimum) offset = mMinimum - mUpperValue;
        else if (l > mMaximum) offset = mMaximum - mLowerValue;
    } else {
        if (l < mMinimum) offset = mMinimum - mLowerValue;
        else if (u > mMaximum) offset = mMaximum - mUpperValue;
    }
    setValues(mLowerValue + offset, mUpperValue + offset);
}

int DoubleRangeSlider::validWidth() const
{
    return width() - scLeftRightMargin * 2;
}

bool DoubleRangeSlider::isInverted() const
{
    return mLowerValue > mUpperValue;
}

QColor DoubleRangeSlider::getBackgroundColorEnabled() const
{
    return isInverted()?mBackgroundColorEnabledInverted:mBackgroundColorEnabled;
}

void DoubleRangeSlider::setRange(double aMinimum, double mMaximum)
{
    setMinimum(aMinimum);
    setMaximum(mMaximum);
}
