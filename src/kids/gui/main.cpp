#include "kids_gui.h"
#include <QApplication>
#include <QDebug>
#include <QtGlobal>

int main(int argc, char *argv[]) {
    // setup logging format
    QString pattern =
        "[%{if-debug}D%{endif}%{if-info}I%{endif}%{if-warning}W%{endif}%{if-"
        "critical}C%{endif}%{if-fatal}F%{endif}] %{message}";
#ifdef QT_DEBUG
    pattern += " (%{function})";
#endif
    qSetMessagePattern(pattern);
    qDebug() << "QT version:" << QT_VERSION_STR
             << "; CXX version:" << __cplusplus;

    // create windows
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
