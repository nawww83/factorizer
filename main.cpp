#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQuickStyle>

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);

    // Нативный стиль Windows не позволяет переопределять background/contentItem
    // в QML, поэтому явно задаём кастомизируемый стиль Basic.
    QQuickStyle::setStyle(QStringLiteral("Basic"));

    app.setOrganizationName("FactorizerApp");
    app.setOrganizationDomain("factorizer.local");
    app.setApplicationName("ParallelFactorizer");

    QQmlApplicationEngine engine;
    QObject::connect(
        &engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []() { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);
    engine.loadFromModule("factorizer", "Main");

    return QGuiApplication::exec();
}
