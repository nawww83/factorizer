#include <QGuiApplication>
#include <QQmlApplicationEngine>

int main(int argc, char *argv[])
{
    QGuiApplication app(argc, argv);

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
