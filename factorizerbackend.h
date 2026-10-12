#pragma once

#include <QObject>
#include <QString>
#include <QtQml/qqmlregistration.h>

class FactorizerBackend : public QObject
{
    Q_OBJECT
    QML_ELEMENT // Регистрирует тип в модуле factorizer для использования в QML

    Q_PROPERTY(bool isRunning READ isRunning NOTIFY isRunningChanged)
    Q_PROPERTY(QString statusText READ statusText NOTIFY statusTextChanged)
    Q_PROPERTY(QString resultText READ resultText NOTIFY resultTextChanged)

public:
    explicit FactorizerBackend(QObject *parent = nullptr);
    ~FactorizerBackend() override;

    bool isRunning() const { return m_isRunning; }
    QString statusText() const { return m_statusText; }
    QString resultText() const { return m_resultText; }

    Q_INVOKABLE void startFactorization(const QString &numberStr, int bitWidth);
    Q_INVOKABLE void stop();
    Q_INVOKABLE QString generateRandomNumber(int bitWidth);

signals:
    void isRunningChanged();
    void statusTextChanged();
    void resultTextChanged();

private:
    void setRunning(bool running);
    void setStatus(const QString &status);
    void setResult(const QString &res);

    bool m_isRunning = false;
    std::atomic<bool> m_userStopped{false};
    QString m_statusText = "Готов к работе";
    QString m_resultText;
};