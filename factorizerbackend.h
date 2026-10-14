#pragma once

#include <QObject>
#include <QString>
#include <QtQml/qqmlregistration.h>

class FactorizerBackend : public QObject
{
    Q_OBJECT
    QML_ELEMENT // Регистрирует тип в модуле factorizer для использования в QML

    Q_PROPERTY(bool darkMode READ isDarkMode WRITE setDarkMode NOTIFY darkModeChanged)
    Q_PROPERTY(bool isRunning READ isRunning NOTIFY isRunningChanged)
    Q_PROPERTY(QString statusText READ statusText NOTIFY statusTextChanged)
    Q_PROPERTY(QString resultText READ resultText NOTIFY resultTextChanged)
    Q_PROPERTY(bool skipEcm READ skipEcm WRITE setSkipEcm NOTIFY skipEcmChanged)

public:
    explicit FactorizerBackend(QObject *parent = nullptr);
    ~FactorizerBackend() override;

    bool isRunning() const { return m_isRunning; }
    QString statusText() const { return m_statusText; }
    QString resultText() const { return m_resultText; }

    bool isDarkMode() const;
    void setDarkMode(bool val);

    // Пропуск алгоритма ECM при факторизации (полезно для полупростых чисел:
    // сразу запускается решето, методы Полларда продолжают работать)
    bool skipEcm() const { return m_skipEcm; }
    void setSkipEcm(bool val);

    Q_INVOKABLE void startFactorization(const QString &numberStr, int bitWidth);
    Q_INVOKABLE void stop();
    Q_INVOKABLE QString generateRandomNumber(int bitWidth);
    Q_INVOKABLE void generateSemiprime(int bitWidth);
    Q_INVOKABLE int getBitLength(const QString &numberStr);

signals:
    void isRunningChanged();
    void statusTextChanged();
    void resultTextChanged();
    void darkModeChanged();
    void skipEcmChanged();
    void semiprimeReady(const QString &number);

private:
    void setRunning(bool running);
    void setStatus(const QString &status);
    void setResult(const QString &res);

    bool m_isRunning = false;
    bool m_darkMode = true;
    bool m_skipEcm = false;
    std::atomic<bool> m_userStopped{false};
    QString m_statusText = "Готов к работе";
    QString m_resultText;
};
