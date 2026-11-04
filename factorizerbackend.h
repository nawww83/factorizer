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
    // Предупреждение о неприменимости выбранного метода (например, SIQS на
    // 512 битах). В отличие от statusText, не затирается прогрессом.
    Q_PROPERTY(QString warningText READ warningText NOTIFY warningTextChanged)
    Q_PROPERTY(QString resultText READ resultText NOTIFY resultTextChanged)
    // Доля выполненной работы 0..1 для индикатора прогресса. Монотонно
    // неубывающая внутри одного запуска (см. setProgress).
    Q_PROPERTY(double progress READ progress NOTIFY progressChanged)
    // Способ факторизации (решето): 0 — авто, 1 — принудительно SIQS.
    Q_PROPERTY(int method READ method WRITE setMethod NOTIFY methodChanged)

public:
    explicit FactorizerBackend(QObject *parent = nullptr);
    ~FactorizerBackend() override;

    bool isRunning() const { return m_isRunning; }
    QString statusText() const { return m_statusText; }
    QString warningText() const { return m_warningText; }
    QString resultText() const { return m_resultText; }
    double progress() const { return m_progress; }

    bool isDarkMode() const;
    void setDarkMode(bool val);

    // Способ факторизации: 0 — авто (решето SIQS до 320 бит, ECM-лестница только
    // выше этого предела), 1 — принудительно SIQS. Явный выбор снимает ограничение
    // по разрядности: пользователь сам решает, что готов ждать (жёсткий предел —
    // 512 бит).
    int method() const { return m_method; }
    void setMethod(int val);

    Q_INVOKABLE void startFactorization(const QString &numberStr, int bitWidth);
    Q_INVOKABLE void stop();
    Q_INVOKABLE QString generateRandomNumber(int bitWidth);
    Q_INVOKABLE void generateSemiprime(int bitWidth);
    Q_INVOKABLE int getBitLength(const QString &numberStr);

signals:
    void isRunningChanged();
    void statusTextChanged();
    void warningTextChanged();
    void resultTextChanged();
    void progressChanged();
    void darkModeChanged();
    void methodChanged();
    void semiprimeReady(const QString &number);

private:
    void setRunning(bool running);
    void setStatus(const QString &status);
    void setWarning(const QString &warning);
    void setResult(const QString &res);
    // Монотонное обновление прогресса (0..1): назад индикатор не идёт.
    void setProgress(double p);

    bool m_isRunning = false;
    bool m_darkMode = true;
    int m_method = 0; // 0 — авто, 1 — принудительно SIQS
    std::atomic<bool> m_userStopped{false};
    QString m_statusText = "Готов к работе";
    QString m_warningText;
    QString m_resultText;
    double m_progress = 0.0;
};
