#include "factorizerbackend.h"
#include <QElapsedTimer>
#include <QSettings>
#include "generic_factorizer.hpp"
#include "u128_utils.h"
#include <bit>
#include <cstdint>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

template <typename T>
struct FactorResult {
    std::map<T, int> factors;       // Сомножитель -> степень
    std::optional<T> remainder;     // Неразложенный остаток (при остановке)
    bool completed = true;
};

// Случайное число полной ширины типа T (старший бит всегда установлен).
template <typename T>
T get_random_full() {
    static thread_local std::mt19937_64 rng(std::random_device{}());

    if constexpr (std::is_same_v<T, ::bignum::u128::U128>) {
        uint64_t lo = rng();
        uint64_t hi = rng() | (1ULL << 63);
        return ::bignum::u128::U128(lo, hi);
    } else if constexpr (requires { typename T::half_type; }) {
        using Half = typename T::half_type;
        Half hi = get_random_full<Half>();
        Half lo = get_random_full<Half>();
        return T(hi, lo);
    } else {
        return T(rng());
    }
}

// Случайное число с ровно bits битами (bits <= 0 — полная ширина типа T).
// UI позволяет произвольную разрядность 64..512, поэтому полное число
// усекается до нужного размера, а старший бит ставится заново —
// иначе разрядность результата не совпала бы с выбранной.
template <typename T>
T get_random_value(int bits = 0) {
    T v = get_random_full<T>();
    const unsigned full = v.bit_width(); // == sizeof(T)*8
    if (bits <= 0 || static_cast<unsigned>(bits) >= full)
        return v;
    v >>= (full - static_cast<unsigned>(bits));
    v |= T(1) << (bits - 1);
    return v;
}

// Конвертер списка множителей std::vector<T> в FactorResult<T> с подсчетом степеней
template <typename T>
FactorResult<T> build_factor_result(const std::vector<T>& raw_factors, const T& original_val, bool user_stopped) {
    FactorResult<T> res;
    for (const auto& f : raw_factors) {
        res.factors[f]++;
    }

    if (user_stopped) {
        res.completed = false;
        if (res.factors.size() == 1 && res.factors.count(original_val)) {
            res.factors.clear();
            res.remainder = original_val;
        } else {
            T prod = T(1);
            for (const auto& [f, p] : res.factors) {
                for (int i = 0; i < p; ++i) prod = prod * f;
            }
            if (prod > T(0) && prod < original_val) {
                res.remainder = ::bignum::siqs::div_q(original_val, prod);
            }
        }
    } else {
        // Конвейер возвращает исходное число без разложения и в том случае,
        // когда все методы исчерпаны (например, полупростое 512 бит: SIQS
        // ограничен 320 битами, GNFS этой реализации — 150). Раньше такой
        // результат помечался completed и выглядел в UI как «простой множитель».
        // Теперь: если число составное — это честный остаток с пометкой,
        // а если простое — разложение тривиально и всё корректно.
        const bool single = res.factors.size() == 1 && res.factors.count(original_val) == 1;
        if (single && !::bignum::generic_ecm::is_prime_generic(original_val, 25)) {
            res.completed = false;
            res.factors.clear();
            res.remainder = original_val;
        } else {
            res.completed = true;
            res.remainder = std::nullopt;
        }
    }

    return res;
}

struct ValidationResult {
    bool ok = false;
    size_t bitLength = 0;
    QString errorMessage;
};

ValidationResult validateInputNumber(const QString &str, int bitWidth) {
    QString trimmed = str.trimmed();
    if (trimmed.isEmpty()) {
        return {false, 0, "Ошибка: поле ввода пустое."};
    }

    int startIdx = 0;
    if (trimmed[0] == '+') {
        startIdx = 1;
    } else if (trimmed[0] == '-') {
        return {false, 0, "Ошибка: отрицательные числа не поддерживаются."};
    }

    if (startIdx >= trimmed.length()) {
        return {false, 0, "Ошибка: некорректный формат числа."};
    }

    for (int i = startIdx; i < trimmed.length(); ++i) {
        if (!trimmed[i].isDigit()) {
            return {false, 0, QString("Ошибка: некорректный символ '%1'. Допустимы только цифры.").arg(trimmed[i])};
        }
    }

    while (startIdx < trimmed.length() && trimmed[startIdx] == '0') {
        startIdx++;
    }

    if (startIdx == trimmed.length()) {
        return {false, 0, "Число 0 не подлежит факторизации."};
    }

    QString numOnly = trimmed.mid(startIdx);
    if (numOnly == "1") {
        return {false, 1, "Число 1 не является ни простым, ни составным."};
    }

    std::vector<uint32_t> limbs;
    limbs.push_back(0);

    for (QChar qc : std::as_const(numOnly)) {
        uint64_t carry = qc.unicode() - '0';
        for (size_t i = 0; i < limbs.size(); ++i) {
            uint64_t cur = static_cast<uint64_t>(limbs[i]) * 10 + carry;
            limbs[i] = static_cast<uint32_t>(cur & 0xFFFFFFFF);
            carry = cur >> 32;
        }
        if (carry > 0) {
            limbs.push_back(static_cast<uint32_t>(carry));
        }
    }

    size_t bitLength = (limbs.size() - 1) * 32 + (32 - std::countl_zero(limbs.back()));

    if (static_cast<int>(bitLength) > bitWidth) {
        int recommended = 0;
        if (bitLength <= 256) recommended = 256;
        else if (bitLength <= 512) recommended = 512;
        else if (bitLength <= 1024) recommended = 1024;
        else if (bitLength <= 2048) recommended = 2048;
        else if (bitLength <= 4096) recommended = 4096;

        QString msg = QString("Ошибка: число требует %1 бит, что превышает выбранную настройку (%2 бит).")
                          .arg(bitLength)
                          .arg(bitWidth);

        if (recommended > 0) {
            msg += QString("\nРекомендуется переключить разрядность на %1 бит.").arg(recommended);
        } else {
            msg += "\nЧисло превышает максимальную поддерживаемую разрядность (4096 бит).";
        }
        return {false, bitLength, msg};
    }

    return {true, bitLength, ""};
}

} // anonymous namespace

FactorizerBackend::FactorizerBackend(QObject *parent)
    : QObject(parent)
{
    QSettings s;
    m_darkMode = s.value("UI/darkMode", true).toBool(); // по умолчанию true
    m_skipEcm = s.value("Algorithms/skipEcm", false).toBool(); // по умолчанию false
    m_gnfsPreferred = s.value("Algorithms/gnfsPreferred", false).toBool(); // по умолчанию false
}

FactorizerBackend::~FactorizerBackend() {
    stop();
}

bool FactorizerBackend::isDarkMode() const
{
    return m_darkMode;
}

void FactorizerBackend::setDarkMode(bool val)
{
    if (m_darkMode != val) {
        m_darkMode = val;
        QSettings s;
        s.setValue("UI/darkMode", val);
        emit darkModeChanged();
    }
}

void FactorizerBackend::setSkipEcm(bool val)
{
    if (m_skipEcm != val) {
        m_skipEcm = val;
        QSettings s;
        s.setValue("Algorithms/skipEcm", val);
        emit skipEcmChanged();
    }
}

void FactorizerBackend::setGnfsPreferred(bool val)
{
    if (m_gnfsPreferred != val) {
        m_gnfsPreferred = val;
        QSettings s;
        s.setValue("Algorithms/gnfsPreferred", val);
        emit gnfsPreferredChanged();
    }
}

void FactorizerBackend::setRunning(bool r) {
    if (m_isRunning != r) {
        m_isRunning = r;
        emit isRunningChanged();
    }
}

void FactorizerBackend::setStatus(const QString &s) {
    m_statusText = s;
    emit statusTextChanged();
}

void FactorizerBackend::setResult(const QString &r) {
    m_resultText = r;
    emit resultTextChanged();
}

QString FactorizerBackend::generateRandomNumber(int bitWidth) {
    using namespace bignum::generic_ecm;

    // Тип подбираем по диапазону: UI разрешает произвольную разрядность
    // 64..512, поэтому точный switch (128/256/512/...) больше не подходит.
    // Число генерируется ровно такой битности, какая выбрана.
    std::string s;
    if (bitWidth <= 128)
        s = get_random_value<U128>(bitWidth).toString();
    else if (bitWidth <= 256)
        s = get_random_value<U256>(bitWidth).toString();
    else if (bitWidth <= 512)
        s = get_random_value<U512>(bitWidth).toString();
    else if (bitWidth <= 1024)
        s = get_random_value<U1024>(bitWidth).toString();
    else if (bitWidth <= 2048)
        s = get_random_value<U2048>(bitWidth).toString();
    else if (bitWidth <= 4096)
        s = get_random_value<U4096>(bitWidth).toString();
    else
        s = "0";
    return QString::fromStdString(s);
}

void FactorizerBackend::generateSemiprime(int bitWidth) {
    if (m_isRunning) return;

    setRunning(true);
    setStatus(QString("Генерация полупростого числа: p×q (%1 бит)...").arg(bitWidth));

    m_userStopped.store(false);
    ::u128::Globals::SetStop(false);

    std::thread([this, bitWidth]() {
        using namespace bignum::generic_ecm;

        auto progress = [this, bitWidth](uint64_t candidates) {
            setStatus(QString("Генерация полупростого числа p×q (%1 бит): проверено кандидатов: %2...")
                          .arg(bitWidth)
                          .arg(candidates));
        };

        // Тип подбираем по диапазону, а запрошенная битность передаётся
        // генератору как есть — он умеет произвольные значения, включая
        // нечётные (65, 137, ...).
        std::string s;
        if (bitWidth <= 128)
            s = generate_semiprime<U128>(bitWidth, nullptr, nullptr, progress).toString();
        else if (bitWidth <= 256)
            s = generate_semiprime<U256>(bitWidth, nullptr, nullptr, progress).toString();
        else if (bitWidth <= 512)
            s = generate_semiprime<U512>(bitWidth, nullptr, nullptr, progress).toString();
        else if (bitWidth <= 1024)
            s = generate_semiprime<U1024>(bitWidth, nullptr, nullptr, progress).toString();
        else if (bitWidth <= 2048)
            s = generate_semiprime<U2048>(bitWidth, nullptr, nullptr, progress).toString();
        else if (bitWidth <= 4096)
            s = generate_semiprime<U4096>(bitWidth, nullptr, nullptr, progress).toString();
        else
            s = "0";

        bool cancelled = ::u128::Globals::LoadStop();
        setRunning(false);

        if (!cancelled && s != "0") {
            setStatus("Генерация завершена: полупростое число p×q готово");
            emit semiprimeReady(QString::fromStdString(s));
        } else if (cancelled) {
            setStatus("Генерация остановлена пользователем");
        } else {
            setStatus("Генерация не удалась");
        }
    }).detach();
}

int FactorizerBackend::getBitLength(const QString &numberStr)
{
    // Используем существующую функцию валидации и подсчета битов (лимит 4096 бит)
    ValidationResult val = validateInputNumber(numberStr, 4096);
    if (!val.ok && val.bitLength == 0) {
        return 0;
    }
    return static_cast<int>(val.bitLength);
}

void FactorizerBackend::stop() {
    m_userStopped.store(true);
    ::u128::Globals::SetStop(true);
    setStatus("Остановка пользователем...");
}

void FactorizerBackend::startFactorization(const QString &numberStr, int bitWidth) {
    if (m_isRunning) return;

    ValidationResult val = validateInputNumber(numberStr, bitWidth);
    if (!val.ok) {
        setResult(val.errorMessage);
        setStatus("Ошибка ввода");
        return;
    }

    setRunning(true);
    setResult("");

    const bool skipEcm = m_skipEcm;
    const bool gnfsPreferred = m_gnfsPreferred;

    // Сразу проговариваем в строке статуса, будет ли вообще применён ECM:
    // он работает только для чисел от 136 бит (меньшие мгновенно раскладывает
    // быстрый SIQS) и целиком отключается переключателем «Пропускать ECM».
    QString ecmNote;
    if (skipEcm)
        ecmNote = ", ECM отключён (Пропускать ECM)";
    else if (val.bitLength <= 135)
        ecmNote = ", ECM не применяется (≤135 бит — быстрый SIQS)";

    setStatus(QString("Выполняется факторизация (%1 бит)%2...")
                  .arg(val.bitLength)
                  .arg(ecmNote));

    m_userStopped.store(false);
    ::u128::Globals::SetStop(false);

    std::string stdStr = numberStr.trimmed().toStdString();
    const int numberBits = static_cast<int>(val.bitLength);

    std::thread([this, stdStr, numberBits, skipEcm, gnfsPreferred]() {
        using namespace bignum::generic_ecm;
        QElapsedTimer timer;
        timer.start();

        bignum::StageProfiler profiler;
        std::string output = "";

        auto format_result = [this](const auto& factor_res, qint64 elapsedMs, const bignum::StageProfiler& prof) {
            std::stringstream ss;
            bool first = true;

            for (const auto& [factor, power] : factor_res.factors) {
                if (!first) ss << " × ";
                ss << factor.toString();
                if (power > 1) ss << "^" << power;
                first = false;
            }

            using T = std::decay_t<decltype(factor_res.factors.begin()->first)>;
            if (factor_res.remainder.has_value() && *factor_res.remainder > T(1)) {
                if (!first) ss << " × ";
                ss << "[остаток: " << factor_res.remainder->toString() << "]";
            }

            if (m_userStopped.load()) {
                ss << "\n\n(Факторизация остановлена пользователем)\n";
            } else if (!factor_res.completed) {
                ss << "\n\n(Факторизация не завершена: методы исчерпаны";
                // Для чисел крупнее практического предела SIQS объясняем причину:
                // в программе нет метода такой мощности, нужен, например, CADO-NFS.
                T unfactored(0);
                if (factor_res.remainder.has_value())
                    unfactored = *factor_res.remainder;
                else if (!factor_res.factors.empty())
                    unfactored = factor_res.factors.begin()->first;
                if (unfactored.bit_width() > ::bignum::generic_ecm::kSiqsMaxBits)
                    ss << ": число крупнее " << ::bignum::generic_ecm::kSiqsMaxBits
                       << " бит — для таких размеров в программе нет достаточно "
                          "мощного метода (нужен, например, CADO-NFS)";
                ss << ")\n";
            }

            ss << "\nОбщее время: " << (elapsedMs / 1000.0) << " сек.\n";
            ss << prof.format_report();
            return ss.str();
        };

        auto update_ui_status = [this](const std::string& msg) {
            if (!m_userStopped.load()) {
                QMetaObject::invokeMethod(this, [this, msg]() {
                    setStatus(QString::fromStdString(msg));
                });
            }
        };

        try {
            // Тип бигинта выбираем по фактической битности числа, а не по
            // выбранной в UI разрядности: теперь она произвольная (64..512),
            // а число должно помещаться в тип с запасом на промежуточные
            // значения методов факторизации.
            auto run = [&](auto tag) {
                using T = typename decltype(tag)::type;
                T val = T::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler, skipEcm, gnfsPreferred);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
            };

            if (numberBits <= 128)
                run(std::type_identity<U128>{});
            else if (numberBits <= 256)
                run(std::type_identity<U256>{});
            else if (numberBits <= 512)
                run(std::type_identity<U512>{});
            else if (numberBits <= 1024)
                run(std::type_identity<U1024>{});
            else if (numberBits <= 2048)
                run(std::type_identity<U2048>{});
            else if (numberBits <= 4096)
                run(std::type_identity<U4096>{});
            else
                output = "Число больше 4096 бит: нет подходящего типа!";
        } catch (...) {
            output = "Ошибка при вычислении!";
        }

        QMetaObject::invokeMethod(this, [this, output]() {
            setResult(QString::fromStdString(output));
            setStatus(m_userStopped.load() ? "Остановлено пользователем" : "Завершено");
            setRunning(false);
        });
    }).detach();
}
