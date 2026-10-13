#include "factorizerbackend.h"
#include "generic_factorizer.hpp"
#include "u128_utils.h"
#include <QElapsedTimer>
#include <thread>
#include <sstream>
#include <vector>
#include <map>
#include <optional>
#include <bit>
#include <cstdint>
#include <random>
#include <string>
#include <utility>

namespace {

template <typename T>
struct FactorResult {
    std::map<T, int> factors;       // Сомножитель -> степень
    std::optional<T> remainder;     // Неразложенный остаток (при остановке)
    bool completed = true;
};

template <typename T>
T get_random_value() {
    static thread_local std::mt19937_64 rng(std::random_device{}());

    if constexpr (std::is_same_v<T, ::bignum::u128::U128>) {
        uint64_t lo = rng();
        uint64_t hi = rng() | (1ULL << 63);
        return ::bignum::u128::U128(lo, hi);
    } else if constexpr (requires { typename T::half_type; }) {
        using Half = typename T::half_type;
        Half hi = get_random_value<Half>();
        Half lo = get_random_value<Half>();
        return T(hi, lo);
    } else {
        return T(rng());
    }
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
        res.completed = true;
        res.remainder = std::nullopt;
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

FactorizerBackend::FactorizerBackend(QObject *parent) : QObject(parent) {}

FactorizerBackend::~FactorizerBackend() {
    stop();
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

    std::string s;
    switch (bitWidth) {
    case 128:  s = get_random_value<U128>().toString(); break;
    case 256:  s = get_random_value<U256>().toString(); break;
    case 512:  s = get_random_value<U512>().toString(); break;
    case 1024: s = get_random_value<U1024>().toString(); break;
    case 2048: s = get_random_value<U2048>().toString(); break;
    case 4096: s = get_random_value<U4096>().toString(); break;
    default:   s = get_random_value<U128>().toString(); break;
    }
    return QString::fromStdString(s);
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
    setStatus(QString("Выполняется факторизация (%1 бит)...").arg(val.bitLength));

    m_userStopped.store(false);
    ::u128::Globals::SetStop(false);

    std::string stdStr = numberStr.trimmed().toStdString();

    std::thread([this, stdStr, bitWidth]() {
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
                ss << "\n\n(Факторизация не завершена: методы исчерпаны)\n";
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
            switch (bitWidth) {
            case 128: {
                U128 val = U128::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            case 256: {
                U256 val = U256::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            case 512: {
                U512 val = U512::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            case 1024: {
                U1024 val = U1024::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            case 2048: {
                U2048 val = U2048::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            case 4096: {
                U4096 val = U4096::fromString(stdStr);
                auto raw_factors = factor_parallel_generic(val, update_ui_status, &profiler);
                profiler.finish();
                auto res = build_factor_result(raw_factors, val, m_userStopped.load());
                output = format_result(res, timer.elapsed(), profiler);
                break;
            }
            }
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
