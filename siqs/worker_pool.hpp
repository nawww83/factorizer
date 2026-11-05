#pragma once

// ---------------------------------------------------------------------------
// siqs/worker_pool.hpp — статистика решета и постоянный пул воркеров.
//
// SieveStats считает полиномы, выживших, состав соотношений и время фаз;
// BatchResult — порция соотношений вместе с её статистикой. WorkerPool
// держит потоки всё решето: батч — это job, который получает каждый воркер
// вместе со своим номером (номер задаёт seed A-набора).
// ---------------------------------------------------------------------------

#include "common.hpp"
#include "relation.hpp"

namespace bignum::siqs_xor {

// --- Статистика и пул воркеров ---------------------------------------------

// Счётчики и таймеры одной порции полиномов. Собираются всегда: восемь
// вызовов steady_clock::now() на полином — доли процента от его стоимости
// (на 300 битах полином стоит ~275 мкс), зато любая оптимизация проверяется
// по одним и тем же числам, а не по ощущениям.
struct SieveStats
{
    u64 polynomials = 0;
    u64 survivors = 0;
    u64 large_prime_hits = 0;
    u64 hits_in_survivors = 0;
    u64 rel_full = 0;
    u64 rel_1lp = 0;
    u64 rel_2lp = 0;
    u64 rel_3lp = 0;
    u64 dropped_big_cofactor = 0;
    u64 dropped_nonsplit = 0;
    double t_prep = 0.0;
    double t_sieve = 0.0;
    double t_scan = 0.0;
    double t_factor = 0.0;

    void add(const SieveStats &o)
    {
        polynomials += o.polynomials;
        survivors += o.survivors;
        large_prime_hits += o.large_prime_hits;
        hits_in_survivors += o.hits_in_survivors;
        rel_full += o.rel_full;
        rel_1lp += o.rel_1lp;
        rel_2lp += o.rel_2lp;
        rel_3lp += o.rel_3lp;
        dropped_big_cofactor += o.dropped_big_cofactor;
        dropped_nonsplit += o.dropped_nonsplit;
        t_prep += o.t_prep;
        t_sieve += o.t_sieve;
        t_scan += o.t_scan;
        t_factor += o.t_factor;
    }

    u64 relations() const { return rel_full + rel_1lp + rel_2lp + rel_3lp; }
};

struct BatchResult
{
    std::vector<RawRelation> relations;
    SieveStats stats;
};

// Постоянный пул воркеров. Раньше std::thread-ы создавались и уничтожались на
// каждый батч (на 300 битах — тысячи запусков за проход), а привязка к P-ядрам
// повторялась каждый раз. Здесь потоки живут всё решето: батч — это job,
// который получает каждый воркер вместе со своим номером. Номер задаёт seed,
// поэтому распределение A-наборов и результат остаются детерминированными.
class WorkerPool
{
public:
    explicit WorkerPool(unsigned n)
    {
        workers_.reserve(n);
        for (unsigned t = 0; t < n; ++t)
            workers_.emplace_back([this, t]() { worker_loop(t); });
    }

    ~WorkerPool() { shutdown(); }

    WorkerPool(const WorkerPool &) = delete;
    WorkerPool &operator=(const WorkerPool &) = delete;

    unsigned size() const { return static_cast<unsigned>(workers_.size()); }

    void run(const std::function<void(unsigned)> &job)
    {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            job_ = &job;
            ++generation_;
            pending_ = workers_.size();
        }
        wake_.notify_all();

        std::unique_lock<std::mutex> lk(mutex_);
        done_.wait(lk, [this]() { return pending_ == 0; });
        job_ = nullptr;
    }

private:
    void worker_loop(unsigned t)
    {
        cpu::pin_current_thread(t);
        u64 seen = 0;
        while (true)
        {
            const std::function<void(unsigned)> *job = nullptr;
            {
                std::unique_lock<std::mutex> lk(mutex_);
                wake_.wait(lk, [this, seen]() { return stopped_ || generation_ != seen; });
                if (stopped_)
                    return;
                seen = generation_;
                job = job_;
            }
            (*job)(t);
            {
                std::lock_guard<std::mutex> lk(mutex_);
                if (--pending_ == 0)
                    done_.notify_one();
            }
        }
    }

    void shutdown()
    {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            stopped_ = true;
        }
        wake_.notify_all();
        for (auto &w : workers_)
            if (w.joinable())
                w.join();
    }

    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::condition_variable done_;
    const std::function<void(unsigned)> *job_ = nullptr;
    u64 generation_ = 0;
    unsigned pending_ = 0;
    bool stopped_ = false;
};

} // namespace bignum::siqs_xor
