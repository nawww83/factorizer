#pragma once

// ---------------------------------------------------------------------------
// siqs_xor.hpp — самодостаточное квадратичное решето (SIQS) с 2LP-циклами.
//
// Портировано из tests/siqs_xor.cpp: одиночный проход по A-полиномам,
// L1/L2-сито по байтовым логарифмам, hit-map для больших простых, лес
// циклов 2LP и разреженный гаусс над GF(2) с XOR-историей строк.
//
// Отличия от эталонного файла:
//   * нет OpenMP — батчи A-полиномов считаются std::thread-воркерами
//     (используется cpu_topology, как в ECM);
//   * весь вывод идёт через status_cb (std::cout в приложении не нужен);
//   * участвует глобальный флаг остановки u128::Globals::LoadStop();
//   * вместо собственного U512 используется фиксированный bignum::UInt<512>.
//
// Код разнесён по частям в siqs/; этот файл — только точка сборки, поэтому
// внешний интерфейс (siqs_xor::factor<T>) не менялся:
//
//   common.hpp            типы, флаг остановки, модулярная арифметика
//   relation.hpp          FactorExp / RawRelation / FullRelation, слияние
//   worker_pool.hpp       SieveStats, BatchResult, постоянный пул потоков
//   large_prime_graph.hpp граф циклов 2LP на остовном лесу
//   gauss.hpp             разреженный гаусс над GF(2)
//   config.hpp            SIQSConfig — параметры одного запуска
//   poly_batch.hpp        update_roots + сито одного A-полинома
//   driver.hpp            автоподбор параметров, run_siqs, factor<T>
// ---------------------------------------------------------------------------

#include "siqs/common.hpp"
#include "siqs/relation.hpp"
#include "siqs/worker_pool.hpp"
#include "siqs/large_prime_graph.hpp"
#include "siqs/gauss.hpp"
#include "siqs/config.hpp"
#include "siqs/poly_batch.hpp"
#include "siqs/driver.hpp"
