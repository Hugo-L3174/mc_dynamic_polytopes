#pragma once
#include <mc_rtc/clock.h>
#include <utility>

namespace mc_rtc
{
/**
 * @brief Returns the exectution time of the provided function (in chrono
 * type system)
 *
 * @param func Function to evaluate
 * @param args Arguments to the function
 *
 * @return Time elapsed
 */
template<typename F, typename... Args>
static mc_rtc::duration_ms timedExecution(F && func, Args &&... args)
{
  auto start = mc_rtc::clock::now();
  func(std::forward<Args>(args)...);
  return std::chrono::duration_cast<mc_rtc::duration_ms>(mc_rtc::clock::now() - start);
}

inline double duration_ms_count(std::chrono::time_point<mc_rtc::clock> start,
                                std::chrono::time_point<mc_rtc::clock> end)
{
  return mc_rtc::duration_ms(end - start).count();
}

inline mc_rtc::duration_ms elapsed_ms(std::chrono::time_point<mc_rtc::clock> start)
{
  return mc_rtc::duration_ms(mc_rtc::clock::now() - start);
}

inline double elapsed_ms_count(std::chrono::time_point<mc_rtc::clock> start)
{
  return mc_rtc::elapsed_ms(start).count();
}

} // namespace mc_rtc
