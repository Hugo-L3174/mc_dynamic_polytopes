#pragma once
#include <mc_rtc/clock.h>
#include <mc_rtc/gui/StateBuilder.h>
#include <mc_rtc/log/Logger.h>
#include <future>
#include <optional>

// TODO: integrate into mc_rtc
namespace mc_dynamic_polytopes
{

/**
 * @brief Helper base class for asynchronous jobs using CRTP.
 *
 * This class provides a generic interface for running asynchronous computations,
 * managing their state, and adding logging and GUI.
 * Child classes must inherit using CRTP and implement required methods.
 *
 * Upon destruction, any associated GUI elements and log entries are automatically removed.
 *
 * The `input` member may be modified when the async task is not running (i.e., before calling startAsync or when
 * running() is false). The input should not be modified while the async task is running (running() is true).
 *
 * @tparam Derived The child class type (CRTP).
 * @tparam Input The input type for the job.
 * @tparam Result The result type produced by the job.
 *
 * @section RequiredMethods Required/Optional Child Class Methods
 * - Result computeJob()
 *   - Implements the actual computation. Called asynchronously. (Required)
 * - void addToLoggerImpl()
 *   - Adds job-specific log entries after the first result is available. (Optional)
 *     @note Will only be called if the user has previously called addToLogger().
 * - void addToGUIImpl()
 *   - Adds job-specific GUI elements after the first result is available. (Optional)
 *     @note Will only be called if the user has previously called addToGUI().
 *
 * @section UsageExample Example Usage
 * @code
 * struct MyAsyncJob : public MakeAsyncJob<MyAsyncJob, MyInput, MyResult>
 * {
 *   MyResult computeJob()
 *   {
 *     MyResult res;
 *     res.value = input.data * 2;
 *     return res;
 *   }
 *
 *   // Optional: add logging
 *   void addToLoggerImpl()
 *   {
 *     logger_->addLogEntry("my_job_value", this, [this]() { return lastResult_->value; });
 *   }
 *
 *   // Optional: add GUI
 *   void addToGUIImpl()
 *   {
 *     gui_->addElement(this, guiCategory_, mc_rtc::gui::Label("Result", [this]() { return lastResult_->value; }));
 *   }
 * };
 *
 * // Usage:
 * MyAsyncJob job;
 * job.input = ...;
 * job.startAsync();
 * // Later, check for completion:
 * if(job.checkResult())
 * {
 *   auto result = *job.lastResult();
 * }
 * // To enable deferred logger/GUI calls:
 * job.addToLogger(logger, "prefix");
 * job.addToGUI(gui, {"Category"});
 * // Upon destruction, GUI and logging entries are automatically removed.
 * @endcode
 */
template<typename Derived, typename Input, typename Result>
struct AsyncJob
{
  struct Timers
  {
    mc_rtc::duration_ms dt_startAsync = mc_rtc::duration_ms::zero();
  };

public:
  Timers timers;
  Input input;

public:
  AsyncJob() = default;
  AsyncJob(const AsyncJob &) = delete;
  AsyncJob & operator=(const AsyncJob &) = delete;
  AsyncJob(AsyncJob &&) = delete;
  AsyncJob & operator=(AsyncJob &&) = delete;

  ~AsyncJob()
  {
    removeFromLogger();
    removeFromGUI();
  }

  void startAsync()
  {
    auto start_async = mc_rtc::clock::now();
    running_ = true;
    futureResult_ = std::async(
        [this]()
        {
          running_ = true;
          return static_cast<Derived *>(this)->computeJob();
        });
    timers.dt_startAsync = mc_rtc::clock::now() - start_async;
  }

  bool running() const noexcept
  {
    return running_;
  }

  bool checkResult()
  {
    if(futureResult_.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
      lastResult_ = futureResult_.get();
      if(logger_ && !inLogger_)
      {
        static_cast<Derived *>(this)->addToLoggerImpl();
        inLogger_ = true;
      }
      if(gui_ && !inGUI_)
      {
        static_cast<Derived *>(this)->addToGUIImpl();
        inGUI_ = true;
      }
      running_ = false;
      return true;
    }
    return false;
  }

  const std::optional<Result> & lastResult() const noexcept
  {
    return lastResult_;
  }

  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix)
  {
    if(logger_) return;
    logger_ = &logger;
    loggerPrefix_ = prefix;
  }

  void addToGUI(mc_rtc::gui::StateBuilder & gui, const std::vector<std::string> & category)
  {
    if(gui_) return;
    gui_ = &gui;
    guiCategory_ = category;
  }

  void removeFromLogger()
  {
    if(logger_)
    {
      logger_->removeLogEntries(this);
    }
  }

  void removeFromGUI()
  {
    if(gui_)
    {
      gui_->removeElements(this);
    }
  }

protected: // bookkeeping for the async job
  bool running_ = false;
  bool inLogger_ = false;
  bool inGUI_ = false;
  std::future<Result> futureResult_;
  std::optional<Result> lastResult_;

  // CRTP: derived must implement this
  // Result computeJob();

  // CRTP: derived must implement these
  void addToLoggerImpl() {}
  void addToGUIImpl() {}

  // Store context for logging/GUI
  mc_rtc::Logger * logger_ = nullptr;
  std::string loggerPrefix_ = "";
  mc_rtc::gui::StateBuilder * gui_ = nullptr;
  std::vector<std::string> guiCategory_;
};

// Helper alias for CRTP inheritance
template<typename Derived, typename Input, typename Result>
using MakeAsyncJob = mc_dynamic_polytopes::AsyncJob<Derived, Input, Result>;
} // namespace mc_dynamic_polytopes
