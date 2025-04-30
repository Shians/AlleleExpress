#pragma once

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <memory>
#include <atomic>
#include <type_traits>
#include <concepts> // C++20

constexpr int DEFAULT_N_THREADS = 2; // Default number of threads

/**
 * @class ThreadPool
 * @brief A modern thread pool implementation for parallel processing using C++20 features
 */
class ThreadPool {
public:
    // Constructor creates the thread pool with specified number of threads
    explicit ThreadPool(size_t num_threads = DEFAULT_N_THREADS) {
        start(num_threads);
    }

    // Destructor stops all threads and joins them
    ~ThreadPool() {
        stop();
    }

    // Enqueue a task to be executed by a thread in the pool
    template<typename F, typename... Args>
        requires std::invocable<F, Args...>
    [[nodiscard]] auto enqueue(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        using return_type = std::invoke_result_t<F, Args...>;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            [func = std::forward<F>(f), ...args = std::forward<Args>(args)]() -> return_type {
                return func(args...);
            }
        );

        std::future<return_type> result = task->get_future();
        {
            std::scoped_lock lock(queue_mutex_);

            // Don't allow enqueueing after stopping the pool
            if (stop_requested_) {
                throw std::runtime_error("Enqueue on stopped ThreadPool");
            }

            tasks_.emplace([task](){ (*task)(); });
        }
        condition_.notify_one();
        return result;
    }

    // Get the number of active threads in the pool
    [[nodiscard]] size_t size() const noexcept {
        return threads_.size();
    }

    // Get the number of tasks waiting in the queue
    [[nodiscard]] size_t queue_size() const {
        std::scoped_lock lock(queue_mutex_);
        return tasks_.size();
    }

private:
    // Start the thread pool with specified number of threads
    void start(size_t num_threads) {
        stop_requested_ = false;
        threads_.reserve(num_threads);

        for (size_t i = 0; i < num_threads; ++i) {
            threads_.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock lock(queue_mutex_);

                        // Wait for task or stop signal
                        condition_.wait(lock, [this] {
                            return stop_requested_ || !tasks_.empty();
                        });

                        // Check if we need to exit due to stop request
                        if (stop_requested_ && tasks_.empty()) {
                            return;
                        }

                        // Get next task
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }

                    // Execute task
                    task();
                }
            });
        }
    }

    // Stop the thread pool
    void stop() {
        {
            std::scoped_lock lock(queue_mutex_);
            stop_requested_ = true;
        }
        condition_.notify_all();

        // Join all threads
        for (auto& thread : threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        // Clear threads vector
        threads_.clear();
    }

    // Thread pool state using standard thread instead of jthread
    std::vector<std::thread> threads_;
    std::queue<std::function<void()>> tasks_;

    // Synchronization
    mutable std::mutex queue_mutex_;
    std::condition_variable condition_;
    bool stop_requested_ = false;
};
