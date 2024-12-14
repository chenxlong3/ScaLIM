#ifndef _TIMER_H
#define _TIMER_H
#include "headers.h"
class Timer
{
private:
    std::chrono::high_resolution_clock::time_point __StartTime, __LastTime, __EndTime;
    const char* __processName;
public:
    Timer(/* args */);
    explicit Timer(const char* processName);
    ~Timer();

    void refresh_time();
    void record_current_time();
    double get_operation_time();
    void log_operation_time();
    void log_operation_time(const char* operation_name);
    void log_operation_time(const char* operation_name, ofstream& ofs);
    // Get the total time from the beginning
	double get_total_time();
    void log_total_time();
    void log_total_time(ofstream& ofs);
    void log_till_now(const char* operation_name);
    void log_till_now(const char* operation_name, ofstream& ofs);
    void log_sub_total_time() const;
};

Timer::Timer(/* args */)
{
    this->__StartTime = std::chrono::high_resolution_clock::now();
    this->__LastTime = this->__StartTime, this->__EndTime = this->__StartTime;
    this->__processName = "Unnamed";
}

Timer::Timer(const char* processName) {
    this->__StartTime = std::chrono::high_resolution_clock::now();
    this->__LastTime = this->__StartTime; this->__EndTime = this->__StartTime;
    this->__processName = processName;
}

Timer::~Timer()
{
}

void Timer::refresh_time() {
    this->__StartTime = std::chrono::high_resolution_clock::now();
    this->__LastTime = this->__StartTime; this->__EndTime = this->__StartTime;
    return;
}

void Timer::record_current_time()
{
    this->__LastTime = this->__EndTime;
    this->__EndTime = std::chrono::high_resolution_clock::now();
}

double Timer::get_operation_time() {
    this->record_current_time();
    std::chrono::duration<double> elapsed = this->__EndTime - this->__LastTime;
    return elapsed.count();
}

void Timer::log_operation_time() {
    const double duration = this->get_operation_time();
    std::cout << "->Time used (sec): " << duration << '\n';
}

void Timer::log_operation_time(const char* operation_name) {
    const double duration = this->get_operation_time();
    std::cout << "->Time used (sec) for operation [" << operation_name << "]: " << duration << '\n';
}

void Timer::log_operation_time(const char* operation_name, ofstream& ofs) {
    const double duration = this->get_operation_time();
    std::cout << "->Time used (sec) for operation [" << operation_name << "]: " << duration << '\n';
    ofs << "->Time used (sec) for operation [" << operation_name << "]: " << duration << '\n';
}

double Timer::get_total_time() {
    this->record_current_time();
    std::chrono::duration<double> elapsed = this->__EndTime - this->__StartTime;
    return elapsed.count();
}

void Timer::log_total_time() {
    const double duration = this->get_total_time();
    std::cout << "--->Time used (sec) for process [" << __processName << "]: " << duration << '\n';
    return;
}

void Timer::log_till_now(const char* operation_name) {
    const double duration = this->get_total_time();
    std::cout << "--->Time used (sec) until operation [" << operation_name << "]: " << duration << '\n';
    return;
}

void Timer::log_till_now(const char* operation_name, ofstream& ofs) {
    const double duration = this->get_total_time();
    std::cout << "--->Time used (sec) until operation [" << operation_name << "]: " << duration << '\n';
    ofs << "--->Time used (sec) until operation [" << operation_name << "]: " << duration << '\n';
    return;
}

void Timer::log_total_time(ofstream& ofs) {
    const double duration = this->get_total_time();
    std::cout << "--->Time used (sec) for process [" << __processName << "]: " << duration << '\n';
    ofs << "--->Time used (sec) for process [" << __processName << "]: " << duration << '\n';
    return;
}

void Timer::log_sub_total_time() const {
    std::chrono::duration<double> elapsed = this->__EndTime - this->__StartTime;
    std::cout << "--->Time used (sec) for process [" << __processName << "]: " << elapsed.count() << '\n';
    return;
}
#endif