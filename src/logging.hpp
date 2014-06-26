#ifndef _LOGGING_H
#define _LOGGING_H


#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

#include <boost/format.hpp>

using namespace std;

static string algorithm_name_for_logging = "LDA";

template<typename T>
T *allocate_1d_with_default(size_t dim, T def) {
    T *r = new T[dim];
    for (size_t i=0; i < dim; ++i) {
        r[i] = def;
    }
    return r;
}

template<typename T>
T **allocate_2d_with_default(size_t dim1, size_t dim2, T def) {
    T **r = new T*[dim1];
    for (size_t i=0; i < dim1; ++i) {
        r[i] = allocate_1d_with_default(dim2, def);
    }
    return r;
}

template<typename T>
void delete_2d(size_t dim_one, T **r) {
    for (size_t i=0; i<dim_one; ++i) {
        delete r[i];
    }
    delete r;
}

struct perf_counter
{
    perf_counter() { reset(); }

    double time() const {
        struct timeval now;
        gettimeofday(&now, NULL);
        return (now.tv_sec - time_.tv_sec) + (now.tv_usec - time_.tv_usec) * 1e-6;
    }

    double time_ms() const {
        return time() * 1e3;
    }

    void reset() {
        gettimeofday(&time_, NULL);
    }

private:
    struct timeval time_;
};

static perf_counter timer;

inline unsigned get_max_rss() {
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  return ru.ru_maxrss;
}

inline std::string human_readable_memory(unsigned max_rss)
{
        assert(max_rss > 0);

        if (max_rss < 1024 * 1024) {
                return str(boost::format("%4dM") % (max_rss / 1024));
    } else {
    return str(boost::format("%4dG") % (max_rss / (1024 * 1024)));
        }
}

inline std::string human_readable_time(double time_in_sec) {
    assert(time_in_sec > 0);
    size_t msec  = size_t(time_in_sec * 1000) % 1000;
    size_t sec   = size_t(time_in_sec);
    size_t hours = sec / 3600;
    size_t mins  = (sec / 60) % 60;
    sec         %= 60;
    return str(boost::format("%3d:%02d:%02d.%03d") % hours % mins % sec % msec);
}


#define LOG(msg) cerr << "[" << str(boost::format("%14s %5s") % human_readable_time(timer.time()) % human_readable_memory(get_max_rss())) << "]\t" << msg << endl;

#endif
