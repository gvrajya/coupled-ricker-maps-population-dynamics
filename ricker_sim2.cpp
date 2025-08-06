#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <map>
#include <atomic>   // For thread-safe counters and flags
#include <chrono>   // For timestamps
#include <thread>   // For the live progress thread
#include <csignal>  // For catching Ctrl+C
#include <sstream>  // Added for std::stringstream
#include <set>      // For tracking completed jobs
#include <tuple>    // For storing parameter sets

// Use OpenMP for parallelization
#include <omp.h>

// --- Global flags for graceful exit ---
std::atomic<bool> computation_finished{false};
std::atomic<bool> exit_flag_raised{false};

// --- Signal handler for Ctrl+C ---
void handle_signal(int signum) {
    if (signum == SIGINT) {
        std::cout << "\nCtrl+C detected. Shutting down gracefully. Please wait..." << std::endl;
        exit_flag_raised = true;
    }
}

// A simple struct to hold the results of one parameter set
struct PointResult {
    double epsilon, local_sigma, global_scale, shock_prob;
    double mean_tau, std_tau;
};

// Global random number generator setup
std::mt19937& get_rng() {
    thread_local std::mt19937 rng(std::random_device{}() + omp_get_thread_num());
    return rng;
}

// Core simulation function (remains the same)
double get_tau_for_single_run(double r, double epsilon, double local_sigma, double global_scale, double shock_prob) {
    const int num_steps = 40000;
    const int transient = 100;
    const double STABILITY_LIMIT = 1e9;
    const int MAX_LAG = 2000;

    std::vector<double> X1(num_steps, 0.0), X2(num_steps, 0.0);
    X1[0] = 0.2; X2[0] = 1.3;

    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    for (int t = 0; t < num_steps - 1; ++t) {
        if (exit_flag_raised) return NAN; // Check for exit signal
        double g_scale_factor = (uniform_dist(get_rng()) < shock_prob) ? (1.0 * global_scale) : (0.1 * global_scale);
        double g = normal_dist(get_rng()) * g_scale_factor;
        auto R_ricker = [&](double x) {
            double clamped_x = std::max(0.0, std::min(10.0, x));
            return clamped_x * std::exp(r * (1.0 - clamped_x));
        };
        auto eta = [&](double x) { return x * std::exp(normal_dist(get_rng()) * local_sigma); };
        double r1 = R_ricker(X1[t]), r2 = R_ricker(X2[t]);
        X1[t + 1] = eta((1.0 - epsilon) * r1 + epsilon * r2) + g;
        X2[t + 1] = eta((1.0 - epsilon) * r2 + epsilon * r1) + g;
        if (!std::isfinite(X1[t + 1]) || !std::isfinite(X2[t + 1]) || std::abs(X1[t + 1]) > STABILITY_LIMIT) return NAN;
    }

    // --- Analysis Part (unchanged) ---
    std::vector<double> X_bar;
    for (int i = transient; i < num_steps; ++i) X_bar.push_back(X1[i] + X2[i]);
    if (X_bar.size() < 2) return NAN;
    std::vector<double> M;
    for (size_t i = 0; i < X_bar.size() - 1; ++i) M.push_back(std::pow(-1, i + 1) * (X_bar[i + 1] - X_bar[i]));
    double m_mean = std::accumulate(M.begin(), M.end(), 0.0) / M.size();
    double m_sq_sum = std::inner_product(M.begin(), M.end(), M.begin(), 0.0);
    double m_std = std::sqrt(m_sq_sum / M.size() - m_mean * m_mean);
    if (m_std < 1e-6) return NAN;
    int max_lag = std::min(MAX_LAG, static_cast<int>(M.size()) - 1);
    double m_norm_sq_sum = 0.0;
    for(double val : M) m_norm_sq_sum += (val - m_mean) * (val - m_mean);
    if (m_norm_sq_sum == 0) return NAN;
    std::vector<double> acf_vals(max_lag);
    for (int k = 0; k < max_lag; ++k) {
        double cross_product = 0.0;
        for (size_t t = 0; t < M.size() - k; ++t) cross_product += (M[t] - m_mean) * (M[t + k] - m_mean);
        acf_vals[k] = cross_product / m_norm_sq_sum;
    }
    int end_idx = max_lag;
    for (int i = 1; i < max_lag; ++i) {
        if (acf_vals[i] <= 0.0) {
            end_idx = i;
            break;
        }
    }
    std::vector<double> x_fit, y_fit;
    for (int lag = 1; lag < end_idx; ++lag) {
        double acf_val = acf_vals[lag];
        if (acf_val > 0.0) {
            x_fit.push_back(static_cast<double>(lag));
            y_fit.push_back(acf_val);
        }
    }
    if (x_fit.size() < 2) return NAN;
    double min_sse = std::numeric_limits<double>::infinity();
    double best_tau = NAN;
    const int NUM_GRID_POINTS = 10000;
    const double MIN_TAU = 0.1;
    const double MAX_TAU = 10000.0;
    double log_min = std::log(MIN_TAU);
    double log_max = std::log(MAX_TAU);
    for (int i = 0; i < NUM_GRID_POINTS; ++i) {
        double log_tau_val = log_min + (log_max - log_min) * static_cast<double>(i) / (NUM_GRID_POINTS - 1);
        double tau_candidate = std::exp(log_tau_val);
        double sse = 0.0;
        for (size_t j = 0; j < x_fit.size(); ++j) {
            double arg = -x_fit[j] / tau_candidate;
            double pred = (arg > -700.0) ? std::exp(arg) : 0.0;
            double resid = y_fit[j] - pred;
            sse += resid * resid;
        }
        if (sse < min_sse) {
            min_sse = sse;
            best_tau = tau_candidate;
        }
    }
    if (!std::isfinite(best_tau) || best_tau <= 0.0) return NAN;
    return best_tau;
}

// Helper to generate a linear space of numbers
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> vec;
    if (num <= 1) { vec.push_back(start); return vec; }
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) vec.push_back(start + delta * i);
    return vec;
}

int main() {
    signal(SIGINT, handle_signal);

    const int N_RUNS_PER_POINT = 25;
    const int EPSILON_POINTS = 20;
    const int LOCAL_SIGMA_POINTS = 15;
    const int GLOBAL_SCALE_POINTS = 8;
    const int SHOCK_PROB_POINTS = 5;
    
    // --- KEY SECTION 1: Reading the existing CSV file ---
    const std::string OUTPUT_FILENAME = "RickerSim_2025-07-18_23-24-39.csv";
    std::set<std::tuple<double, double, double, double>> completed;
    std::ifstream infile(OUTPUT_FILENAME);
    if (infile) {
        std::string line;
        std::getline(infile, line); // Skip header row

        while (std::getline(infile, line)) {
            std::stringstream ss(line);
            std::string token;
            double params[4];
            try {
                // Read the four parameter values from the line
                for (int i = 0; i < 4; ++i) {
                    if (!std::getline(ss, token, ',')) throw std::runtime_error("Incomplete line");
                    params[i] = std::stod(token);
                }
                // Add the successfully read parameter set to our 'completed' set
                completed.insert(std::make_tuple(params[0], params[1], params[2], params[3]));
            } catch (const std::exception& e) {
                // This catches errors if a line is malformed, making the reader more robust.
                // You can uncomment the line below to see warnings for skipped lines.
                // std::cerr << "Warning: Skipping malformed line in CSV: " << line << std::endl;
            }
        }
        infile.close();
    }

    // --- Parameter Ranges ---
    auto epsilon_range = linspace(0.05, 0.15, EPSILON_POINTS);
    auto local_sigma_range = linspace(0.05, 0.15, LOCAL_SIGMA_POINTS);
    auto global_scale_range = linspace(0.0, 0.05, GLOBAL_SCALE_POINTS);
    auto shock_prob_range = linspace(0.0, 0.02, SHOCK_PROB_POINTS);

    // Create a vector of all parameter combinations
    struct ParamSet { double eps, ls, gs, sp; };
    std::vector<ParamSet> all_params;
    for (double eps : epsilon_range) for (double ls : local_sigma_range) for (double gs : global_scale_range) for (double sp : shock_prob_range) {
        all_params.push_back({eps, ls, gs, sp});
    }
    
    const long long total_params = all_params.size();
    std::atomic<long long> params_done{completed.size()}; // Start counter from already completed size

    std::cout << "--- 🚀 Resuming 4D Parameter Sweep ---" << std::endl;
    std::cout << "Total parameter sets: " << total_params << std::endl;
    std::cout << "Already completed: " << completed.size() << std::endl;
    std::cout << "Appending to: " << OUTPUT_FILENAME << std::endl;

    // Open file in append mode
    std::ofstream outfile(OUTPUT_FILENAME, std::ios::app);
    
    std::thread progress_thread([&]() {
        while (!computation_finished) {
            std::cout << "\rProgress: " << params_done << " / " << total_params 
                      << " parameter sets (" << std::fixed << std::setprecision(1) 
                      << (100.0 * params_done / total_params) << "%)" << std::flush;
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
    });
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < all_params.size(); ++i) {
        if (exit_flag_raised) continue;

        auto p = all_params[i];

        // --- KEY SECTION 2: Skipping the work ---
        // Before running the simulation, check if the parameter set exists in our 'completed' set.
        if (completed.count(std::make_tuple(p.eps, p.ls, p.gs, p.sp))) {
            // If it exists, we skip this iteration of the loop.
            continue;
        }

        // The growth rate 'r' is held constant for this 4D sweep
        const double r_const = 2.3;

        std::vector<double> results_for_point;
        for (int k = 0; k < N_RUNS_PER_POINT; ++k) {
            results_for_point.push_back(get_tau_for_single_run(r_const, p.eps, p.ls, p.gs, p.sp));
        }
        
        std::vector<double> valid_results;
        for(double res : results_for_point) if(std::isfinite(res)) valid_results.push_back(res);

        double mean_tau = NAN, std_tau = NAN;
        if (valid_results.size() >= static_cast<size_t>(N_RUNS_PER_POINT / 2)) {
            mean_tau = std::accumulate(valid_results.begin(), valid_results.end(), 0.0) / valid_results.size();
            double sq_sum = std::inner_product(valid_results.begin(), valid_results.end(), valid_results.begin(), 0.0);
            std_tau = std::sqrt(sq_sum / valid_results.size() - mean_tau * mean_tau);
        }
        
        #pragma omp critical
        {
            if(std::isfinite(mean_tau)) {
                 outfile << p.eps << "," << p.ls << "," << p.gs << ","
                         << p.sp << "," << mean_tau << "," << std_tau << "\n";
            }
        }
        params_done++;
    }

    computation_finished = true;
    progress_thread.join();

    // Final status message
    std::cout << "\n---------------------------------------------------" << std::endl;
    if (exit_flag_raised) {
        std::cout << "--- 🛑 Computation stopped early. Partial data saved to " << OUTPUT_FILENAME << " ---" << std::endl;
    } else {
        std::cout << "--- ✅ Simulations complete. Full data saved to " << OUTPUT_FILENAME << " ---" << std::endl;
    }
    
    return 0;
}