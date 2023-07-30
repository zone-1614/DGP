#pragma once

#include <filesystem>
#include <iostream>
#include <random>
#include <functional>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <vector>
#include <queue>

#include <pmp/MatVec.h>
#include <pmp/visualization/MeshViewer.h>
#include <pmp/io/io.h>
#include <imgui.h>

namespace zone {

namespace fs = std::filesystem;
static fs::path current_path = fs::current_path();

std::ostream& operator<<(std::ostream& os, pmp::Point p) {
    return os << p[0] << ", " << p[1] << ", " << p[2];
}

void elapse(std::function<void()> f) {
    auto t0 = std::chrono::steady_clock::now();
    f();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = t1 - t0;
    std::cout << "cost " << std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() << "ms" << std::endl;
}

// return a random double in [0.0, 1.0)
double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

pmp::Point random_vector(double min, double max) {
    pmp::Point p;
    p[0] = random_double(min, max);
    p[1] = random_double(min, max);
    p[2] = random_double(min, max);
    return p;
}

pmp::Point random_unit_vector() {
    while (true) {
        auto p = random_vector(-1.0, 1.0);
        if (pmp::norm(p) >= 1.0) continue;
        return pmp::normalize(p);
    }
}

double gaussian(double x, double sigma) {
    return std::exp(- (x * x) / (sigma * sigma));
}

class MyViewer : public pmp::MeshViewer {
public:
    MyViewer(std::string title) 
        : pmp::MeshViewer(title.c_str(), 900, 600) {
        std::filesystem::path model_path(zone::current_path.string() + "/model");
        for (auto& file : std::filesystem::directory_iterator(model_path)) {
            models_.push_back(file.path().stem().string());
        }
        std::string first = zone::current_path.string() + "/model/" + models_.front() + ".obj";
        load_mesh(first.c_str());
    }

protected:
    void process_imgui() {
        pmp::MeshViewer::process_imgui();
        
        ImGui::SetNextItemWidth(150);
        ImGui::PushItemWidth(150);

        ImGui::Spacing();
        if (ImGui::CollapsingHeader("model", ImGuiTreeNodeFlags_DefaultOpen)) {
            static int item_current_idx = 0;
            const char* combo_preview_value = (models_[item_current_idx]).c_str();
            if (ImGui::BeginCombo("models", combo_preview_value))
            {
                for (int n = 0; n < models_.size(); n++)
                {
                    const bool is_selected = (item_current_idx == n);
                    if (ImGui::Selectable((models_[n]).c_str(), is_selected)) {
                        std::string filename = zone::current_path.string() + "/model/" + models_[n] + ".obj";
                        load_mesh(filename.c_str());
                        item_current_idx = n;
                    }

                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            if (ImGui::Button("reload mesh")) {
                load_mesh(filename_.c_str());
                set_draw_mode("Hidden Line");
            }
        }

        ImGui::Spacing();
        if (ImGui::CollapsingHeader("Save as", ImGuiTreeNodeFlags_DefaultOpen)) {
            static char filename[128] = "save.obj";
            ImGui::InputText("filename", filename, IM_ARRAYSIZE(filename));
            if (ImGui::Button("save")) {
                // save mesh to file
                std::filesystem::path p(zone::current_path.string() + "/model/" + filename);
                pmp::write(mesh_, p);
            }
        }
    }

private:
    std::vector<std::string> models_;
};

class ThreadPool {
public:
    ThreadPool(int n_threads = 16)
        : stop(false)
        // , n_running(0) 
    {
        for (int i = 0; i < n_threads; i++) {
            threads.push_back(std::thread([this] {
                std::function<void()> task;
                while (true) {
                    {// acquire lock
                        std::unique_lock lck(this->mtx);
                        while (!this->stop && this->tasks.empty()) {
                            this->cond.wait(lck);
                        }
                        if (this->stop)
                            return;
                        task = this->tasks.front();
                        this->tasks.pop();
                    }// release lock
                    // this->n_running++;
                    task();
                    // this->n_running--;
                    wait_cond.notify_one();
                }
            }));
        }
    }
    void enqueue(std::function<void()> f) {
        {// acquire lock
            std::unique_lock lck(mtx);
            tasks.push(f);
        }// release lock
        cond.notify_one();
    }
    void waitAll() {
        std::unique_lock lck(mtx);
        wait_cond.wait(lck, [this] {
            // return this->tasks.empty() && this->n_running == 0;
            return this->tasks.empty();
        });
    }
    ~ThreadPool() {
        stop = true;
        cond.notify_all();
        for (auto& t : threads)
            t.join();
    }
private:
    std::mutex mtx;
    std::condition_variable cond;
    std::condition_variable wait_cond; // 用于等待所有线程执行完毕
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> tasks;
    // std::atomic_uint n_running;
    bool stop;
};

}
