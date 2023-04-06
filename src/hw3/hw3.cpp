#include <iostream>

#include <pmp/MatVec.h>
#include <pmp/visualization/MeshViewer.h>
#include <pmp/io/io.h>
#include <pmp/algorithms/Normals.h>
#include <imgui.h>

#include "../util.h"
using namespace zone;

namespace zone {

class DenoisingViewer : public zone::MyViewer {
public:
    DenoisingViewer() 
        : MyViewer("denoising viewer") { }

protected:
    void process_imgui() {
        MyViewer::process_imgui();
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("hw3", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::SliderInt("normal iter times", &normal_iter_times, 0, 1000, "%d");
            ImGui::SliderInt("vertex iter times", &vertex_iter_times, 0, 1000, "%d");
            ImGui::SliderFloat("sigma s", &sigma_s, 0.0, 1.0, "%.3f", 1.0f);
            ImGui::SliderFloat("sigma r", &sigma_r, 0.0, 1.0, "%.3f", 1.0f);
            if (ImGui::Button("Add noise")) {
                add_noise();
            }
            if (ImGui::Button("Denoising")) {
                denoising();
            }
            update_mesh();
        }
    }

private:
    void add_noise() {
        // 计算平均的边长
        const int n = mesh_.n_edges();
        double avg_edge_length = 0.0;
        for (auto e : mesh_.edges()) {
            const auto v0 = mesh_.vertex(e, 0);
            const auto v1 = mesh_.vertex(e, 1);
            const auto p0 = mesh_.position(v0);
            const auto p1 = mesh_.position(v1);
            avg_edge_length += pmp::norm(p0 - p1);
        }
        avg_edge_length /= (double)n;

        // 用平均边长的 20% 作为扰动向量的长度
        double l = 0.2 * avg_edge_length;
        for (auto v : mesh_.vertices()) {
            auto p = mesh_.position(v);            
            set_point(v, p + random_unit_vector() * l);
        }
    }

    // 只考虑封闭的mesh
    void denoising() {
        // mesh 的体积, 用于滤波后的 scale
        double old_volume = compute_volume();
        // 初始化面的法向
        auto new_fnormals = mesh_.face_property<pmp::Normal>("new_fnormal");
        auto old_fnormals = mesh_.face_property<pmp::Normal>("old_fnormal");
        for (auto f : mesh_.faces()) {
            old_fnormals[f] = pmp::Normals::compute_face_normal(mesh_, f);
        }
        std::cout << "init face normal" << std::endl;

        int n = mesh_.n_faces();
        double error = 0.0;
        for (int i = 0; i < normal_iter_times; i++) {
            error = 0.0;
            // 更新面的法向
            for (auto fi : mesh_.faces()) {
                // ni: normal of fi
                auto ni = old_fnormals[fi];
                // ci: center of fi
                auto h0 = mesh_.halfedge(fi);
                auto h1 = mesh_.next_halfedge(h0);
                auto h2 = mesh_.next_halfedge(h1);
                auto pi = mesh_.position(mesh_.to_vertex(h0));
                auto qi = mesh_.position(mesh_.to_vertex(h1));
                auto ri = mesh_.position(mesh_.to_vertex(h2));
                auto ci = (pi + qi + ri) / 3.0;

                double total_weight = 0.0;
                pmp::Point new_ni{0.0, 0.0, 0.0}; // n_{i}^{t+1}
                for (auto h : mesh_.halfedges(fi)) {
                    auto oh = mesh_.opposite_halfedge(h);
                    auto fj = mesh_.face(oh);
                    // compute Aj: area of fj
                    auto p = mesh_.position(mesh_.to_vertex(oh));
                    auto q = mesh_.position(mesh_.from_vertex(oh));
                    auto r = mesh_.position(mesh_.to_vertex(mesh_.next_halfedge(oh)));
                    auto pq = q - p, pr = r - p;
                    auto Aj = 0.5 * pmp::norm(pmp::cross(pq, pr));

                    // Ws 
                    auto cj = (p + q + r) / 3.0;
                    auto Ws = gaussian(pmp::norm(ci - cj), sigma_s);
                    // Wr
                    auto nj = old_fnormals[fj]; // nj: normal of fj
                    auto Wr = gaussian(pmp::norm(ni - nj), sigma_r);

                    double weight = Aj * Ws * Wr;
                    total_weight += weight;
                    new_ni += weight * nj;
                }
                new_fnormals[fi] = new_ni;
            }

            // 将法向归一化
            for (auto f : mesh_.faces()) {
                new_fnormals[f] = pmp::normalize(new_fnormals[f]);
            }

            // 计算误差
            for (auto f : mesh_.faces()) {
                error += pmp::dot(new_fnormals[f], old_fnormals[f]);
                old_fnormals[f] = new_fnormals[f];
            }
            error /= (double)n;
            std::cout << "normal error: " << error << std::endl;

        }

        for (int i = 0; i < vertex_iter_times; i++) {
            // 更新顶点
            for (auto vi : mesh_.vertices()) {
                // number of face around vi
                int Ni = 0;
                pmp::Point delta{0.0, 0.0, 0.0};
                auto xi = mesh_.position(vi);
                for (auto fj : mesh_.faces(vi)) {
                    Ni++;
                    auto nj = new_fnormals[fj]; // normal of fj
                    // cj: center of fj
                    auto h0 = mesh_.halfedge(fj);
                    auto h1 = mesh_.next_halfedge(h0);
                    auto h2 = mesh_.next_halfedge(h1);
                    auto p0 = mesh_.position(mesh_.to_vertex(h0));
                    auto p1 = mesh_.position(mesh_.to_vertex(h1));
                    auto p2 = mesh_.position(mesh_.to_vertex(h2));
                    auto cj = (p0 + p1 + p2) / 3.0;
                    delta += (pmp::dot(nj, cj - xi) * nj);
                }
                delta /= (double)Ni;
                set_point(vi, xi + delta);
            }

            // 计算误差
            error = 0.0;
            for (auto fk : mesh_.faces()) {
                auto nk = new_fnormals[fk];
                for (auto h : mesh_.halfedges(fk)) {
                    auto xi = mesh_.position(mesh_.to_vertex(h)), xj = mesh_.position(mesh_.from_vertex(h));
                    error += pmp::dot(xi - xj, nk) * pmp::dot(xi - xj, nk);
                }
            }
            std::cout << "vertex error " << error << std::endl;
        }

        // 5. 缩放, 保持体积不变
        double new_volume = compute_volume();
        double k = old_volume / new_volume;
        for (auto v : mesh_.vertices()) {
            set_point(v, mesh_.position(v) * k);
        }
    }

    double compute_volume() {
        double sum = 0.0;
        for (auto f : mesh_.faces()) {
            auto h0 = mesh_.halfedge(f);
            auto h1 = mesh_.next_halfedge(h0);
            auto h2 = mesh_.next_halfedge(h1);
            auto p0 = mesh_.position(mesh_.to_vertex(h0));
            auto p1 = mesh_.position(mesh_.to_vertex(h1));
            auto p2 = mesh_.position(mesh_.to_vertex(h2));
            Eigen::Matrix3d m;
            m << p0[0], p0[1], p0[2],
                 p1[0], p1[1], p1[2],
                 p2[0], p2[1], p2[2];
            sum += m.determinant();
        }
        return std::abs(sum) / 6.0;
    }

private:
    float sigma_s = 0.1, sigma_r = 0.1;
    int normal_iter_times = 10, vertex_iter_times = 10;
};


}
int main() {
    zone::DenoisingViewer window;
    window.run();
}