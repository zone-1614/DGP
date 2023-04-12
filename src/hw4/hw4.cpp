#include <iostream>
#include <vector>

#include <Eigen/Sparse>

#include "../util.h"

namespace zone {

class ParameterazationViewer : public MyViewer {
public:
    ParameterazationViewer() : MyViewer("Parameterazation Viewer") { 
        set_draw_mode("Hidden Line");
    }

protected:
    void process_imgui() {
        MyViewer::process_imgui();

        ImGui::Spacing();
        if (ImGui::CollapsingHeader("hw4", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::Button("execute")) {
                tutte_embedding();
                update_mesh();
                view_all();
            }
        }
    }

private:
    void tutte_embedding() {
        // 找到boundary
        std::vector<pmp::Halfedge> boundary;
        pmp::Halfedge first, current;

        // 找第一条边界
        for (auto h : mesh_.halfedges()) {
            if (mesh_.is_boundary(h)) {
                boundary.push_back(h);
                first = h;
                break;
            }
        }

        current = first;
        // 沿着这条边界找所有的边界
        do {
            auto v = mesh_.to_vertex(current);
            for (auto hv : mesh_.halfedges(v)) {
                if (!mesh_.is_boundary(hv))
                    continue;
                else if (mesh_.opposite_halfedge(hv) == current)
                    continue;
                
                boundary.push_back(hv);
                current = hv;
            }
        } while (current != first);

        std::vector<pmp::Vertex> boundary_v;
        boundary_v.resize(boundary.size());
        for (int i = 0; i < boundary.size(); i++) {
            boundary_v[i] = mesh_.to_vertex(boundary[i]);
        }
        // 把boundary上的点映射到圆上
        const int n = boundary_v.size();
        const float step = 2.0f * EIGEN_PI / (float) n;
        float theta = 0.0f;
        for (int i = 0; i < n; i++) {
            pmp::Point p{ std::cosf(theta), std::sinf(theta), 0.0f };
            theta += step;
            auto v = boundary_v[i];
            set_point(v, p);
        }

        // 解方程组
        std::vector<Eigen::Triplet<float>> tri;
        Eigen::MatrixXf B(mesh_.n_vertices(), 2);
        B.setZero();
        for (auto v : mesh_.vertices()) {
            int i = v.idx() ;

            if (mesh_.is_boundary(v)) {
                tri.push_back({ i, i, 1.0f });
                auto p = mesh_.position(v);
                B(i, 0) = p[0];
                B(i, 1) = p[1];
            } else {
                float ii = 0.0f;
                for (auto vv : mesh_.vertices(v)) {
                    int j = vv.idx();
                    tri.push_back({ i, j, -1.0f });
                    ii += 1.0f;
                }
                tri.push_back({ i, i, ii });
            }
        }
        Eigen::SparseMatrix<float> A(mesh_.n_vertices(), mesh_.n_vertices());
        A.setFromTriplets(tri.begin(), tri.end());
        
        Eigen::SparseLU<Eigen::SparseMatrix<float>> solver(A);
        Eigen::MatrixXf X = solver.solve(B);

        for (auto v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) {
                continue;
            } else {
                pmp::Point p{ X(v.idx(), 0), X(v.idx(), 1), 0.0 };
                set_point(v, p);
            }
        }
    }
};

}

int main() {
    zone::ParameterazationViewer window;
    window.run();
}