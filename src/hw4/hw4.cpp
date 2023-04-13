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
        // get boundary vertices
        auto boundary_v = get_boundary_vertices();

        // map the boundary_v to a unit circle
        const int n = boundary_v.size();
        const float step = 2.0f * EIGEN_PI / (float) n;
        float theta = 0.0f;
        for (int i = 0; i < n; i++) {
            pmp::Point p{ std::cosf(theta), std::sinf(theta), 0.0f };
            theta += step;
            auto v = boundary_v[i];
            set_point(v, p);
        }

        // build and solve the linear system
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

        // update the coordinate
        for (auto v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) {
                continue;
            } else {
                pmp::Point p{ X(v.idx(), 0), X(v.idx(), 1), 0.0 };
                set_point(v, p);
            }
        }
    }

    std::vector<pmp::Vertex> get_boundary_vertices() {
        std::vector<pmp::Vertex> boundary_v;
        pmp::Halfedge first, current;

        // find first boundary halfedge
        for (auto h : mesh_.halfedges()) {
            if (mesh_.is_boundary(h)) {
                boundary_v.push_back(mesh_.to_vertex(h));
                first = h;
                current = h;
                break;
            }
        }

        // find all boundary vertex
        do {
            for (auto hv : mesh_.halfedges(mesh_.to_vertex(current))) {
                if (!mesh_.is_boundary(hv)) {
                    continue;
                } else if (mesh_.opposite_halfedge(hv) == current) {
                    continue;
                } else {
                    current = hv;
                    boundary_v.push_back(mesh_.to_vertex(hv));
                }
            }
        } while (current != first);

        return boundary_v;
    }
};

}

int main() {
    zone::ParameterazationViewer window;
    window.run();
}