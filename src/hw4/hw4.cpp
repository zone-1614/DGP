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
            if (ImGui::Button("tutte")) {
                tutte_embedding();
                update_mesh();
                view_all();
            }
            if (ImGui::Button("LSCM")) {
                LSCM();
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
            mesh_.position(v) = p;
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
                mesh_.position(v) = p;
            }
        }
    }

    void LSCM() {
        // 1. get to boundary vertex
        auto boundary = get_boundary_vertices();
        // two pinned boundary vertex 
        auto pv1 = boundary.front(), pv2 = boundary[boundary.size() / 2];

        // the pinned vertex index
        int pinned_idx1 = pv1.idx(), pinned_idx2 = pv2.idx();
        if (pinned_idx1 > pinned_idx2) {
            std::swap(pinned_idx1, pinned_idx2);
        }

        // pinned coord
        Eigen::Vector4f p;
        p << 0.0f, 1.0f, 0.0f, 1.0f;

        // 2. compute the auxiliary value. face area, project vector of halfedge
        compute_auxiliary();
        auto area = mesh_.face_property<float>("my_area");
        auto project_vec = mesh_.halfedge_property<pmp::Point>("project_vector");

        // Sparse Matrices
        int nv = mesh_.n_vertices(), nf = mesh_.n_faces();
        std::vector<Eigen::Triplet<Eigen::scomplex>> tri_M;
        std::vector<Eigen::Triplet<float>> tri_A, tri_B;
        tri_M.reserve(nf * 3);
        tri_A.reserve(nf * 3);
        Eigen::SparseMatrix<float> A(2 * nf, 2 * nv - 4), B(2 * nf, 4);

        // fill the tri_M
        for (auto f : mesh_.faces()) {
            auto sqrt_area = std::sqrt(area[f]);
            int fidx = f.idx();
            for (auto h : mesh_.halfedges(f)) {
                auto v = mesh_.to_vertex(mesh_.next_halfedge(h));
                int vidx = v.idx();
                auto pv = project_vec[h];
                float real = pv[0] / sqrt_area, imag = pv[1] / sqrt_area;
                // fill the tri_M
                tri_M.push_back({ fidx, vidx, { real, imag } });
            }
        }
        
        // use tri_M to build tri_A and tri_B
        for (auto tri : tri_M) {
            float real = tri.value().real(), imag = tri.value().imag();
            int row = tri.row(), col = tri.col();
            if (col == pinned_idx1) {
                tri_B.push_back({ row, 0, real });
                tri_B.push_back({ row + nf, 2, real });
                tri_B.push_back({ row, 2, -imag });
                tri_B.push_back({ row + nf, 0, imag });
            } else if (col == pinned_idx2) {
                tri_B.push_back({ row, 1, real });
                tri_B.push_back({ row + nf, 3, real });
                tri_B.push_back({ row, 3, -imag });
                tri_B.push_back({ row + nf, 1, imag });
            } else {
                if (col > pinned_idx2) col--;
                if (col > pinned_idx1) col--;
                tri_A.push_back({ row, col, real });
                tri_A.push_back({ row + nf, col + nv - 2, real });
                tri_A.push_back({ row, col + nv - 2, -imag });
                tri_A.push_back({ row + nf, col, imag });
            }
        }
        A.setFromTriplets(tri_A.begin(), tri_A.end());
        B.setFromTriplets(tri_B.begin(), tri_B.end());
        tri_M.clear();
        tri_A.clear();
        tri_B.clear();

        Eigen::VectorXf b = -B * p;
        
        std::cout << "solving the matrix" << std::endl;
        // solve
        A.makeCompressed();
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solver(A);
        Eigen::VectorXf x = solver.solve(b);

        // update coordinate
        for (auto v : mesh_.vertices()) {
            int idx = v.idx();
            if (idx == pinned_idx1) {
                pmp::Point p{ 0.0f, 0.0f ,0.0f };
                mesh_.position(v) = p;
            } else if (idx == pinned_idx2) {
                pmp::Point p{ 1.0f, 1.0f ,0.0f };
                mesh_.position(v) = p;
            } else {
                if (idx > pinned_idx2) idx--;
                if (idx > pinned_idx1) idx--;
                pmp::Point p{ x(idx), x(idx + nv - 2), 0.0f };
                mesh_.position(v) = p;
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

    // compute all faces area and the vertex--halfedge project vector 
    void compute_auxiliary() {
        auto area = mesh_.face_property<float>("my_area");
        auto project_vector = mesh_.halfedge_property<pmp::Point>("project_vector");
        for (auto f : mesh_.faces()) {
            auto h0 = mesh_.halfedge(f);
            auto h1 = mesh_.next_halfedge(h0);
            auto h2 = mesh_.next_halfedge(h1);
            auto v0 = mesh_.to_vertex(h2);
            auto v1 = mesh_.to_vertex(h0);
            auto v2 = mesh_.to_vertex(h1);
            auto p0 = mesh_.position(v0);
            auto p1 = mesh_.position(v1);
            auto p2 = mesh_.position(v2);
            auto e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
            auto ne0 = pmp::norm(e0), ne2 = pmp::norm(e2);

            area[f] = pmp::norm(pmp::cross(e0, e1));
            pmp::Point x0{ 0.0f, 0.0f, 0.0f }, x1{ 0.0f, 0.0f, 0.0f }, x2{ 0.0f, 0.0f, 0.0f };
            x1[0] = ne0;
            x2[0] = pmp::dot(e0, -e2) / ne0;
            x2[1] = pmp::norm(pmp::cross(e0, -e2)) / ne0;
            
            project_vector[h0] = x1 - x0;
            project_vector[h1] = x2 - x1;
            project_vector[h2] = x0 - x2;
        }
    }

};

}

int main() {
    zone::ParameterazationViewer window;
    window.run();
}
