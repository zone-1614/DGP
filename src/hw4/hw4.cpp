#include <iostream>
#include <vector>

#include <Eigen/Sparse>

#include <pmp/algorithms/DifferentialGeometry.h>

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
            ImGui::SliderInt("ARAP max iter times", &ARAP_max_iter_times, 1, 50);
            if (ImGui::Button("ARAP")) {
                ARAP();
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

    void ARAP() {
        // 1. execute tutte embedding before ARAP
        tutte_embedding(); 
        std::cout << "tutte" << std::endl;
        // each faces area (multiply by 2.0)
        auto area = compute_all_area();
        std::cout << "compute area" << std::endl;
        // each halfedge's cot_weight
        auto cot = compute_all_cot();
        std::cout << "compute cot" << std::endl;
        // init uv as xy
        auto ARAP_uv = init_uv();
        std::cout << "init uv" << std::endl;

        // 2. Global matrix G is constant
        Eigen::SparseMatrix<float> G = build_global_matrix();
        std::cout << "build_global_matrix" << std::endl;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver(G);

        // // Matrix X also constant
        // Eigen::MatrixX2f X = build_X();

        // 3. iteration
        float error = 10.0f;
        auto jacobian = mesh_.face_property<Eigen::Matrix2f>("jacobian");
        auto L = mesh_.face_property<Eigen::Matrix2f>("ARAP_L");
        int k = 0;
        while (k++ < ARAP_max_iter_times) {
            // Local phase
            // for each triangle, namely each face, using u,x to compute Jacobian
            for (auto f : mesh_.faces()) {
                jacobian[f] = compute_jacobian(f);
                // SVD
                Eigen::JacobiSVD<Eigen::Matrix2f> svd(jacobian[f], Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Matrix2f U = svd.matrixU(), V = svd.matrixV();
                // Eigen::Vector2f Sigma = svd.singularValues();
                Eigen::Matrix2f Lf = U * V.transpose();

                // keep det(Lf) = 1
                if (Lf.determinant() < 0.0f) {
                    if (Lf(0, 0) > 0.0f) {
                        Lf.col(1) *= -1.0f;
                    } else {
                        Lf.col(0) *= -1.0f;
                    }
                }
                L[f] = Lf;
            }

            // Global phase
            // build matrix B
            Eigen::MatrixX2f B = build_B();

            // solve GU = B,  U = [u, v]
            Eigen::MatrixX2f U = solver.solve(B);
            for (auto v : mesh_.vertices()) {
                int i = v.idx();
                pmp::Point uv{ U(i, 0), U(i, 1), 0.0f };
                ARAP_uv[v] = uv;
            }

            // update the error
            error = 0.0f;
            for (auto h : mesh_.halfedges()) {
                float cot_ij = cot[h];
                auto vi = mesh_.from_vertex(h), vj = mesh_.to_vertex(h);
                auto uvi = ARAP_uv[vi], uvj = ARAP_uv[vj];
                auto xyi = mesh_.position(vi), xyj = mesh_.position(vj);
                Eigen::Vector2f UV;
                UV << uvi[0] - uvj[0], uvi[1] - uvj[1];
                Eigen::Vector2f X;
                X << xyi[0] - xyj[0], xyi[1] - xyj[1];
                auto f = mesh_.face(h);
                if (!f.is_valid()) {
                    continue;
                }
                Eigen::Vector2f E = UV - L[f] * X;
                error += cot_ij * E.squaredNorm();
            }
            std::cout << "error: " << error << ", iter: " << k << " times" << std::endl;
        }

        // update vertex position
        for (auto v : mesh_.vertices()) {
            mesh_.position(v) = ARAP_uv[v];
        }
    }

    pmp::FaceProperty<float> compute_all_area() {
        auto area = mesh_.face_property<float>("ARAP_area");
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
            area[f] = pmp::norm(pmp::cross(e0, e1));
        }
        return area;
    }
    
    pmp::HalfedgeProperty<float> compute_all_cot() {
        auto cot = mesh_.halfedge_property<float>("ARAP_cot");
        for (auto f : mesh_.faces()) {
            auto h0 = mesh_.halfedge(f);
            auto h1 = mesh_.next_halfedge(h0);
            auto h2 = mesh_.next_halfedge(h1);
            auto v0 = mesh_.to_vertex(h0);
            auto v1 = mesh_.to_vertex(h1);
            auto v2 = mesh_.to_vertex(h2);
            auto p0 = mesh_.position(v0);
            auto p1 = mesh_.position(v1);
            auto p2 = mesh_.position(v2);

            cot[h0] = pmp::cotan(p0 - p1, p2 - p1);
            cot[h1] = pmp::cotan(p1 - p2, p0 - p2);
            cot[h2] = pmp::cotan(p1 - p0, p2 - p0);
        }
        return cot;
    }

    pmp::VertexProperty<pmp::Point> init_uv() {
        auto ARAP_uv = mesh_.vertex_property<pmp::Point>("ARAP_uv");
        for (auto v : mesh_.vertices()) {
            ARAP_uv[v] = mesh_.position(v);
        }
        return ARAP_uv;
    }

    Eigen::SparseMatrix<float> build_global_matrix() {
        int nv = mesh_.n_vertices();
        std::vector<Eigen::Triplet<float>> G_tri;
        G_tri.reserve(nv * 5);
        Eigen::SparseMatrix<float> G(nv, nv);

        auto cot = mesh_.halfedge_property<float>("ARAP_cot");
        for (auto v : mesh_.vertices()) {
            float d = 0.0f; // diagnose element
            int i = v.idx();
            for (auto h : mesh_.halfedges(v)) {
                int j = mesh_.to_vertex(h).idx();
                float w = 0.0f;
                w += cot[h];
                w += cot[mesh_.opposite_halfedge(h)];
                G_tri.push_back({ i, j, -w });
                d += w;
            }
            G_tri.push_back({ i, i, d });
        }

        G.setFromTriplets(G_tri.begin(), G_tri.end());
        G.makeCompressed();
        return G;
    }

    // Eigen::MatrixX2f build_X() {
    //     int nv = mesh_.n_vertices();
    //     Eigen::MatrixX2f X;
    //     for (auto v : mesh_.vertices()) {
    //         int i = v.idx();
    //         auto p = mesh_.position(v);
    //         X(i, 0) = p[0];
    //         X(i, 1) = p[1];
    //     }
    //     return X;
    // }
    Eigen::MatrixX2f build_B() {
        int nv = mesh_.n_vertices();
        Eigen::MatrixX2f B = Eigen::MatrixX2f::Zero(nv, 2);
        auto L_ = mesh_.face_property<Eigen::Matrix2f>("ARAP_L");
        auto cot = mesh_.halfedge_property<float>("ARAP_cot");
        for (auto v : mesh_.vertices()) {
            int i = v.idx();
            auto xi = mesh_.position(v);
            Eigen::Vector2f t;
            t.setZero();
            for (auto h : mesh_.halfedges(v)) {
                Eigen::Matrix2f w;
                w.setZero();
                auto f = mesh_.face(h);
                auto oh = mesh_.opposite_halfedge(h);
                auto of = mesh_.face(oh);
                auto vj = mesh_.to_vertex(h);
                auto xj = mesh_.position(vj);
                if (f.is_valid() && h.is_valid()) {
                    w += L_[f] * cot[h];
                }
                if (of.is_valid() && oh.is_valid()) {
                    w += cot[oh] * L_[of];
                }
                Eigen::Vector2f x;
                x << xi[0] - xj[0], xi[1] - xj[1];
                t += w * x;
            }
            B(i, 0) = t(0);
            B(i, 1) = t(1);
        }
        return B;
    }

    Eigen::Matrix2f compute_jacobian(const pmp::Face& f) {
        auto area = mesh_.face_property<float>("ARAP_area");
        auto ARAP_uv = mesh_.vertex_property<pmp::Point>("ARAP_uv");
        auto hi = mesh_.halfedge(f);
        auto hj = mesh_.next_halfedge(hi);
        auto hk = mesh_.next_halfedge(hj);
        auto vi = mesh_.to_vertex(hk);
        auto vj = mesh_.to_vertex(hi);
        auto vk = mesh_.to_vertex(hj);
        auto xyi = mesh_.position(vi);
        auto xyj = mesh_.position(vj);
        auto xyk = mesh_.position(vk);

        auto uvi = ARAP_uv[vi], uvj = ARAP_uv[vj], uvk = ARAP_uv[vk];
        Eigen::Matrix2f j;
        j << (xyj[1] - xyk[1]) * uvi[0] + (xyk[1] - xyi[1]) * uvj[0] + (xyi[1] - xyj[1]) * uvk[0],
             (xyj[0] - xyk[0]) * uvi[0] + (xyk[0] - xyi[0]) * uvj[0] + (xyi[0] - xyj[0]) * uvk[0],
             (xyj[1] - xyk[1]) * uvi[0] + (xyk[1] - xyi[1]) * uvj[0] + (xyi[1] - xyj[1]) * uvk[0],
             (xyj[0] - xyk[0]) * uvi[0] + (xyk[0] - xyi[0]) * uvj[0] + (xyi[0] - xyj[0]) * uvk[0];
        return j;
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

private:
    int ARAP_max_iter_times = 10;
};

}

int main() {
    zone::ParameterazationViewer window;
    window.run();
}
