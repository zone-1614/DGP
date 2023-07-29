#include <iostream>
#include <vector>
#include <set>
#include <thread>
#include <omp.h>

#include <Eigen/Sparse>

#include <pmp/algorithms/DifferentialGeometry.h>

#include "../util.h"

namespace zone {

class DeformationViewer : public MyViewer {
public:
    DeformationViewer() : MyViewer("Deformation Viewer"), state(0)
    {
        set_draw_mode("Hidden Line");
    }

    void process_imgui() {
        MyViewer::process_imgui();
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("hw6", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::Button("camera")) {
                state = 0;
            }
            if (ImGui::Button("select fixed")) {
                state = 1;
            }
            if (ImGui::Button("select handle")) {
                state = 2;
            }
            if (ImGui::Button("clear fixed")) {
                std::cout << "clear fixed" << std::endl;
                fixed.clear();
            }
            if (ImGui::Button("clear handle")) {
                std::cout << "clear handle" << std::endl;
                handle.clear();
            }
            ImGui::InputFloat3("delta x y z", delta_p.data());
            if (ImGui::Button("deformation")) {
                deformation();
                update_mesh();
                view_all();
            }
        }
    }

    void motion(double xpos, double ypos) override {
        if (state == 0) {
            TrackballViewer::motion(xpos, ypos);
        }
    }

    void mouse(int button, int action, int mods) override {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            double x, y;
            cursor_pos(x, y);
            auto v = pick_vertex(static_cast<int>(x), static_cast<int>(y));
            if (!v.is_valid()) 
                return ;
            // 不能重复添加点, 所以用set, 而不是vector
            // 如果该点已经在fixed或者handle中, 就不再添加
            bool exists_in_fixed = fixed.count(v);
            bool exists_in_handle = handle.count(v);

            if (exists_in_fixed || exists_in_handle)
                return;

            if (state == 1) {
                fixed.insert(v);
                auto p = mesh_.position(v);
                std::cout << "add fixed point, idx: " << v.idx() << ", coord: " << p << std::endl;
            } else if (state == 2) {
                handle.insert(v);
                auto p = mesh_.position(v);
                std::cout << "add handle point, idx: " << v.idx() << ", coord: " << p << std::endl;
            }
        }
    }

    void deformation() {
        // 线程数
        // 把handle移动到目的位置
        for (auto v : handle) {
            mesh_.position(v) += delta_p;
        }

        auto cot = compute_cot_weight();
        // 论文里面的eij, 因为是常量, 所以在一开始全部算出来, 避免重复计算
        auto eij = compute_eij();
        const auto origin_pos = mesh_.positions(); // 原来的坐标, 后面需要用到

        // global 阶段需要用的变量, 提前声明了, 防止重复创建
        const int nv = mesh_.n_vertices();
        Eigen::MatrixX3f uv = Eigen::MatrixX3f::Zero(nv, 3);
        Eigen::MatrixX3f b  = Eigen::MatrixX3f::Zero(nv, 3);

        // global阶段的拉普拉斯矩阵也是常量, 提前创建了
        std::vector<Eigen::Triplet<float>> L_tri;
        L_tri.reserve(nv * 5);
        for (auto v : mesh_.vertices()) {
            // 对于handle中的点, 拉普拉斯矩阵对应的行只有一个1
            if (handle.count(v)) {
                L_tri.emplace_back(v.idx(), v.idx(), 1.0f);
                continue;
            }
            float sum = 0.0f;
            for (auto h : mesh_.halfedges(v)) {
                // 跳过边界
                if (mesh_.is_boundary(h)) {
                    continue;
                }
                float w = cot[h] + cot[mesh_.opposite_halfedge(h)];
                sum += w;
                L_tri.emplace_back(v.idx(), mesh_.to_vertex(h).idx(), -w);
            }
            L_tri.emplace_back(v.idx(), v.idx(), sum);
        }
        Eigen::SparseMatrix<float> LaplacianMatrix(nv, nv);
        LaplacianMatrix.setFromTriplets(L_tri.begin(), L_tri.end());
        LaplacianMatrix.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<float>> solver(LaplacianMatrix);

        auto Lt = mesh_.vertex_property<Eigen::Matrix3f>("Lt");
        for (int i = 0; i < 10; i++) { // 迭代次数
            // local 阶段, 求Lt
            for (auto v : mesh_.vertices()) {
                auto from_p = mesh_.position(v);
                // 对每一个顶点, 求出jacobian
                Eigen::Matrix3f J = Eigen::Matrix3f::Zero();
                for (auto h : mesh_.halfedges(v)) {
                    auto to_v = mesh_.to_vertex(h);
                    auto to_p = mesh_.position(to_v);
                    auto eij_ = from_p - to_p;
                    auto weight = cot[h] + cot[mesh_.opposite_halfedge(h)];
                    Eigen::Vector3f E_(eij_[0], eij_[1], eij_[2]);
                    Eigen::Vector3f E(eij[h][0], eij[h][1], eij[h][2]);
                    J += weight * (E * E_.transpose());
                }
                // 对Jacobian做SVD, 得到R = VU'
                Eigen::JacobiSVD<Eigen::Matrix3f> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Matrix3f U = svd.matrixU();
                Eigen::Matrix3f V = svd.matrixV();
                Eigen::Matrix3f R = V * U.transpose();
                // 保证 det(R) > 0, 这样才是旋转, 要不然可能是翻转
                if (R.determinant() < 0) {
                    U.col(2) *= -1.0f;
                    R = V * U.transpose();
                }
                Lt[v] = R;
            }

            // global 阶段, 利用Lt去优化能量, 得到新的坐标, 下一个local阶段又根据新的坐标得到新的Lt
            for (auto v : mesh_.vertices()) {
                Eigen::Vector3f b_ = Eigen::Vector3f::Zero();
                for (auto h : mesh_.halfedges(v)) {
                    auto to_v = mesh_.to_vertex(h);
                    Eigen::Vector3f E(eij[h][0], eij[h][1], eij[h][2]);
                    Eigen::Matrix3f JR = Lt[v] + Lt[to_v];
                    auto weight = (cot[h] + cot[mesh_.opposite_halfedge(h)]) / 2.0f;
                    b_ += weight * (JR * E);
                }
                b.row(v.idx()) += b_;
            }

            // handle和fixed的b改回原来的
            for (auto v : handle) {
                Eigen::Vector3f origin_b;
                origin_b[0] = origin_pos[v.idx()][0];
                origin_b[1] = origin_pos[v.idx()][1];
                origin_b[2] = origin_pos[v.idx()][2];
                b.row(v.idx()) = origin_b;
            }
            for (auto v : fixed) {
                Eigen::Vector3f origin_b;
                origin_b[0] = origin_pos[v.idx()][0];
                origin_b[1] = origin_pos[v.idx()][1];
                origin_b[2] = origin_pos[v.idx()][2];
                b.row(v.idx()) = origin_b;
            }
            uv.col(0) = solver.solve(b.col(0));
            uv.col(1) = solver.solve(b.col(1));
            uv.col(2) = solver.solve(b.col(2));

            for (auto v : mesh_.vertices()) {
                mesh_.position(v) = pmp::Point(uv(v.idx(), 0), uv(v.idx(), 1), uv(v.idx(), 2));
            }
        }
    }

    auto compute_cot_weight() -> pmp::HalfedgeProperty<float> {
        auto cot = mesh_.halfedge_property<float>("cot_weight", 1.0f);
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
            cot[h0] = pmp::cotan(p0 - p2, p1 - p2);
            cot[h1] = pmp::cotan(p1 - p0, p2 - p0);
            cot[h2] = pmp::cotan(p0 - p1, p2 - p1);
        }
        return cot;
    }

    auto compute_eij() -> pmp::HalfedgeProperty<pmp::Point> {
        auto eij = mesh_.halfedge_property<pmp::Point>("eij");
        for (auto h : mesh_.halfedges()) {
            auto to = mesh_.to_vertex(h), from = mesh_.from_vertex(h);
            auto to_p = mesh_.position(to), from_p = mesh_.position(from);
            eij[h] = from_p - to_p;
        }
        return eij;
    }

    int state; // 0 表示正常模式(移动相机), 1表示选固定点, 2表示选handle, 3表示
    std::set<pmp::Vertex> fixed;
    std::set<pmp::Vertex> handle;
    pmp::Point delta_p{0,0,0};
};

}

int main() {
    zone::DeformationViewer viewer;
    viewer.run();
}