#include <iostream>

#include <pmp/MatVec.h>
#include <pmp/visualization/MeshViewer.h>
#include <pmp/io/read_obj.h>
#include <imgui.h>

#include "../util.h"

class CurvatureViewer : public pmp::MeshViewer {
public:
    CurvatureViewer() 
        : pmp::MeshViewer("curvature viewer", 900, 600) {
        std::filesystem::path model_path(zone::current_path.string() + "/model");
        for (auto& file : std::filesystem::directory_iterator(model_path)) {
            models.push_back(file.path().stem().string());
        }
        std::string first = zone::current_path.string() + "/model/" + models.front() + ".obj";
        load_mesh(first.c_str());
    }

protected:
    void process_imgui() {
        pmp::MeshViewer::process_imgui();

        ImGui::Spacing();
        if (ImGui::CollapsingHeader("model", ImGuiTreeNodeFlags_DefaultOpen)) {
            static int item_current_idx = 0;
            const char* combo_preview_value = (models[item_current_idx]).c_str();
            if (ImGui::BeginCombo("models", combo_preview_value))
            {
                for (int n = 0; n < models.size(); n++)
                {
                    const bool is_selected = (item_current_idx == n);
                    if (ImGui::Selectable((models[n]).c_str(), is_selected)) {
                        std::string filename = zone::current_path.string() + "/model/" + models[n] + ".obj";
                        this->load_mesh(filename.c_str());
                        item_current_idx = n;
                    }

                    if (is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            if (ImGui::Button("reload mesh")) {
                load_mesh(filename_.c_str());
                set_draw_mode("Smooth Shading");
            }
        }

        ImGui::Spacing();
        if (ImGui::CollapsingHeader("curvature type", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::RadioButton("mean curvature", &curvature_type, 0);
            ImGui::RadioButton("Gaussian curvature", &curvature_type, 1);
            if (ImGui::Button("execute")) {
                if (curvature_type == 0) {
                    visualize_mean_curvature();
                } else {
                    visualize_Gaussian_curvature();
                }
            }
        }
    }

private:
    // 平均权拉普拉斯  效果不好, 不用了
    pmp::Point avg_laplacian(pmp::Vertex v) {
        auto pos = mesh_.position(v);
        pmp::Point lap{0.0, 0.0, 0.0};
        for (auto vv : mesh_.vertices(v)) {
            lap += (mesh_.position(vv) - pos);
        }
        return lap;
    }

    // cot 权拉普拉斯
    pmp::Point laplacian(pmp::Vertex v) {
        auto pos = mesh_.position(v);
        pmp::Point lap{0.0, 0.0, 0.0};
        for (auto h : mesh_.halfedges(v)) {
            double weight = 0.0;
            auto oh = mesh_.opposite_halfedge(h);
            auto p = mesh_.position(mesh_.to_vertex(h));
            auto op = mesh_.position(mesh_.to_vertex(oh));
            
            if (!mesh_.is_boundary(h)) {
                auto pp = mesh_.position(mesh_.to_vertex(mesh_.next_halfedge(h)));
                auto l1 = p - pp, l2 = op - pp;
                double cot_h = (pmp::dot(l1, l2)) / pmp::norm(pmp::cross(l1, l2));
                cot_h = std::clamp(cot_h, 0.0, 11.43);
                weight += cot_h;
            }

            if (!mesh_.is_boundary(oh)) {
                auto pp = mesh_.position(mesh_.to_vertex(mesh_.next_halfedge(oh)));
                auto l1 = p - pp, l2 = op - pp;
                double cot_oh = (pmp::dot(l1, l2)) / pmp::norm(pmp::cross(l1, l2));
                cot_oh = std::clamp(cot_oh, 0.0, 11.43);
                weight += cot_oh;
            }

            lap += weight * (p - pos);
        }
        lap /= 2.0 * local_averaging_region(v);
        return lap;
    }

    double local_averaging_region(pmp::Vertex v) {
        double region = 0.0;
        for (auto h : mesh_.halfedges(v)) {
            // 默认是 ccw
            if (mesh_.is_boundary(h)) continue;
            /*
                     p
                  /     \
                h0       h2
                /          \
               q --- h1 --- r
            */
            auto h0 = h;
            auto h1 = mesh_.next_halfedge(h0);
            auto h2 = mesh_.next_halfedge(h1);

            // 点
            // auto p = mesh_.from_vertex(h0);
            // auto q = mesh_.from_vertex(h1);
            // auto r = mesh_.from_vertex(h2);
            // to_vertex 性能更快
            auto p = mesh_.position(mesh_.to_vertex(h2));
            auto q = mesh_.position(mesh_.to_vertex(h0));
            auto r = mesh_.position(mesh_.to_vertex(h1));

            // 边
            auto pq = q - p, pr = r - p, qr = r - p;

            // 判断点p对应的角(以下简称角p)是不是钝角, 是的话用对边中点替代外心
            // 如果角q或者角r为钝角, 取对边中点代替外心之后, 变成一个三角形
            // pq 与 pr 点乘之后为:  |pq| * |pr| * cos p, 记为cosp
            double cosp = pmp::dot(pq, pr), cosq = -pmp::dot(pq, pr), cosr = pmp::dot(pr, qr); 
            double area = pmp::norm(pmp::cross(pq, pr)); // 三角形面积
            if (cosp < 0.0) {
                region += 0.25 * area;
            } else if (cosq < 0.0 || cosr < 0.0) {
                region += 0.125 * area;
            } else {
                // 三个角都不是钝角, 求 cot q, cot r
                double cotq = cosq / area, cotr = cosr / area;
                // 把 q r clamp 到 [5, 175] 度之间
                cotq = std::clamp(cotq, -11.43, 11.43);
                cotr = std::clamp(cotr, -11.43, 11.43);
                area += 0.125 * cotq * pmp::sqrnorm(pq);
                area += 0.125 * cotr * pmp::sqrnorm(pr);
            }
        }
        return region;
    }

    double mean_curvature(pmp::Vertex v) {
        return 0.5 * pmp::norm(laplacian(v));
    }

    double Gaussian_curvature(pmp::Vertex v) {
        if (mesh_.is_boundary(v)) 
            return 0.0;

        double total_angle = 0.0;
        auto p0 = mesh_.position(v);
        for (auto h : mesh_.halfedges(v)) {
            auto p1 = mesh_.position(mesh_.to_vertex(h));
            auto next_h = mesh_.ccw_rotated_halfedge(h);
            auto p2 = mesh_.position(mesh_.to_vertex(next_h));
            // normalize
            auto l1 = p1 - p0, l2 = p2 - p0;
            double cos_theta = pmp::dot(l1, l2) / ( pmp::norm(l1) * pmp::norm(l2) );
            total_angle += acos(cos_theta);
        }
        return (2 * EIGEN_PI - total_angle) / local_averaging_region(v);
    }

    void visualize_mean_curvature() {
        // 计算每一个顶点的平均曲率
        // 不要直接用最大最小的曲率, 应该去掉最大的k%和最小的k%(这个k取10或者5), 要不然偏差太大了
        auto mean_curv = mesh_.vertex_property<double>("v:mean_curv");
        std::vector<double> mean_curv_;
        const int n = mesh_.n_vertices();
        mean_curv_.reserve(n);
        for (auto v : mesh_.vertices()) {
            auto curv = mean_curvature(v);
            mean_curv[v] = curv;
            mean_curv_.push_back(curv);
        }

        std::sort(mean_curv_.begin(), mean_curv_.end());
        double k = 10.0 / 100.0;
        int kk = k * n;
        const double min_curv = mean_curv_[kk], max_curv = mean_curv_[n - 1 - kk];
        const double d_curv = max_curv - min_curv;
        // 将曲率映射为纹理
        auto texture = mesh_.vertex_property<pmp::TexCoord>("v:tex");
        for (auto v : mesh_.vertices()) {
            double x = (mean_curv[v] - min_curv) / d_curv;
            x = std::clamp(x, 0.0, 1.0); // 防止出现小于0.0和大于1.0的数
            texture[v] = pmp::TexCoord(x, 0.0);
        }

        mesh_.use_cold_warm_texture();
        update_mesh();
        set_draw_mode("Texture");
    }

    void visualize_Gaussian_curvature() {
        // 计算每一个顶点的高斯曲率
        auto Gaussian_curv = mesh_.vertex_property<double>("v:Gaussian_curv");
        std::vector<double> Gaussian_curv_;
        const int n = mesh_.n_vertices();
        Gaussian_curv_.reserve(n);
        for (auto v : mesh_.vertices()) {
            auto curv = Gaussian_curvature(v);
            Gaussian_curv[v] = curv;
            Gaussian_curv_.push_back(curv);
        }

        std::sort(Gaussian_curv_.begin(), Gaussian_curv_.end());
        double k = 5.0 / 100.0;
        int kk = k * n;
        const double min_curv = Gaussian_curv_[kk], max_curv = Gaussian_curv_[n - 1 - kk];
        const double d_curv = max_curv - min_curv;
        // 将曲率映射为纹理
        auto texture = mesh_.vertex_property<pmp::TexCoord>("v:tex");
        for (auto v : mesh_.vertices()) {
            double x = (Gaussian_curv[v] - min_curv) / d_curv;
            x = std::clamp(x, 0.0, 1.0); // 防止出现小于0.0和大于1.0的数
            texture[v] = pmp::TexCoord(x, 0.0);
        }

        mesh_.use_cold_warm_texture();
        update_mesh();
        set_draw_mode("Texture");
    }

private:
    std::vector<std::string> models;
    int curvature_type = 0;
};

int main() {
    CurvatureViewer window;
    window.run();
}