#include <vector>
#include <algorithm>
#include <tuple>

#include <pmp/SurfaceMesh.h>
#include <pmp/io/read_obj.h>
#include <spdlog/spdlog.h>

#include "../util.h"

class prim {
public:
    prim(std::string filename = "Balls.obj") {
        // 读 mesh
        auto fullname = zone::current_path.string() + "/model/" + filename;
        auto path = zone::fs::path(fullname);
        pmp::read_obj(mesh_, path);        

        // 把 mesh 的顶点存到vector中, 并添加属性index
        auto index = mesh_.add_vertex_property<int>("v:my_index");
        int i = 0;
        for (auto v : mesh_.vertices()) {
            vertices_.push_back(v);
            index[v] = i++;
        }
    }

    // 默认以index 0 开始
    std::vector<std::vector<int>> run() {
        // init
        auto index = mesh_.get_vertex_property<int>("v:my_index");
        std::vector<std::vector<int>> ans(vertices_.size());
        std::vector<std::tuple<int, double, bool>> prev_dis_vis(vertices_.size(), { -1, DBL_MAX, false });
        prev_dis_vis[0] = { 0, 0, true };
        for (auto vv : mesh_.vertices(vertices_[0])) {
            prev_dis_vis[index[vv]] = { 0, distance(vv, vertices_[0]), false };
        }

        // loop
        for (int a = 1; a < vertices_.size(); a++) {
            double min_dis = DBL_MAX;
            double min_index = 1;
            for (int i = 1; i < vertices_.size(); i++) {
                if (std::get<2>(prev_dis_vis[i])) continue;
                if (min_dis > std::get<1>(prev_dis_vis[i])) {
                    min_dis = std::get<1>(prev_dis_vis[i]);
                    min_index = i;
                }
            }

            std::get<2>(prev_dis_vis[min_index]) = true;
            int prev = std::get<0>(prev_dis_vis[min_index]);
            ans[prev].push_back(min_index);

            for (auto vv : mesh_.vertices(vertices_[min_index])) {
                int idx = index[vv];
                if (std::get<2>(prev_dis_vis[idx])) continue;
                double new_distance = distance(vv, vertices_[min_index]);
                if (new_distance < std::get<1>(prev_dis_vis[idx])) {
                    prev_dis_vis[idx] = { min_index, new_distance, false };
                }
            }
        }
        return ans;
    }
private: 
    double distance(const pmp::Vertex& v1, const pmp::Vertex& v2) {
        return pmp::distance(mesh_.position(v1), mesh_.position(v2));
    }

private:
    pmp::SurfaceMesh mesh_;
    std::vector<pmp::Vertex> vertices_;
};


int main(int argc, char** argv) {
    std::string filename;
    if (argc == 1) {
        filename = "Balls.obj";
    } else {
        filename = argv[1];
    }
    prim d(filename);
    auto ans = d.run();

    for (auto vec : ans) {
        for (auto v : vec) {
            std::cout << v << "  ";
        }
        std::cout << std::endl;
    }
}

