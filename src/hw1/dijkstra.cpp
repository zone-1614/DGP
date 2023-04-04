#include <vector>
#include <algorithm>

#include <pmp/SurfaceMesh.h>
#include <pmp/io/read_obj.h>
#include <spdlog/spdlog.h>

#include "../util.h"

class dijkstra {
public: 
    dijkstra(std::string filename = "Bunny_head.obj") {
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

    // 输入一个点的index, 返回它到所有点的最短路径
    std::vector<double> run(int idx) {
        if (idx >= vertices_.size()) {
            return { };
        }
        // 储存 distance 和 visited, 初始化距离为 DBL_MAX, 并且未访问
        std::vector<std::pair<double, bool>> dis_vis(vertices_.size(), { DBL_MAX, false });
        dis_vis[idx] = { 0.0, true };
        auto index = mesh_.get_vertex_property<int>("v:my_index");
        // 取得源点的handle, 初始化距离
        auto v = vertices_[idx];
        for (auto vv : mesh_.vertices(v)) {
            int vv_idx = index[vv];
            dis_vis[vv_idx].first = distance(v, vv);
        }

        // 循环 size - 1 次, 得到所有最短路径
        for (int a = 0; a < dis_vis.size() - 1; a++) {
            // 找到当前距离最近且 vis 为 false 的顶点的下标
            int min_index = 0;
            double min_dis = DBL_MAX;
            for (int i = 0; i < dis_vis.size(); i++) {
                if (dis_vis[i].second) continue;
                if (dis_vis[i].first < min_dis) {
                    min_index = i;
                    min_dis = dis_vis[i].first;
                }
            }

            // 把找到的点的vis状态标为 true
            dis_vis[min_index].second = true;

            // 更新其他点的距离
            for (int i = 0; i < dis_vis.size(); i++) {
                if (dis_vis[i].second) continue;
                double new_distance = min_dis + distance(vertices_[i], vertices_[min_index]);
                if (new_distance < dis_vis[i].first) {
                    dis_vis[i].first = new_distance;
                }
            }
        }

        std::vector<double> ans(dis_vis.size());
        std::transform(dis_vis.begin(), dis_vis.end(), ans.begin(), [](const std::pair<double, bool>& pair) {
            return pair.first;
        });
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
    dijkstra d(filename);
    auto dis = d.run(4);

    for (int i = 0; i < dis.size(); i++) {
        spdlog::info("index {} to index 4, distance {}", i, dis[i]);
    }
}
