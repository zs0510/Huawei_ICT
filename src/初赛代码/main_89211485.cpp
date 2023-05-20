#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include <ctime>
#include <unordered_set>
#include <unordered_map>
#include <limits.h>

using namespace std;

const int edges_add_max = 20000;// 可以增加的最大边数

struct Edge {
    int id, distance, nid0, nid1;
    vector<bool> channel_used;// 记录通道是否被使用
    Edge() {}
    Edge(int d, int _n0, int _n1): distance(d), nid0(_n0), nid1(_n1) {}
};

struct Node {
    int id;
    vector<int> eids;
    Node() {}
};

struct Operation {
    int id, nid0, nid1, channel_id = -1;
    vector<int> edges, amplifiers;// 业务经过的边ID, 以及经过放大器的节点ID
    Operation() {}
};

int num_of_nodes, num_of_edges_initial, num_of_operations, num_of_channels, distance_attenuation_max;
vector<Node> nodes;// num_of_nodes 是节点个数, ID 从[0, num_of_nodes-1]
vector<Edge> edges;// num_of_edges 是边个数, ID 从[0, num_of_edges-1]
vector<Operation> operations;// num_of_operations 是业务个数, ID 从[0, num_of_operations-1]
unordered_map<int, unordered_map<int, int>> dists_min;// 记录两个节点之间直接连接的最短距离: 新生成的边的长度
unordered_map<int, unordered_map<int, vector<int>>> paths_of_nodes;// 以节点序列记录两个节点之间的路径

vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id = -1);// 返回源点与目标点之间的节点路径(channel_id 为 -1 时不考虑通道占用情况)
bool is_channel_id_ok(vector<vector<int>>& path_edges_total, int channel_id, vector<int>& path_edges);// 检查这些边的通道 pid 是否被占用
int get_eid_channel_id_ok(int nid0, int nid1, int channel_id);// 检查节点(直接相连)之间能否通过这个通道进行连接
int get_channel_id_empty_max(vector<vector<int>>& path_edges);
int add_edge(int nid0, int nid1);// 新增两个节点之间的边
vector<vector<int>> get_path_edges_total_via_nodes(vector<int>& path_nodes);// 返回节点序列中所有可能的边
void sort_operations();


int main() {
    // 输入
    cin >> num_of_nodes >> num_of_edges_initial >> num_of_operations >> num_of_channels >> distance_attenuation_max;
    nodes.resize(num_of_nodes);
    edges.resize(num_of_edges_initial);
    operations.resize(num_of_operations);
    for (int i = 0; i < num_of_nodes; ++i) {
        nodes[i].id = i;
    }
    for (int i = 0; i < num_of_edges_initial; ++i) {// 读入 num_of_edges 行边
        cin >> edges[i].nid0 >> edges[i].nid1 >> edges[i].distance;
        edges[i].id = i;
        edges[i].channel_used = vector<bool>(num_of_channels, false);
        nodes[edges[i].nid0].eids.push_back(edges[i].id);
        nodes[edges[i].nid1].eids.push_back(edges[i].id);
        // 记录最短距离
        if (dists_min[edges[i].nid0].count(edges[i].nid1)) {
            dists_min[edges[i].nid0][edges[i].nid1]
            = dists_min[edges[i].nid1][edges[i].nid0]
                    = min(edges[i].distance, dists_min[edges[i].nid0][edges[i].nid1]);
        } else {
            dists_min[edges[i].nid0][edges[i].nid1]
            = dists_min[edges[i].nid1][edges[i].nid0]
                    = edges[i].distance;
        }

    }

    for (int i = 0; i < num_of_operations; ++i) {// 读入 num_of_operations 行业务
        cin >> operations[i].nid0 >> operations[i].nid1;
        operations[i].id = i;
    }

#ifdef __APPLE__
    // 输出建图信息
    cout << "\n";
    for (int i = 0; i < num_of_nodes; ++i) {
        auto& node = nodes[i];
        cout << "Node ID: " << node.id << ", ";
        for (int j = 0; j < node.eids.size(); ++j) {
            auto& edge = edges[node.eids[j]];
            cout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
        }
        cout << "\n";
    }
#endif

    // 对每个节点的邻接边进行排序, 把距离短的放在前面
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }

    // 对业务进行排序
    // 1. 把路径可选择的少的排前面
    // 2. 把长度长的排在前面
    sort_operations();


    // 处理
    for (auto& op : operations) {
//        auto& op = operations[op_id];
        int nid_source = op.nid0, nid_target = op.nid1;

        vector<int> path_nodes;// 这条路径通过的节点

        if (paths_of_nodes[nid_source].count(nid_target)) {
            path_nodes = paths_of_nodes[nid_source][nid_target];
        } else {
            // 初始化
            // 使用广搜得到可以到达的路径(不区分通道是否被使用完)
            path_nodes = bfs_find_path_nodes(nid_source, nid_target);
            paths_of_nodes[nid_source][nid_target] = path_nodes;
            auto tmp = path_nodes;
            reverse(tmp.begin(), tmp.end());
            paths_of_nodes[nid_target][nid_source] = tmp;
        }
        vector<vector<int>> path_edges_total = get_path_edges_total_via_nodes(path_nodes);// 节点之间的所有的边

        // 检查路径是否可达
        int path_channel_id = -1;// 本次业务使用的通道 ID
        vector<int> path_edges;// 这条路径通过的边

        // 检查可用的通道编号, 找到第一个可用的
        for (int channel_id = 0; channel_id < num_of_channels; ++channel_id) {
            if (is_channel_id_ok(path_edges_total, channel_id, path_edges)) {
                path_channel_id = channel_id;
                break;
            }
        }

        if (path_channel_id == -1) {
            // 当前节点序列不满足要求
#ifdef __APPLE__
            cout << "path_channel_id is -1\n";
//            for (int i = 0; i < path_edges_total.size(); ++i) {
//                cout << "(";
//                for (int j = 0; j < path_edges_total[i].size(); ++j) {
//                    cout << path_edges_total[i][j];
//                    if (j + 1 != path_edges_total[i].size()) cout << ",";
//                }
//                cout << ")";
//                if (i + 1 == path_edges_total.size()) cout << "\n";
//                else cout << "->";
//            }
#endif
            // solution 1: 重新组织一下节点序列, 以检查看是否可以在不新增边的前提下完成需求
            for (int channel_id = 0; channel_id < num_of_channels; ++channel_id) {
                auto path_nodes_tmp = bfs_find_path_nodes(nid_source, nid_target, channel_id);
                if (!path_nodes_tmp.empty()) {
                    auto path_edges_total_tmp = get_path_edges_total_via_nodes(path_nodes_tmp);
                    if (is_channel_id_ok(path_edges_total_tmp, channel_id, path_edges)) {
                        // 确认没问题, 再统一修改
                        path_nodes = path_nodes_tmp;
                        path_edges_total = path_edges_total_tmp;
                        path_channel_id = channel_id;
                        break;
                    }
                }
            }

            // solution 2: 无法通过现有边完成业务时, 新增边以满足要求
            if (path_channel_id == -1) {
                // 简单方法, 把所有边都新增一条
                path_edges.clear();
                path_channel_id = get_channel_id_empty_max(path_edges_total);// 统计每条边的通道使用情况, 选出需要增加边最少的一种情况来增加边

                for (int i = 1; i < path_nodes.size(); ++i) {
                    // 检查是否存在通道为0且空闲的边, 尽量少建新边
                    int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id);
                    if (eid_old != -1) {
                        path_edges.emplace_back(eid_old);
                        continue;
                    }
                    int eid_new = add_edge(path_nodes[i - 1], path_nodes[i]);
                    path_edges.emplace_back(eid_new);
                }
            }

        }
        // 请保证 path_edges 及 pid_edges 是可行的
#ifdef __APPLE__
        cout << "nid_source = " << nid_source << ", nid_target = " << nid_target << ", path = ";
        for (int i = 0; i < path_nodes.size(); ++i) {
            cout << path_nodes[i];
            if (i + 1 == path_nodes.size()) cout << "\n";
            else cout << "->";
        }
#endif
        // 将这些边的通道标为已占用
        for (auto& eid : path_edges) {
            auto& edge = edges[eid];
#ifdef __APPLE__
            if (edge.channel_used[path_channel_id]) {
                cerr << "Error: Reusing the same channel\n";
                continue;
            }
#endif
            edge.channel_used[path_channel_id] = true;
        }
        // 放置放大器
        vector<int> path_amplifiers;
        int dist_sum = edges[path_edges.front()].distance;
        for (int i = 1; i < path_edges.size(); ++i) {
            // dist_sum 表示当前到达这个节点的线路长度
            auto& edge = edges[path_edges[i]];
            if (dist_sum + edge.distance >= distance_attenuation_max) {
                path_amplifiers.push_back(path_nodes[i]);
                dist_sum = edge.distance;
            } else {
                dist_sum += edge.distance;
            }
        }
        op.channel_id = path_channel_id;
        op.edges = path_edges;
        op.amplifiers = path_amplifiers;
    }


    // 输出
    // 输出加边信息
    cout << edges.size() - num_of_edges_initial << endl;
    for (int i = num_of_edges_initial; i < edges.size(); ++i) {
        cout << edges[i].nid0 << " " << edges[i].nid1 << endl;
    }
    // 输出业务通道信息
    sort(operations.begin(), operations.end(), [&](Operation& op0, Operation& op1) {
        return op0.id < op1.id;
    });
    for (int i = 0; i < num_of_operations; ++i) {
        auto& operation = operations[i];
        // 通道ID, 边数目, 放大器数目
        cout << operation.channel_id << " " << operation.edges.size() << " " << operation.amplifiers.size();
        for (auto& eid : operation.edges) {
            cout << " " << eid;
        }
        for (auto& nid : operation.amplifiers) {
            cout << " " << nid;
        }
        cout << endl;
    }

    return 0;
}

vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id) {
    // TODO: 使用双向搜索来加速 https://oi-wiki.org/search/bidirectional/
    unordered_map<int, int> prev;
    queue<int> que;
    que.push(nid_source);
    prev[nid_source] = -1;
    while (!que.empty() && !prev.count(nid_target)) {
        int nid = que.front();
        que.pop();
        for (auto& eid : nodes[nid].eids) {
            auto& edge = edges[eid];
            if (channel_id != -1 && edge.channel_used[channel_id]) continue;// 考虑通道使用情况
            int nid_next = edge.nid0 + edge.nid1 - nid;
            if (prev.count(nid_next)) continue;
            que.push(nid_next);
            prev[nid_next] = nid;
        }
    }

    if (!prev.count(nid_target)) {// 无法在当前通道编号下连接 nid_source 与 nid_target
        return {};
    }

    vector<int> path_nodes;
    int nid = nid_target;
    while (nid != -1) {
        path_nodes.push_back(nid);
        nid = prev[nid];
    }

    reverse(path_nodes.begin(), path_nodes.end());

#ifdef __APPLE__
    cout << "path: ";
    for (int i = 0; i < path_nodes.size(); ++i) {
        cout << path_nodes[i];
        if (i + 1 == path_nodes.size()) cout << "\n";
        else cout << "->";
    }
#endif
    return path_nodes;
}

bool is_channel_id_ok(vector<vector<int>>& path_edges_total, int channel_id, vector<int>& path_edges) {

    path_edges.clear();
    int i = 0;
    while (i < path_edges_total.size()) {
        bool is_empty = false;
        for (auto& eid : path_edges_total[i]) {
            auto& edge = edges[eid];
            if (edge.channel_used[channel_id]) continue;
            is_empty = true;
            path_edges.emplace_back(eid);
            break;
        }
        if (!is_empty) break;
        ++i;
    }

    if (i != path_edges_total.size()) {
        path_edges.clear();
    }

    return path_edges.size() == path_edges_total.size();

}

int get_eid_channel_id_ok(int nid0, int nid1, int channel_id) {
    int res = -1;
    for (auto& eid : nodes[nid0].eids) {
        auto& edge = edges[eid];
        if (edge.channel_used[channel_id]) continue;
        int nid_neighbor = edge.nid0 + edge.nid1 - nid0;
        if (nid_neighbor != nid1) continue;
        res = eid;
        break;
    }
    return res;
}

int get_channel_id_empty_max(vector<vector<int>>& path_edges) {

    //  记录每个通道的空闲数目, 注意: 在一对节点之间, 一个 通道ID 只能被增加一次
    vector<int> counts_empty_channel(num_of_channels, 0), indices_channel(num_of_channels);
    for (int i = 0; i < num_of_channels; ++i) {
        indices_channel[i] = i;
    }

    for (auto& e_ids : path_edges) {
        vector<bool> is_empty(num_of_channels, false);
        for (auto& e_id : e_ids) {
            auto& edge = edges[e_id];
            for (int i = 0; i < num_of_channels; ++i) {
                if (is_empty[i]) continue;// 已找到空位, 不用再寻找
                if (!edge.channel_used[i]) {
                    is_empty[i] = true;// 本条边的本通道可供使用
                }
            }
        }
        for (int i = 0; i < num_of_channels; ++i) {
            if (is_empty[i]) {
                ++counts_empty_channel[i];
            }
        }
    }

    sort(indices_channel.begin(), indices_channel.end(), [&](int& a, int& b) {
       return counts_empty_channel[a] > counts_empty_channel[b];
    });

    return indices_channel.front();

}

int add_edge(int nid0, int nid1) {

    Edge edge;
    edge.id = edges.size();
    edge.nid0 = nid0;
    edge.nid1 = nid1;
    edge.channel_used = vector<bool>(num_of_channels, false);
    edge.distance = dists_min[edge.nid0][edge.nid1];
    nodes[edge.nid0].eids.push_back(edge.id);
    nodes[edge.nid1].eids.push_back(edge.id);
    edges.emplace_back(edge);

    return edge.id;

}

vector<vector<int>> get_path_edges_total_via_nodes(vector<int>& path_nodes) {
    vector<vector<int>> path_edges_total;
    // 把节点形式的路径转换为边形式的路径
    for (int i = 1; i < path_nodes.size(); ++i) {
        int nid_prev = path_nodes[i - 1], nid_curr = path_nodes[i];
        vector<int> path_edges_tmp;
        auto& node_prev = nodes[nid_prev];
        for (auto& eid : node_prev.eids) {
            auto& edge = edges[eid];
            if (edge.nid0 + edge.nid1 - nid_prev == nid_curr) {
                path_edges_tmp.push_back(eid);
            }
        }
        path_edges_total.emplace_back(path_edges_tmp);
    }
    if (path_edges_total.size() + 1 != path_nodes.size()) {
#ifdef __APPLE__
        cerr << "path_edges_total.size() + 1 != path_nodes.size()" << endl;
        cerr << "path: ";
        for (int i = 0; i < path_nodes.size(); ++i) {
            cerr << path_nodes[i];
            if (i + 1 == path_nodes.size()) cerr << "\n";
            else cerr << "->";
        }
#endif
    }
    return path_edges_total;
}

void sort_operations() {
    // 对业务进行排序
    // 1. 把路径可选择的少的排前面
    // 2. 把长度长的排在前面
    vector<int> path_min(num_of_operations, INT_MAX), path_len(num_of_operations);

    for (int op_id = 0; op_id < num_of_operations; ++op_id) {
        vector<int> path_nodes;// 这条路径通过的节点
        auto& op = operations[op_id];
        int nid_source = op.nid0, nid_target = op.nid1;
        if (paths_of_nodes[nid_source].count(nid_target)) {
            path_nodes = paths_of_nodes[nid_source][nid_target];
        } else {
            // 初始化
            // 使用广搜得到可以到达的路径(不区分通道是否被使用完)
            path_nodes = bfs_find_path_nodes(nid_source, nid_target);
            paths_of_nodes[nid_source][nid_target] = path_nodes;
            auto tmp = path_nodes;
            reverse(tmp.begin(), tmp.end());
            paths_of_nodes[nid_target][nid_source] = tmp;
        }
        path_len[op_id] = path_nodes.size();
        for (int i = 1; i < path_nodes.size(); ++i) {
            int nid_prev = path_nodes[i - 1];
            int nid_curr = path_nodes[i];
            int cnt = 0;
            for (auto& eid : nodes[nid_prev].eids) {
                auto& edge = edges[eid];
                int nid_neighbor = edge.nid0 + edge.nid1 - nid_prev;
                if (nid_neighbor == nid_curr) {
                    ++cnt;
                }
            }
            path_min[op_id] = min(path_min[op_id], cnt);
        }
    }

    sort(operations.begin(), operations.end(), [&](Operation& op0, Operation& op1) {
        if (path_min[op0.id] != path_min[op1.id]) {
            return path_min[op0.id] < path_min[op1.id];// 把路径可选项更小的放在前面
        }
        return path_len[op0.id] > path_len[op1.id];// 把路径长的业务放在前面
    });


}