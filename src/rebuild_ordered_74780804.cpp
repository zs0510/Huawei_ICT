#pragma GCC optimize(2)
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
#include <random>
#include <queue>

using namespace std;

const int edges_add_max = 20000;// 可以增加的最大边数

struct Edge {
    int id;// 边ID
    int distance;// 边长度
    int nid0;// 边的一个端点
    int nid1;// 边的另一个端点
    int count_channel_used = 0;
    bool valid = true;// true 表示有效, false 表示无效(被删除)
    vector<int> channel2operation;// channel2operation[i] = j, j = -1 表示 通道i 未使用, 否则表示 通道i 被 业务j 使用
    Edge() {}
    Edge(int d, int _n0, int _n1): distance(d), nid0(_n0), nid1(_n1) {}
};

struct Node {
    int id;// 节点ID
    vector<int> eids;// 节点的临接边
    Node() {}
};

struct Operation {
    int id;// 业务ID
    int nid0;// 业务的源节点
    int nid1;// 业务的目标节点
    int channel_id = -1;// 业务使用的通道编号
    vector<int> nodes;// 业务经过的节点序列ID
    vector<int> edges;// 业务经过的边序列ID
    vector<int> amplifiers;// 放置放大器
    Operation() {}
};

int num_of_nodes, num_of_edges_initial, num_of_operations, num_of_channels, distance_attenuation_max;
vector<Node> nodes;// num_of_nodes 是节点个数, ID 从[0, num_of_nodes-1]
vector<Edge> edges;// num_of_edges 是边个数, ID 从[0, num_of_edges-1]
vector<Operation> operations;// num_of_operations 是业务个数, ID 从[0, num_of_operations-1]
unordered_map<int, unordered_map<int, int>> dists_min;// 记录两个节点之间直接连接的最短距离: 新生成的边的长度
unordered_map<int, unordered_map<int, vector<int>>> paths_of_nodes;// 以节点序列记录两个节点之间的路径

vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id = -1, unordered_set<int> baned_nodes = {});// 返回源点与目标点之间的节点路径(channel_id 为 -1 时不考虑通道占用情况)
void dfs_find_path_edges(const int& nid_target,
                         int nid_current,
                         vector<int>& channel_empty,
                         vector<int>& path_edges,
                         unordered_set<int>& visited,
                         int& count_added_min,
                         int& channel_id_min,
                         vector<int>& path_edges_min);
pair<int, vector<int>> find_cheapest_path_cid_nodes(int nid_source, int nid_target, int cost_base);// 在已知必须要加边的情况下, 选出加边最少的那条路径
pair<int, vector<int>> find_cheapest_path_nodes(int nid_source, int nid_target, int cid, int cost_cheapest);
bool is_channel_id_ok(vector<vector<int>>& path_edges_total, int channel_id, vector<int>& path_edges);// 检查这些边的通道 pid 是否被占用
int get_eid_channel_id_ok(int nid0, int nid1, int channel_id);// 检查节点(直接相连)之间能否通过这个通道进行连接
pair<int, int> get_channel_id_empty_max__cost_cheapest(vector<vector<int>>& path_edges);
int add_edge(int nid0, int nid1);// 新增两个节点之间的边
vector<vector<int>> get_path_edges_total_via_nodes(vector<int>& path_nodes);// 返回节点序列中所有可能的边
void sort_operations();
void place_amplifier(Operation& op);
void optimization_select_shortest_edge();
void optimization_transfer_operation();
void optimization_path_rebuild();


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
        edges[i].channel2operation = vector<int>(num_of_channels, -1);
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
        for (auto& eid : node.eids) {
            auto& edge = edges[eid];
            cout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
        }
//        for (int j = 0; j < node.eids.size(); ++j) {
//            auto& edge = edges[node.eids[j]];
//            cout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
//        }
        cout << "\n";
    }
#endif

    // 对每个节点的邻接边进行排序, 把距离短的放在前面
    // TODO: 关键边? 或许应该把不怎么被使用的边放在前面, 把被很多业务使用的边放在后面
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }

    // 对业务进行排序

    sort_operations();

    //  随机排序
//    random_device rd;
//    mt19937 g(rd());
//    shuffle(operations.begin(), operations.end(), g);

    // 自定义随机种子排序
//    srand(1234567);
//    random_shuffle(operations.begin(), operations.end());

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

        int path_channel_id = -1;// 本次业务使用的通道 ID
        vector<int> path_edges;// 这条路径通过的边

        // 检查最优节点个数的路径是否可达
        vector<vector<int>> path_edges_total = get_path_edges_total_via_nodes(path_nodes);// 节点之间的所有的边
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
            vector<int> path_nodes_tmp_shortest, path_edges_tmp_shortest;// 记录最短的可行路径
            for (int channel_id = 0; channel_id < num_of_channels; ++channel_id) {
                // 检查是否可以使用 通道ID 去连接源点与目标点
                auto path_nodes_tmp = bfs_find_path_nodes(nid_source, nid_target, channel_id);
                if (!path_nodes_tmp.empty()) {
                    auto path_edges_total_tmp = get_path_edges_total_via_nodes(path_nodes_tmp);
                    vector<int> path_edges_tmp;
                    if (is_channel_id_ok(path_edges_total_tmp, channel_id, path_edges_tmp)) {
                        // 确认没问题, 再统一修改
                        if (path_channel_id == -1 || path_nodes_tmp.size() < path_nodes_tmp_shortest.size()) {
                            path_channel_id = channel_id;
                            path_nodes_tmp_shortest = path_nodes_tmp;
                            path_edges_tmp_shortest = path_edges_tmp;
                        }
                    }
                }
            }

            if (path_channel_id != -1) {
                // 在不新增边的前提下, 找到了可达路径
                path_nodes = path_nodes_tmp_shortest;
                path_edges = path_edges_tmp_shortest;

            } else {
                // solution 2: 无法通过现有边完成业务时, 新增边以满足要求
                // TODO: 不简单地仅仅依赖初始遍历得到的节点序列, 选择一个加边更少的节点序列?
                // 重新生成一个节点序列, 这个节点序列可以通过添加较少的边完成需求

                // 统计每条边的通道使用情况, 计算当前节点序列需要添加的最少边个数, 用作基础解
                auto solution_base = get_channel_id_empty_max__cost_cheapest(path_edges_total);
                path_channel_id = solution_base.first;
//                // 尝试找出比基础解更优的解
//                auto res = find_cheapest_path_cid_nodes(nid_source, nid_target, solution_base.second);
//                if (!res.second.empty()) {// 找到了比基础解更优的解
//                    path_channel_id = res.first;
//                    path_nodes = res.second;
//                }

                // 搜索一个局部连接序列, 这个序列可以在不添加新边的前提下连接 path_nodes[i - 1], path_nodes[i]
                for (int i = 1; i < path_nodes.size(); ++i) {
                    int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id);
                    if (eid_old == -1) {
                        unordered_set<int> baned_nodes(path_nodes.begin(), path_nodes.end());
                        auto path_nodes_tmp_without_add_new_edge = bfs_find_path_nodes(path_nodes[i - 1], path_nodes[i], path_channel_id, baned_nodes);
                        if (!path_nodes_tmp_without_add_new_edge.empty()) {
                            path_nodes.insert(path_nodes.begin() + i, path_nodes_tmp_without_add_new_edge.begin() + 1, path_nodes_tmp_without_add_new_edge.end() - 1);
                            i += (path_nodes_tmp_without_add_new_edge.size() - 2);
                        }
                    }
                }
                // 针对通道被占用的边, 新建边来满足业务需求
                for (int i = 1; i < path_nodes.size(); ++i) {
                    // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
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
        // 请保证 path_nodes, path_edges 及 path_channel_id 是可行的
#ifdef __APPLE__
        cout << "success! nid_source = " << nid_source << ", nid_target = " << nid_target << ", path = ";
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
            if (edge.channel2operation[path_channel_id] != -1) {
                cerr << "Error: Reusing the same channel\n";
                continue;
            }
#endif
            edge.channel2operation[path_channel_id] = op.id;
            ++edge.count_channel_used;
        }

        op.channel_id = path_channel_id;
        op.nodes = path_nodes;
        op.edges = path_edges;


    }

//    optimization_transfer_operation();// 将通道使用率低的边上的业务转移到其它边上从而移除掉该边, 但目前未奏效

    optimization_path_rebuild();// 在保证可行性的前提下重新铺设线路, 这个能降低成本, 但偶尔会超时? 74847907 -> 74785404

    {
        // 放大器的优化
        // 初始化: 放置放大器
        for (auto& op : operations) {
            place_amplifier(op);
        }

        // 后处理: 在只更换边的情况下, 检查是否存在路径长度更短的连接方案
        optimization_select_shortest_edge();
    }


    // 输出
    // 输出加边信息
    cout << edges.size() - num_of_edges_initial << endl;// 新增边的个数
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

vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id, unordered_set<int> baned_nodes) {
    // TODO: 使用双向搜索来加速 https://oi-wiki.org/search/bidirectional/
    // 当能通过当前通道连接这两个节点时, 返回经过的节点序列(从源点到目标点)
    // 否则返回一个空序列
    // channel_id == -1 表示不考虑通道编号
    unordered_map<int, int> prev;
    queue<int> que;
    que.push(nid_source);
    prev[nid_source] = -1;
    baned_nodes.erase(nid_source);
    baned_nodes.erase(nid_target);
    while (!que.empty() && !prev.count(nid_target)) {
        int nid = que.front();
        que.pop();
        for (auto& eid : nodes[nid].eids) {
            auto& edge = edges[eid];
            if (!edge.valid) continue;// 被标记为删除的边
            if (channel_id != -1 && edge.channel2operation[channel_id] != -1) continue;// 考虑通道使用情况
            int nid_next = edge.nid0 + edge.nid1 - nid;
            if (prev.count(nid_next) || baned_nodes.count(nid_next)) {
                // TODO: 选择度更大的前驱节点使用? 如何避免出现环?
                continue;
            }
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
    // 当能通过当前通道铺设线路时, 返回 true, 并且 path_edges 不空
    // 否则, 返回 false, path_edges 为空
    path_edges.clear();
    int i = 0;
    while (i < path_edges_total.size()) {
        bool is_empty = false;// 记录节点之间的边本通道是否为空(未使用)
        for (auto& eid : path_edges_total[i]) {
            auto& edge = edges[eid];
            if (!edge.valid) continue;// 被标记为删除的边
            if (edge.channel2operation[channel_id] != -1) continue;
            is_empty = true;
            path_edges.emplace_back(eid);
            break;
        }
        if (!is_empty) break;
        ++i;
    }

    // 无法通过当前通道编号去铺设线路
    if (i != path_edges_total.size()) {
        path_edges.clear();
    }

    return path_edges.size() == path_edges_total.size();

}

int get_eid_channel_id_ok(int nid0, int nid1, int channel_id) {
    // 获取两个节点之间的边, 该边的 channel_id 未被使用
    int res = -1;
    for (auto& eid : nodes[nid0].eids) {
        auto& edge = edges[eid];
        if (!edge.valid) continue;// 被标记为删除的边
        if (edge.channel2operation[channel_id] != -1) continue;
        int nid_neighbor = edge.nid0 + edge.nid1 - nid0;
        if (nid_neighbor != nid1) continue;
        res = eid;
        break;
    }
    return res;
}

pair<int, int> get_channel_id_empty_max__cost_cheapest(vector<vector<int>>& path_edges) {

    //  记录每个通道的空闲数目, 注意: 在一对节点之间, 一个 通道ID 只能被增加一次
    vector<int> counts_empty_channel(num_of_channels, 0), indices_channel(num_of_channels);
    for (int i = 0; i < num_of_channels; ++i) {
        indices_channel[i] = i;
    }

    for (auto& e_ids : path_edges) {
        vector<bool> is_empty(num_of_channels, false);
        for (auto& e_id : e_ids) {
            auto& edge = edges[e_id];
            if (!edge.valid) continue;// 被标记为删除的边
            for (int i = 0; i < num_of_channels; ++i) {
                if (is_empty[i]) continue;// 已找到空位, 不用再寻找
                if (edge.channel2operation[i] == -1) {
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

    return make_pair(indices_channel.front(), counts_empty_channel[indices_channel.front()]);

}

int add_edge(int nid0, int nid1) {

    Edge edge;
    edge.id = edges.size();
    edge.nid0 = nid0;
    edge.nid1 = nid1;
    edge.valid = true;
    edge.channel2operation = vector<int>(num_of_channels, -1);
    edge.distance = dists_min[edge.nid0][edge.nid1];
    nodes[edge.nid0].eids.push_back(edge.id);
    nodes[edge.nid1].eids.push_back(edge.id);

    edges.emplace_back(edge);

//    sort(nodes[edge.nid0].eids.begin(), nodes[edge.nid0].eids.end(), [&](int& a, int& b) {
//        return edges[a].distance < edges[b].distance;
//    });
//
//    sort(nodes[edge.nid1].eids.begin(), nodes[edge.nid1].eids.end(), [&](int& a, int& b) {
//        return edges[a].distance < edges[b].distance;
//    });

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
            if (!edge.valid) continue;// 被标记为删除的边
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
    // path_min 记录路径的最窄处宽度, path_len 记录路径的最短长度
    vector<int> path_min(num_of_operations, INT_MAX), path_len(num_of_operations);
    vector<long long> path_multiplicative(num_of_operations, 1);

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
            int cnt = 0;// 记录这两个节点之间能通过多少边连接
            for (auto& eid : nodes[nid_prev].eids) {
                auto& edge = edges[eid];
                int nid_neighbor = edge.nid0 + edge.nid1 - nid_prev;
                if (nid_neighbor == nid_curr) {
                    ++cnt;
                }
            }
            path_min[op_id] = min(path_min[op_id], cnt);
            path_multiplicative[op_id] *= (long long)(cnt);
        }
    }

    sort(operations.begin(), operations.end(), [&](Operation& op0, Operation& op1) {
        if (path_min[op0.id] != path_min[op1.id]) {
            return path_min[op0.id] < path_min[op1.id];// 把路径可选项更小的放在前面
        }
        if (path_multiplicative[op0.id] != path_multiplicative[op1.id]) {
            return path_multiplicative[op0.id] < path_multiplicative[op1.id];
        }
        return path_len[op0.id] > path_len[op1.id];// 把路径长的业务放在前面
    });


}

void place_amplifier(Operation& op) {


    auto& path_nodes = op.nodes;
    auto& path_edges = op.edges;

    // 放置放大器
    vector<int> path_amplifiers;
    int dist_sum = edges[path_edges.front()].distance;
    for (int i = 1; i < path_edges.size(); ++i) {
        // dist_sum 表示当前到达这个节点的线路长度
        auto& edge = edges[path_edges[i]];
        if (dist_sum + edge.distance >= distance_attenuation_max) {
            path_amplifiers.push_back(path_nodes[i]);// 其实是这条边的左端点, 注意: 边个数与节点个数差 1
            dist_sum = edge.distance;
        } else {
            dist_sum += edge.distance;
        }
    }
    op.amplifiers = path_amplifiers;

}

void optimization_select_shortest_edge() {

    sort(operations.begin(), operations.end(), [&](Operation& a, Operation& b) {
        return a.amplifiers.size() > b.amplifiers.size();
    });
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }
    for (auto& op : operations) {
        if (op.amplifiers.empty()) break;
        // 记录当前的连接方案
        vector<int> path_edges_backup = op.edges;
        vector<int> path_amplifiers_backup = op.amplifiers;

        auto& path_nodes = op.nodes;
        auto& path_edges = op.edges;
        for (int i = 1; i < path_nodes.size(); ++i) {
            auto nid0 = path_nodes[i - 1], nid1 = path_nodes[i];
            // 把每条边都换成最短的
            for (auto& eid : nodes[nid0].eids) {
                auto& edge = edges[eid];
                if (!edge.valid) continue;
                if (edge.nid0 + edge.nid1 - nid0 != nid1) continue;
                if (edge.channel2operation[op.channel_id] != -1) continue;
                if (edge.distance < edges[path_edges[i - 1]].distance) {
                    path_edges[i - 1] = eid;
                    break;
                }
            }
        }
        place_amplifier(op);
        if (path_edges_backup.size() <= op.amplifiers.size()) {
            // 不需要改变, 恢复
            op.edges = path_edges_backup;
            op.amplifiers = path_amplifiers_backup;
        } else {
            // 改变可以减少放大器的个数, 更新
            // 先把旧边的通道释放出来, 再占用新边的通道
            for (auto& eid : path_edges_backup) {
                edges[eid].channel2operation[op.channel_id] = -1;
                --edges[eid].count_channel_used;
            }
            for (auto& eid : path_edges) {
                edges[eid].channel2operation[op.channel_id] = op.id;
                ++edges[eid].count_channel_used;
            }
        }
    }

}

void optimization_transfer_operation() {

    // 优化: 以通道使用率低的边为基础, 把业务转移到其它边上来移除这条边
    vector<int> edge_indices(edges.size());
    for (int i = 0; i < edge_indices.size(); ++i) {
        edge_indices[i] = i;
    }
    // 把通道使用率低的边放前面
    sort (edge_indices.begin(), edge_indices.end(), [&](int& eid0, int& eid1) {
        return edges[eid0].count_channel_used < edges[eid1].count_channel_used;
    });
    for (auto& eid : edge_indices) {
        if (eid < num_of_edges_initial) continue;// 初始输入的边, 无法被删去的边
        auto& edge = edges[eid];
        for (int cid = 0; cid < num_of_channels && edge.count_channel_used != 0; ++cid) {
            if (edge.channel2operation[cid] == -1) continue;// 未被使用的通道, 无需转移
            auto op_id = edge.channel2operation[cid];
            auto& op = operations[op_id];
            auto& path_nodes = op.nodes;
            auto& path_edges = op.edges;
            bool can_move_to_other_edges = false;// 这个通道上的业务是否可以转移到别的边上?
            for (int i = 1; i < path_nodes.size(); ++i) {
                if (path_edges[i - 1] == eid) {// 找到本边在该业务上的位置
                    unordered_set<int> baned_nodes(op.nodes.begin(), op.nodes.end());
                    auto path_nodes_tmp_without_add_new_edge = bfs_find_path_nodes(path_nodes[i - 1], path_nodes[i], cid, baned_nodes);
                    if (!path_nodes_tmp_without_add_new_edge.empty()) {
                        vector<int> path_edges_tmp_without_add_new_edge;
                        for (int j = 1; j < path_nodes_tmp_without_add_new_edge.size(); ++j) {
                            int eid_tmp_without_add_new_edge = get_eid_channel_id_ok(path_nodes_tmp_without_add_new_edge[j - 1], path_nodes_tmp_without_add_new_edge[j], cid);
                            if (eid_tmp_without_add_new_edge != -1) {
                                path_edges_tmp_without_add_new_edge.push_back(eid_tmp_without_add_new_edge);
                            } else break;
                        }
                        if (path_nodes_tmp_without_add_new_edge.size() - 1 == path_edges_tmp_without_add_new_edge.size()) {
                            // 更新边, 插入顶点
                            path_nodes.insert(path_nodes.begin() + i, path_nodes_tmp_without_add_new_edge.begin() + 1, path_nodes_tmp_without_add_new_edge.end() - 1);
                            // 移除旧边, 插入新边
                            path_edges[i - 1] = path_edges_tmp_without_add_new_edge[0];
                            path_edges.insert(path_edges.begin() + i, path_edges_tmp_without_add_new_edge.begin() + 1, path_edges_tmp_without_add_new_edge.end());
                            for (auto& eid_tmp : path_edges_tmp_without_add_new_edge) {
                                auto& edge_tmp = edges[eid_tmp];
                                ++edge_tmp.count_channel_used;
                                edge_tmp.channel2operation[cid] = op_id;
                            }
                            can_move_to_other_edges = true;
                        }

                    }
                    break;
                }
            }
            if (!can_move_to_other_edges) break;
            // 这条边上的业务能转移到别的边上
            --edge.count_channel_used;
            edge.channel2operation[cid] = -1;
        }
        if (edge.count_channel_used == 0) {
            // 这条边可以被移除
            edge.valid = false;
        }
    }

    // 移除掉 未被使用的, 后续添加进来 的边

    int eid_new = 0;
    for (int eid = 0; eid < edges.size(); ++eid) {
        auto& edge = edges[eid];
        if (edge.valid) {
            edge.id = eid_new;
            edge_indices[eid] = eid_new;// 建立边 旧ID 与 新ID 的映射
            ++eid_new;
        } else {
            edge_indices[eid] = edges.size();// 无效边: 将其ID置为一个较大值, 便于后续移除
        }
    }

    // 移除掉旧边
    sort(edges.begin(), edges.end(), [&](Edge& e0, Edge& e1) {
        return e0.id < e1.id;
    });
    edges.resize(eid_new);// 把冗余的边自动删除

    // 更新所有业务的边ID
    for (auto& op : operations) {
        for (auto& eid : op.edges) {
            eid = edge_indices[eid];
        }
    }

}

void dfs_find_path_edges(const int& nid_target,
                         int nid_current,
                         vector<int>& channel_empty,
                         vector<int>& path_edges,
                         unordered_set<int>& visited,
                         int& count_added_min,
                         int& channel_id_min,
                         vector<int>& path_edges_min) {

    //  计算当前连接情况下, 需要新增多少条边
    int count_added = path_edges.size() + 1, channel_id = -1;
    for (int cid = 0; cid < num_of_channels; ++cid) {
        if ((int)(path_edges.size()) - channel_empty[cid] < count_added) {
            count_added = (int)(path_edges.size()) - channel_empty[cid];
            channel_id = cid;
        }

    }

    if (count_added >= count_added_min) {
        return;
    }

    if (nid_current == nid_target) {
        count_added_min = count_added;
        channel_id_min = channel_id;
        path_edges_min = path_edges;
        return;
    }

    auto& node_current = nodes[nid_current];
    visited.insert(nid_current);

    for (auto& eid : node_current.eids) {
        auto& edge = edges[eid];
        int nid_next = edge.nid0 + edge.nid1 - nid_current;
        if (visited.count(nid_next)) continue;

        for (int cid = 0; cid < num_of_channels; ++cid) {
            if (edge.channel2operation[cid] == -1) {
                ++channel_empty[cid];
            }
        }
        path_edges.push_back(eid);

        dfs_find_path_edges(nid_target, nid_current, channel_empty, path_edges, visited,
                            count_added_min, channel_id_min, path_edges_min);

        for (int cid = 0; cid < num_of_channels; ++cid) {
            if (edge.channel2operation[cid] == -1) {
                --channel_empty[cid];
            }
        }
        path_edges.pop_back();

    }

    visited.erase(nid_current);


}

void optimization_path_rebuild() {
    // 把所有节点的边按长度升序排
    for (auto& node : nodes) {
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }

    for (auto& op : operations) {
        place_amplifier(op);
    }
    sort(operations.begin(), operations.end(), [&] (Operation& op0, Operation& op1) {
        return op0.amplifiers.size() > op1.amplifiers.size();
    });

    for (auto& op : operations) {
        // 先暂时释放通道
        for (auto& eid : op.edges) {
            edges[eid].channel2operation[op.channel_id] = -1;
        }
        // 重新寻找一个可行解
        int path_channel_id = -1;
        vector<int> path_nodes_tmp_shortest, path_edges_tmp_shortest;// 记录最短的可行路径
        for (int channel_id = 0; channel_id < num_of_channels; ++channel_id) {
            // 检查是否可以使用 通道ID 去连接源点与目标点
            auto path_nodes_tmp = bfs_find_path_nodes(op.nid0, op.nid1, channel_id);
            if (!path_nodes_tmp.empty()) {
                auto path_edges_total_tmp = get_path_edges_total_via_nodes(path_nodes_tmp);
                vector<int> path_edges_tmp;
                if (is_channel_id_ok(path_edges_total_tmp, channel_id, path_edges_tmp)) {
                    // 确认没问题, 再统一修改
                    if (path_channel_id == -1 || path_nodes_tmp.size() < path_nodes_tmp_shortest.size()) {
                        path_channel_id = channel_id;
                        path_nodes_tmp_shortest = path_nodes_tmp;
                        path_edges_tmp_shortest = path_edges_tmp;
                    }
                }
            }
        }

        if (path_channel_id != -1) {
            // 切换到新通道
            op.channel_id = path_channel_id;
            op.nodes = path_nodes_tmp_shortest;
            op.edges = path_edges_tmp_shortest;
        }
        // 更新占有通道信息
        for (auto& eid : op.edges) {
            edges[eid].channel2operation[op.channel_id] = op.id;
        }

    }

}

pair<int, vector<int>> find_cheapest_path_cid_nodes(int nid_source, int nid_target, int cost_base) {
    int cid_cheapest = -1, cost_cheapest = cost_base;
    vector<int> path_nodes_cheapest;
    for (int cid = 0; cid < num_of_channels; ++cid) {
        auto res = find_cheapest_path_nodes(nid_source, nid_target, cid, cost_cheapest);
        if (res.first < cost_cheapest) {
            cost_cheapest = res.first;
            path_nodes_cheapest = res.second;
            cid_cheapest = cid;
            break;// 加速, 找到一个就退出
        }
    }

    return make_pair(cid_cheapest, path_nodes_cheapest);
}

pair<int, vector<int>> find_cheapest_path_nodes(int nid_source, int nid_target, int cid, int cost_cheapest) {

    unordered_map<int, int> prev_nid, prev_cost;
    // tuple<代价, 当前节点ID, 前驱节点ID>
    priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>, greater<tuple<int, int, int>>> heap;// 小顶堆
    heap.push(make_tuple(0, nid_source, -1));
    while (!heap.empty()) {
        auto [cost_current, nid_current, nid_prev] = heap.top();
        heap.pop();
        prev_nid[nid_current] = nid_prev;
        prev_cost[nid_current] = cost_current;
        if (cost_current >= cost_cheapest || nid_current == nid_target) break;
        auto& node = nodes[nid_current];
        for (auto& eid : node.eids) {
            auto& edge = edges[eid];
            int nid_next = edge.nid0 + edge.nid1 - nid_current;
            if (prev_nid.count(nid_next)) continue;
            // 若被占用, 则代价+1, 否则代价不变
            int cost_next = cost_current + (edge.channel2operation[cid] != -1);
            heap.push(make_tuple(cost_next, nid_next, nid_current));
        }

    }

    if (!prev_cost.count(nid_target)) {
        return make_pair(cost_cheapest + 1, vector<int>());
    }

    vector<int> path_nodes;
    int nid = nid_target;
    while (nid != -1) {
        path_nodes.push_back(nid);
        nid = prev_nid[nid];
    }

    reverse(path_nodes.begin(), path_nodes.end());

    return make_pair(prev_cost[nid_target], path_nodes);

}