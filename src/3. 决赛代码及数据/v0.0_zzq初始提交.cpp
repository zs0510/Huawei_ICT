//#pragma GCC optimize(2)
#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>
#include <string>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include <limits.h>

#include <ctime>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;

/****************************************************以下是基本变量定义******************************************************/ 
static const int INF = 1e9 + 7;
static const int process_time = 3;           // 设定程序运行的最长时间 MAX_TIME / 2
static bool test = false;

struct Edge {
    int id;                                  // 边ID
    int distance;                            // 边长度
    int nid0;                                // 边的一个端点
    int nid1;                                // 边的另一个端点
    int count_channel_used = 0;
    bool valid = true;                      // true 表示有效, false 表示无效(被删除)
    vector<int> channel2operation;          // channel2operation[i] = j, j = -1 表示 通道i 未使用, 否则表示 通道i 被 业务j 使用
    vector<int> channel2path;               // channel2operation[i] = j, j = -1 表示 通道i 未使用, 否则表示 通道i 被 路径j 使用
    Edge() {}
    Edge(int d, int _n0, int _n1): distance(d), nid0(_n0), nid1(_n1) {}
};

struct Node {
    int id;                                 // 节点ID
    vector<int> eids;                       // 节点的临接边
    Node() {}
};

struct Operation {
    int id;                                 // 业务ID
    int nid0;                               // 业务的源节点
    int nid1;                               // 业务的目标节点
    int n_paths;                            // 本业务的路径数目
    unordered_set<int> channel_id;          // 本业务使用的通道编号
    unordered_set<int> edges_used;          // 已经被这条业务使用过的边
    vector<int> paths;                      // 这条业务的路径
    Operation() {}
};

struct Path {
    int id;
    int nid0;
    int nid1;
    int operation_id;                      // 这个路径所属于的业务
    int channel_id = -1;
    vector<int> nodes;
    vector<int> edges;
    vector<int> amplifiers;

};

/****************************************************以下是全局变量的定义******************************************************/ 
stringstream scout;
int num_of_nodes;
int num_of_edges_initial;
int num_of_operations;
int num_of_paths;
int num_of_channels;
int distance_attenuation_max;
std::chrono::steady_clock::time_point time_begin;
vector<Node> nodes;                                                     // num_of_nodes 是节点个数, ID 从[0, num_of_nodes-1]
vector<Edge> edges;                                                     // num_of_edges 是边个数, ID 从[0, num_of_edges-1]
vector<Operation> operations;                                           // num_of_operations 是业务个数, ID 从[0, num_of_operations-1]
vector<Path> paths;                                                     // num_of_paths 是路径个数, ID 从[0, num_of_paths-1]
unordered_map<int, unordered_map<int, int>> dists_min;                  // 记录两个节点之间直接连接的最短距离: 新生成的边的长度
unordered_map<int, unordered_map<int, vector<int>>> paths_of_nodes;     // 以节点序列记录两个节点之间的路径

/****************************************************以下是函数声明******************************************************/ 
vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id, unordered_set<int>& baned_nodes);   // 返回源点与目标点之间的节点路径(channel_id 为 -1 时不考虑通道占用情况)
vector<vector<int>> get_path_edges_total_via_nodes(vector<int>& path_nodes);                                        // 返回节点序列中所有可能的边
pair<int, vector<int>> find_cheapest_path_cid_nodes(int nid_source, int nid_target, int cost_base);                 // 在已知必须要加边的情况下, 选出加边最少的那条路径
pair<int, vector<int>> find_cheapest_path_nodes(int nid_source,
                                                int nid_target,
                                                int cid,
                                                int cost_cheapest,
                                                unordered_set<int>& baned_edges);

pair<int, int> get_channel_id_empty_max__cost_cheapest(vector<vector<int>>& path_edges);
int get_eid_channel_id_ok(int nid0, int nid1, int channel_id, unordered_set<int>& baned_edges);                     // 检查节点(直接相连)之间能否通过这个通道进行连接
int add_edge(int nid0, int nid1);                                                                                   // 新增两个节点之间的边
void sort_data();
void place_amplifier(Path& op);
void optimization_select_shortest_edge();
void optimization_transfer_operation();
void optimization_path_rebuild();
void clear_data();
bool is_channel_id_ok(vector<vector<int>>& path_edges_total, int channel_id, vector<int>& path_edges);              // 检查这些边的通道 pid 是否被占用
bool pave();                                                                                                        // 铺设线路
void print_map();
void tcout(stringstream& scout);

static inline bool timeout() {                                                                           // 记录从程序起始到现在过去了多少时间
    using namespace std::chrono;
    std::chrono::steady_clock::time_point time_current = std::chrono::steady_clock::now();
    duration<double> delta = duration_cast<duration<double>>(time_current - time_begin);
    return delta.count() >= process_time;
}

/****************************************************以下是主处理函数******************************************************/ 
void init_data() {
    // 输入
    cin >> num_of_nodes >> num_of_edges_initial >> num_of_operations >> num_of_paths >> num_of_channels >> distance_attenuation_max;
    nodes.resize(num_of_nodes);
    edges.resize(num_of_edges_initial);
    operations.resize(num_of_operations);
    paths.resize(num_of_paths);
    for (int i = 0; i < num_of_nodes; ++i) {
        nodes[i].id = i;
    }
    for (int i = 0; i < num_of_edges_initial; ++i) {// 读入 num_of_edges 行边
        cin >> edges[i].nid0 >> edges[i].nid1 >> edges[i].distance;
        edges[i].id = i;
        edges[i].channel2operation = vector<int>(num_of_channels, -1);
        edges[i].channel2path = vector<int>(num_of_channels, -1);
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

    int path_id = 0;
    for (int i = 0; i < num_of_operations; ++i) {// 读入 num_of_operations 行业务
        cin >> operations[i].nid0 >> operations[i].nid1 >> operations[i].n_paths;
        operations[i].id = i;
        for (int j = 0; j < operations[i].n_paths; ++j) {
            operations[i].paths.push_back(path_id);
            paths[path_id].id = path_id;
            paths[path_id].nid0 = operations[i].nid0;
            paths[path_id].nid1 = operations[i].nid1;
            paths[path_id].operation_id = operations[i].id;
            ++path_id;
        }
    }

    print_map();
}

/**
 * function: handle_data
 * param: null
 * input: null
 * Description: 路径铺设主入口程序，最重要的处理程序
 * ouput: null
 */
void handle_data() {
    // 处理
    auto nodes_backup = nodes;
    auto edges_backup = edges;
    auto operations_backup = operations;
    auto paths_backup = paths;
    auto paths_best = paths;
    auto edges_best = edges;
    int edges_new_smallest = INF;

    while (true) {
        nodes = nodes_backup;
        edges = edges_backup;
        operations = operations_backup;
        paths = paths_backup;

        sort_data();// 对业务进行排序

        pave();// 铺设线路

        int edges_new = edges.size() - num_of_edges_initial;
//        cout << "\nedges_new = " << edges_new << endl;
        if (edges_new < edges_new_smallest) {
            edges_new_smallest = edges_new;
            paths_best = paths;
            edges_best = edges;
        }
        // 程序达到设定运行时间，退出
        if (timeout()) {
            break;
        }
//        break;
    }

    paths = paths_best;
    edges = edges_best;

    // optimization_transfer_operation();// 将通道使用率低的边上的业务转移到其它边上从而移除掉该边, 但目前未奏效

//    optimization_path_rebuild();// 在保证可行性的前提下重新铺设线路, 这个能降低成本, 但偶尔会超时? 74847907 -> 74785404

    {
        // 放大器的优化
        // 初始化: 放置放大器
        for (auto& op : operations) {
            for (auto& pid : op.paths) {
                place_amplifier(paths[pid]);
            }
        }

        // 后处理: 在只更换边的情况下, 检查是否存在路径长度更短的连接方案
//        optimization_select_shortest_edge();
    }
}

// 输出数据处理完的结果
void output_data() {
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
        for (auto& pid : operation.paths) {
            auto& path = paths[pid];
            // 通道ID, 边数目, 放大器数目
            cout << path.channel_id << " " << path.edges.size() << " " << path.amplifiers.size();
            for (auto& eid : path.edges) {
                cout << " " << eid;
            }
            for (auto& nid : path.amplifiers) {
                cout << " " << nid;
            }
            cout << endl;
        }
    }
}

/****************************************************以下是主函数入口******************************************************/ 
int main() {
    // 记录起始时间
    time_begin = std::chrono::steady_clock::now();

    init_data();
    handle_data();
    output_data();

    return 0;
}

/****************************************************以下是打印函数的定义******************************************************/ 
void print_map() {
    // 输出建图信息
    scout << "\n建图信息如下: \n";
    for (int i = 0; i < num_of_nodes; ++i) {
        auto& node = nodes[i];
        scout << "Node ID: " << node.id << ", ";
        for (auto& eid : node.eids) {
            auto& edge = edges[eid];
            scout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
        }
//        for (int j = 0; j < node.eids.size(); ++j) {
//            auto& edge = edges[node.eids[j]];
//            cout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
//        }
        scout << "\n";
    }
    scout << "建图信息输出完毕\n";
    scout << "\n业务信息如下: \n";
    for (int pid = 0; pid < num_of_paths; ++pid) {
        auto& path = paths[pid];
        scout << "path.id = " << path.id << ", nid0 = " << path.nid0 << ", nid1 = " << path.nid1 << "\n";
    }
    tcout(scout);
}

void tcout(stringstream& scout) {
    if (test) {
        cout << scout.str() << '\n';
    }
    // 清空字符串流
    scout.str("");
}

/****************************************************以下是主要搜索函数的定义******************************************************/ 
vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id, unordered_set<int>& baned_edges) {
    // TODO: 使用双向搜索来加速 https://oi-wiki.org/search/bidirectional/
    // 当能不添加新边, 能通过当前通道连接这两个节点时, 返回经过的节点序列(从源点到目标点)
    // 否则返回一个空序列
    // channel_id == -1 表示不考虑通道编号
    unordered_map<int, int> prev;
    queue<int> que;
    que.push(nid_source);
    prev[nid_source] = -1;
    while (!que.empty() && !prev.count(nid_target)) {
        int nid = que.front();
        que.pop();
        for (auto& eid : nodes[nid].eids) {
            auto& edge = edges[eid];
            if (baned_edges.count(eid)) continue;// 已被本业务其它路径使用
            if (!edge.valid) continue;// 被标记为删除的边
            if (channel_id != -1 && edge.channel2operation[channel_id] != -1) continue;// 考虑通道使用情况
            int nid_next = edge.nid0 + edge.nid1 - nid;
            if (prev.count(nid_next)) {
                // TODO: 选择度更大的前驱节点使用? 如何避免出现环?
                continue;
            }
            que.push(nid_next);
            prev[nid_next] = nid;
        }
    }

    if (!prev.count(nid_target)) {     // 无法在当前通道编号下连接 nid_source 与 nid_target
        return {};
    }

    vector<int> path_nodes;
    int nid = nid_target;
    while (nid != -1) {
        path_nodes.push_back(nid);
        nid = prev[nid];
    }

    reverse(path_nodes.begin(), path_nodes.end());

    return path_nodes;
}

bool is_channel_id_ok(vector<vector<int>>& path_edges_total, int channel_id, vector<int>& path_edges) {
    // 当能通过当前通道铺设线路时, 返回 true, 并且 path_edges 不空
    // 否则, 返回 false, path_edges 为空
    path_edges.clear();
    int i = 0;
    while (i < path_edges_total.size()) {
        // 记录节点之间的边本通道是否为空(未使用)
        bool is_empty = false;                      
        for (auto& eid : path_edges_total[i]) {
            auto& edge = edges[eid];
            // 被标记为删除的边
            if (!edge.valid) continue;              
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

int get_eid_channel_id_ok(int nid0, int nid1, int channel_id, unordered_set<int>& baned_edges) {
    // 获取两个节点之间的边, 该边的 channel_id 未被使用
    int res = -1;
    for (auto& eid : nodes[nid0].eids) {
        auto& edge = edges[eid];
        if (baned_edges.count(eid)) continue;
        // 被标记为删除的边
        if (!edge.valid) continue;
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
            // 被标记为删除的边
            if (!edge.valid) continue;
            for (int i = 0; i < num_of_channels; ++i) {
                // 已找到空位, 不用再寻找
                if (is_empty[i]) continue;
                if (edge.channel2operation[i] == -1) {
                    // 本条边的本通道可供使用
                    is_empty[i] = true;
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
    // 返回加边最少的 通道编号 以及 加边数目
    return make_pair(indices_channel.front(), path_edges.size() - counts_empty_channel[indices_channel.front()]);

}

int add_edge(int nid0, int nid1) {

    Edge edge;
    edge.id = edges.size();
    edge.nid0 = nid0;
    edge.nid1 = nid1;
    edge.valid = true;
    edge.channel2operation = vector<int>(num_of_channels, -1);
    edge.channel2path = vector<int>(num_of_channels, -1);
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

void clear_data() {
//    for (auto& op : operations) {
//        op.channel_id = -1;
//        op.nodes.clear();
//        op.edges.clear();
//        op.amplifiers.clear();
//    }
//    for (auto& edge : edges) {
//        edge.count_channel_used = 0;
//        edge.valid = true;
//        edge.channel2operation = vector<int>(num_of_channels, -1);
//    }
//    paths_of_nodes.clear();

}

void sort_data() {

    // 对每个节点的邻接边进行排序, 把距离短的放在前面
    // TODO: 关键边? 或许应该把不怎么被使用的边放在前面, 把被很多业务使用的边放在后面
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }

    // 对业务进行排序

    //  随机排序
   random_device rd;
   mt19937 g(rd());

   sort(operations.begin(), operations.end(), [&](Operation& op1, Operation& op2) {
       if (op2.n_paths == 2) return false;
       if (op1.n_paths == 2) return true;
       return false;
   });

   int index_of_last_two = -1;
   for (int i = 0; i < num_of_operations; ++i) {
       if (operations[i].n_paths == 2) {
           index_of_last_two = i;
       } else if (operations[i].n_paths != 2) {
           break;
       }
   }

   if (index_of_last_two == -1) {
       shuffle(operations.begin(), operations.end(), g);
   } else {
       shuffle(operations.begin(), operations.begin() + index_of_last_two + 1, g);
       shuffle(operations.begin() + index_of_last_two + 1, operations.end(), g);
   }

#ifdef __APPLE__
   for (auto& op : operations) {
       cout << op.n_paths << "->";
   }
   cout << endl;
#endif

   return;

}

void place_amplifier(Path& path) {

    auto& path_nodes = path.nodes;
    auto& path_edges = path.edges;

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
    path.amplifiers = path_amplifiers;

}

void optimization_select_shortest_edge() {

//    sort(operations.begin(), operations.end(), [&](Operation& a, Operation& b) {
//        return a.amplifiers.size() > b.amplifiers.size();
//    });
//    for (int nid = 0; nid < num_of_nodes; ++nid) {
//        auto& node = nodes[nid];
//        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
//            return edges[eid0].distance < edges[eid1].distance;
//        });
//    }
//    for (auto& op : operations) {
//        if (op.amplifiers.empty()) break;
//        // 记录当前的连接方案
//        vector<int> path_edges_backup = op.edges;
//        vector<int> path_amplifiers_backup = op.amplifiers;
//
//        auto& path_nodes = op.nodes;
//        auto& path_edges = op.edges;
//        for (int i = 1; i < path_nodes.size(); ++i) {
//            auto nid0 = path_nodes[i - 1], nid1 = path_nodes[i];
//            // 把每条边都换成最短的
//            for (auto& eid : nodes[nid0].eids) {
//                auto& edge = edges[eid];
//                if (!edge.valid) continue;
//                if (edge.nid0 + edge.nid1 - nid0 != nid1) continue;
//                if (edge.channel2operation[op.channel_id] != -1) continue;
//                if (edge.distance < edges[path_edges[i - 1]].distance) {
//                    path_edges[i - 1] = eid;
//                    break;
//                }
//            }
//        }
//        place_amplifier(op);
//        if (path_edges_backup.size() <= op.amplifiers.size()) {
//            // 不需要改变, 恢复
//            op.edges = path_edges_backup;
//            op.amplifiers = path_amplifiers_backup;
//        } else {
//            // 改变可以减少放大器的个数, 更新
//            // 先把旧边的通道释放出来, 再占用新边的通道
//            for (auto& eid : path_edges_backup) {
//                edges[eid].channel2operation[op.channel_id] = -1;
//                --edges[eid].count_channel_used;
//            }
//            for (auto& eid : path_edges) {
//                edges[eid].channel2operation[op.channel_id] = op.id;
//                ++edges[eid].count_channel_used;
//            }
//        }
//    }

}

void optimization_transfer_operation() {

//    // 优化: 以通道使用率低的边为基础, 把业务转移到其它边上来移除这条边
//    vector<int> edge_indices(edges.size());
//    for (int i = 0; i < edge_indices.size(); ++i) {
//        edge_indices[i] = i;
//    }
//    // 把通道使用率低的边放前面
//    sort (edge_indices.begin(), edge_indices.end(), [&](int& eid0, int& eid1) {
//        return edges[eid0].count_channel_used < edges[eid1].count_channel_used;
//    });
//    for (auto& eid : edge_indices) {
//        if (eid < num_of_edges_initial) continue;// 初始输入的边, 无法被删去的边
//        auto& edge = edges[eid];
//        for (int cid = 0; cid < num_of_channels && edge.count_channel_used != 0; ++cid) {
//            if (edge.channel2operation[cid] == -1) continue;// 未被使用的通道, 无需转移
//            auto op_id = edge.channel2operation[cid];
//            auto& op = operations[op_id];
//            auto& path_nodes = op.nodes;
//            auto& path_edges = op.edges;
//            bool can_move_to_other_edges = false;// 这个通道上的业务是否可以转移到别的边上?
//            for (int i = 1; i < path_nodes.size(); ++i) {
//                if (path_edges[i - 1] == eid) {// 找到本边在该业务上的位置
//                    unordered_set<int> baned_nodes(op.nodes.begin(), op.nodes.end());
//                    auto path_nodes_tmp_without_add_new_edge = bfs_find_path_nodes(path_nodes[i - 1], path_nodes[i], cid, baned_nodes);
//                    if (!path_nodes_tmp_without_add_new_edge.empty()) {
//                        vector<int> path_edges_tmp_without_add_new_edge;
//                        for (int j = 1; j < path_nodes_tmp_without_add_new_edge.size(); ++j) {
//                            int eid_tmp_without_add_new_edge = get_eid_channel_id_ok(path_nodes_tmp_without_add_new_edge[j - 1], path_nodes_tmp_without_add_new_edge[j], cid);
//                            if (eid_tmp_without_add_new_edge != -1) {
//                                path_edges_tmp_without_add_new_edge.push_back(eid_tmp_without_add_new_edge);
//                            } else break;
//                        }
//                        if (path_nodes_tmp_without_add_new_edge.size() - 1 == path_edges_tmp_without_add_new_edge.size()) {
//                            // 更新边, 插入顶点
//                            path_nodes.insert(path_nodes.begin() + i, path_nodes_tmp_without_add_new_edge.begin() + 1, path_nodes_tmp_without_add_new_edge.end() - 1);
//                            // 移除旧边, 插入新边
//                            path_edges[i - 1] = path_edges_tmp_without_add_new_edge[0];
//                            path_edges.insert(path_edges.begin() + i, path_edges_tmp_without_add_new_edge.begin() + 1, path_edges_tmp_without_add_new_edge.end());
//                            for (auto& eid_tmp : path_edges_tmp_without_add_new_edge) {
//                                auto& edge_tmp = edges[eid_tmp];
//                                ++edge_tmp.count_channel_used;
//                                edge_tmp.channel2operation[cid] = op_id;
//                            }
//                            can_move_to_other_edges = true;
//                        }
//
//                    }
//                    break;
//                }
//            }
//            if (!can_move_to_other_edges) break;
//            // 这条边上的业务能转移到别的边上
//            --edge.count_channel_used;
//            edge.channel2operation[cid] = -1;
//        }
//        if (edge.count_channel_used == 0) {
//            // 这条边可以被移除
//            edge.valid = false;
//        }
//    }
//
//    // 移除掉 未被使用的, 后续添加进来 的边
//
//    int eid_new = 0;
//    for (int eid = 0; eid < edges.size(); ++eid) {
//        auto& edge = edges[eid];
//        if (edge.valid) {
//            edge.id = eid_new;
//            edge_indices[eid] = eid_new;// 建立边 旧ID 与 新ID 的映射
//            ++eid_new;
//        } else {
//            edge_indices[eid] = edges.size();// 无效边: 将其ID置为一个较大值, 便于后续移除
//        }
//    }
//
//    // 移除掉旧边
//    sort(edges.begin(), edges.end(), [&](Edge& e0, Edge& e1) {
//        return e0.id < e1.id;
//    });
//    edges.resize(eid_new);// 把冗余的边自动删除
//
//    // 更新所有业务的边ID
//    for (auto& op : operations) {
//        for (auto& eid : op.edges) {
//            eid = edge_indices[eid];
//        }
//    }

}

void optimization_path_rebuild() {
    // 把所有节点的边按长度升序排
//    for (auto& node : nodes) {
//        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
//            return edges[eid0].distance < edges[eid1].distance;
//        });
//    }
//
//    for (auto& op : operations) {
//        place_amplifier(op);
//    }
//    sort(operations.begin(), operations.end(), [&] (Operation& op0, Operation& op1) {
//        return op0.amplifiers.size() > op1.amplifiers.size();
//    });
//
//    for (auto& op : operations) {
//        // 先暂时释放通道
//        for (auto& eid : op.edges) {
//            edges[eid].channel2operation[op.channel_id] = -1;
//        }
//        // 重新寻找一个可行解
//        int path_channel_id = -1;
//        vector<int> path_nodes_tmp_shortest, path_edges_tmp_shortest;// 记录最短的可行路径
//        for (int channel_id = 0; channel_id < num_of_channels; ++channel_id) {
//            // 检查是否可以使用 通道ID 去连接源点与目标点
//            auto path_nodes_tmp = bfs_find_path_nodes(op.nid0, op.nid1, channel_id);
//            if (!path_nodes_tmp.empty()) {
//                auto path_edges_total_tmp = get_path_edges_total_via_nodes(path_nodes_tmp);
//                vector<int> path_edges_tmp;
//                if (is_channel_id_ok(path_edges_total_tmp, channel_id, path_edges_tmp)) {
//                    // 确认没问题, 再统一修改
//                    if (path_channel_id == -1 || path_nodes_tmp.size() < path_nodes_tmp_shortest.size()) {
//                        path_channel_id = channel_id;
//                        path_nodes_tmp_shortest = path_nodes_tmp;
//                        path_edges_tmp_shortest = path_edges_tmp;
//                    }
//                }
//            }
//        }
//
//        if (path_channel_id != -1) {
//            // 切换到新通道
//            op.channel_id = path_channel_id;
//            op.nodes = path_nodes_tmp_shortest;
//            op.edges = path_edges_tmp_shortest;
//        }
//        // 更新占有通道信息
//        for (auto& eid : op.edges) {
//            edges[eid].channel2operation[op.channel_id] = op.id;
//        }
//
//    }

}

pair<int, vector<int>> find_cheapest_path_nodes(int nid_source,
                                                int nid_target,
                                                int cid,
                                                int cost_cheapest,
                                                unordered_set<int>& baned_edges) {

    // TODO: 在加边数量相同的情况下, 选择经过边的数量最小的路径
    // get<0> 表示代价
    unordered_map<int, int> prev_nid, prev_cost;
    // tuple<当前路径加边的数目, 经过的边数目, pair<当前节点ID, 前驱节点ID>>
    priority_queue<tuple<int, int, pair<int, int>>, vector<tuple<int, int, pair<int, int>>>, greater<tuple<int, int, pair<int, int>>>> heap;// 小顶堆
    heap.push(make_tuple(0, 0, make_pair(nid_source, -1)));
    while (!heap.empty()) {
        auto [cost_current, length_current, pair_nodes] = heap.top();
        heap.pop();
        if (cost_current > cost_cheapest) break;
        int nid_current = get<0>(pair_nodes);
        int nid_prev = get<1>(pair_nodes);
        if (prev_nid.count(nid_current)) continue;// 已找到最少的加边路径
        prev_nid[nid_current] = nid_prev;
        prev_cost[nid_current] = cost_current;
        if (nid_current == nid_target) {
            break;
        }
        auto& node = nodes[nid_current];
        for (auto& eid : node.eids) {
//            if (baned_edges.count(eid)) continue;// 这条边已被这个业务其它路径使用
            auto& edge = edges[eid];
            int nid_next = edge.nid0 + edge.nid1 - nid_current;
            if (prev_nid.count(nid_next)) continue;// 已找到最少的加边路径
            // 若被占用, 则代价+1, 否则代价不变
            int cost_next = cost_current + (edge.channel2operation[cid] != -1);
            if (baned_edges.count(eid)) ++cost_next;// 这条边被禁用, 代价 + 1
            if (cost_next <= cost_cheapest) {
                heap.push(make_tuple(cost_next, length_current + 1, make_pair(nid_next, nid_current)));
            }

        }

    }

    if (!prev_cost.count(nid_target)) {
        return make_pair(-1, vector<int>());
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

bool pave() {

    for (auto& op : operations) {
        for (int i = 0; i < op.paths.size(); ++i) {
            auto& path = paths[op.paths[i]];
            int nid_source = path.nid0, nid_target = path.nid1;
            int path_channel_id = -1;
            vector<int> path_nodes;
            vector<int> path_edges;
            if (i == 1 && op.paths.size() == 2) {// 硬约束: 这条路径必须得与前一条路径使用的通道相同
                path_channel_id = (*op.channel_id.begin());
                vector<int> path_nodes_tmp = bfs_find_path_nodes(nid_source, nid_target, path_channel_id, op.edges_used);
                if (path_nodes_tmp.empty()) {
                    // 需要加边来满足要求
#ifdef __APPLE__
                    cout << "未找到相同通道的可用线路\n";
#endif
                    auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, path_channel_id, INF, op.edges_used);
                    path_nodes = get<1>(pair_cost_nodes);
                } else {
#ifdef __APPLE__
                    cout << "找到了相同通道的可用线路\n";
#endif
                    path_nodes = path_nodes_tmp;
                }

            } else {

                unordered_set<int> baned_edges_emtpy;
                path_nodes = bfs_find_path_nodes(nid_source, nid_target, -1, baned_edges_emtpy);
                // 得到一个基础解, 用于剪枝
                path_channel_id = rand() % num_of_channels;
                int cost_cheapest = 0;
                for (int i = 1; i < path_nodes.size(); ++i) {
                    // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
                    int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id, op.edges_used);
                    if (eid_old == -1) {
                        ++cost_cheapest;
                    }
                }
                int cid = rand() % num_of_channels;
                auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, cost_cheapest, op.edges_used);
                auto cost_tmp = get<0>(pair_cost_nodes);
                auto path_nodes_tmp = get<1>(pair_cost_nodes);
                if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
                    cost_cheapest = cost_tmp;
                    path_nodes = path_nodes_tmp;
                    path_channel_id = cid;
                }


//                for (int cid = 0; cid < num_of_channels; ++cid) {
//                    vector<int> path_nodes_tmp = bfs_find_path_nodes(nid_source, nid_target, cid, op.edges_used);
//                    if (!path_nodes_tmp.empty()) {
//                        // 找到第一条可用的线路
//                        path_nodes = path_nodes_tmp;
//                        path_channel_id = cid;
//                        break;
//                    }
//                }

//                if (path_channel_id == -1) {
//                    unordered_set<int> baned_edges_emtpy;
//                    path_nodes = bfs_find_path_nodes(nid_source, nid_target, -1, baned_edges_emtpy);
//
//                    path_channel_id = 0;
//                    int num_of_edges_smallest = 9999999999;
//                    for (int cid = 0; cid < num_of_channels; ++cid) {
//                        int num_of_new_edges = 0;
//                        for (int i = 1; i < path_nodes.size(); ++i) {
//                            // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
//                            int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id, op.edges_used);
//                            if (eid_old == -1) {
//                                ++num_of_new_edges;
//                            }
//                        }
//                        if (num_of_new_edges < num_of_edges_smallest) {
//                            path_channel_id = cid;
//                            num_of_edges_smallest = num_of_new_edges;
//                        }
//                    }
//                }
//
//                if (path_channel_id == -1) {
//                    int cost_cheapest = 9999999999;
//                    for (int cid = 0; cid < num_of_channels; ++cid) {
//                        auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, op.edges_used);
//                        auto cost_tmp = get<0>(pair_cost_nodes);
//                        auto path_nodes_tmp = get<1>(pair_cost_nodes);
//                        if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
//                            path_channel_id = cid;
//                            path_nodes = path_nodes_tmp;
//                            cost_cheapest = cost_tmp;
//                        }
//                    }
//                }

//                if (path_channel_id == -1) {
//                    int cost_cheapest = 9999999999;
//                    int cid = rand() % num_of_channels;
//                    for (int cid = 0; cid < num_of_channels; ++cid) {
//                        auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, op.edges_used);
//                        auto cost_tmp = get<0>(pair_cost_nodes);
//                        auto path_nodes_tmp = get<1>(pair_cost_nodes);
//                        if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
//                            path_channel_id = cid;
//                            path_nodes = path_nodes_tmp;
//                            cost_cheapest = cost_tmp;
//                        }
//                    }
//                }

            }


//            unordered_set<int> baned_edges_emtpy;
//            path_nodes = bfs_find_path_nodes(nid_source, nid_target, -1, baned_edges_emtpy);
//            if (i == 1 && op.paths.size() == 2) {
//                path_channel_id = (*op.channel_id.begin());
//            } else {
//                path_channel_id = rand() % num_of_channels;
//            }


            // 针对通道被占用的边, 新建边来满足业务需求
            for (int i = 1; i < path_nodes.size(); ++i) {
                // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
                int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id, op.edges_used);
                if (eid_old != -1) {
                    path_edges.emplace_back(eid_old);
                    continue;
                }
                int eid_new = add_edge(path_nodes[i - 1], path_nodes[i]);
                path_edges.emplace_back(eid_new);
            }

            // 请保证 path_nodes, path_edges 及 path_channel_id 是可行的
            for (auto& eid : path_edges) {
                auto& edge = edges[eid];
#ifdef __APPLE__
                if (edge.channel2operation[path_channel_id] != -1) {
                    cerr << "Error: Reusing the same channel\n";
                    continue;
                }
#endif
                edge.channel2operation[path_channel_id] = op.id;
                edge.channel2path[path_channel_id] = path.id;
                ++edge.count_channel_used;
            }

#ifdef __APPLE__
            cout << "path.id = " << path.id << ": ";
            for (int i = 0; i < path_nodes.size(); ++i) {
                cout << path_nodes[i];
                if (i + 1 != path_nodes.size()) cout << "->";
            }
            cout << "\n";
#endif

            path.edges = path_edges;
            path.nodes = path_nodes;
            path.channel_id = path_channel_id;

            op.edges_used.insert(path_edges.begin(), path_edges.end());// 记录这条业务使用的边
            op.channel_id.insert(path_channel_id);
        }

    }

    return true;

}