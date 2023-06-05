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

// TEST_HDU01 确定是否输出调试信息, 提交时请注释下行代码!
//#define TEST_HDU01

using namespace std;

/****************************************************以下是基本变量定义******************************************************/
static const int INF = 1e9 + 7;              // 取一个极大的质数, 便于取模
static const double process_time_max = 90.0;  // 设定程序运行的最长时间(秒)

class ChannelStatue {                        // 自定义128位的数字用于表示通道使用情况, 位为1表示空闲, 位为0表示在用
public:
    ChannelStatue() : high(0), low(0){}
    ChannelStatue(int32_t l) : high(-(l < 0)), low(l) {}
    ChannelStatue(int64_t l) : high(-(l < 0)), low(l) {}
    ChannelStatue(uint32_t l) : high(0), low(l) {}
    ChannelStatue(uint64_t l) : high(0), low(l) {}
    ChannelStatue(uint64_t h, uint64_t l) : high(h), low(l) {}

    bool     operator == (const ChannelStatue& rhs)const { return high == rhs.high && low == rhs.low; }
    bool     operator != (const ChannelStatue& rhs)const { return high != rhs.high || low != rhs.low; }
    bool     operator < (const ChannelStatue& rhs)const { return (high == rhs.high) ? low < rhs.low : high < rhs.high; }
    bool     operator < (const int64_t& rhs)const { return *this < ChannelStatue(rhs); }
    bool     operator !()const                    { return !(high != 0 || low != 0); }
    ChannelStatue  operator -()const                    { return ++ChannelStatue(~high, ~low); }
    ChannelStatue  operator ~()const                    { return ChannelStatue(~high, ~low); }

    ChannelStatue& operator++()    { high += (++low == 0); return *this; }
    ChannelStatue& operator--()    { high -= (low-- == 0); return *this; }
    ChannelStatue  operator++(int) { auto tmp = *this; ++(*this); return tmp; }
    ChannelStatue  operator--(int) { auto tmp = *this; --(*this); return tmp; }

    ChannelStatue& operator |= (const ChannelStatue& rhs) { high |= rhs.high; low |= rhs.low; return *this; }
    ChannelStatue& operator &= (const ChannelStatue& rhs) { high &= rhs.high; low &= rhs.low; return *this; }
    ChannelStatue& operator ^= (const ChannelStatue& rhs) { high ^= rhs.high; low ^= rhs.low; return *this; }

    ChannelStatue& operator += (const ChannelStatue& rhs) { const uint64_t old = low; low += rhs.low;  high += rhs.high + (low < old); return *this; }
    ChannelStatue& operator -= (const ChannelStatue& rhs) { return *this += -rhs; }

    friend ChannelStatue operator + (const ChannelStatue& l, const ChannelStatue& r)   { return ChannelStatue(l) += r; }
    friend ChannelStatue operator + (const ChannelStatue& l, const uint64_t& r)   { return ChannelStatue(l) += ChannelStatue(r); }
    friend ChannelStatue operator + (const ChannelStatue& l, const uint32_t& r)   { return ChannelStatue(l) += ChannelStatue(r); }
    friend ChannelStatue operator + (const ChannelStatue& l, const int32_t& r)   { return ChannelStatue(l) += ChannelStatue(r); }
    friend ChannelStatue operator + (const uint64_t& l, const ChannelStatue& r)   { return ChannelStatue(l) += r; }
    friend ChannelStatue operator - (const ChannelStatue& l, const ChannelStatue& r)   { return ChannelStatue(l) -= r; }
    friend ChannelStatue operator | (const ChannelStatue& l, const ChannelStatue& r)   { return ChannelStatue(l) = (r); }
    friend ChannelStatue operator & (const ChannelStatue& l, const ChannelStatue& r)   { return ChannelStatue(l) &= r; }
    friend ChannelStatue operator & (const ChannelStatue& l, const uint64_t& r)   { return ChannelStatue(l) &= ChannelStatue(r); }
    friend ChannelStatue operator ^ (const ChannelStatue& l, const ChannelStatue& r)   { return ChannelStatue(l) ^= r; }
    friend bool    operator >  (const ChannelStatue& l, const ChannelStatue& r)  { return r < l; }
    friend bool    operator >  (const ChannelStatue& l, const int64_t& r)  { return ChannelStatue(r) < l; }
    friend bool    operator >  (const int64_t& l, const ChannelStatue& r)  { return r < ChannelStatue(l); }

    friend bool    operator >=  (const ChannelStatue& l, const ChannelStatue& r) { return l == r || l > r; }
    friend bool    operator >=  (const ChannelStatue& l, const int64_t& r) { return l >= ChannelStatue(r); }
    friend bool    operator >=  (const int64_t& l, const ChannelStatue& r) { return ChannelStatue(l) >= r; }
    friend bool    operator <=  (const ChannelStatue& l, const ChannelStatue& r) { return l == r || l < r; }
    friend bool    operator <=  (const ChannelStatue& l, const int64_t& r) { return l <= ChannelStatue(r); }
    friend bool    operator <=  (const int64_t& l, const ChannelStatue& r) { return ChannelStatue(l) <= r; }

    inline bool is_nonzero() { return high != 0 || low != 0; }
    inline int get_bit(int pos);
    inline bool set_bit_one(int pos);
    inline bool set_bit_zero(int pos);

    uint64_t high;
    uint64_t low;
};

struct Edge {
    int id;                                  // 边ID
    int distance;                            // 边长度
    int jump;                                // 边跳数
    int nid0;                                // 边的一个端点
    int nid1;                                // 边的另一个端点
    int base_eid = -1;                       // 新增边由 base_eid 拷贝而来
    int count_channel_used = 0;
    double priority = 0.0;
    vector<int> channel2operation;           // channel2operation[i] = j, j = -1 表示 通道i 未使用, 否则表示 通道i 被 业务j 使用
    vector<int> channel2path;                // channel2operation[i] = j, j = -1 表示 通道i 未使用, 否则表示 通道i 被 路径j 使用
    ChannelStatue channel_statue = ChannelStatue(0, 0);                  // 状态压缩, 使用 channel_statue 表示
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
    unordered_map<int, unordered_set<int>> cid2pid; // 本业务使用的通道编号及使用其的路径编号
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
int num_of_paths;                                                       // 最大业务数量
int num_of_channels;                                                    // 最大通道号
int distance_attenuation_max;                                           // 最大衰减距离
int jump_attenuation_max;                                               // 最大跳数
ChannelStatue channel_statue_empty = 0;                                 // 通道完全未使用的边的通道状态信息
std::chrono::steady_clock::time_point time_begin;                       // 记录程序开始的时刻
vector<Node> nodes;                                                     // num_of_nodes 是节点个数, ID 从[0, num_of_nodes-1]
vector<Edge> edges;                                                     // edges.size() 是边个数, ID 从[0, edges.size()-1]
vector<Operation> operations;                                           // num_of_operations 是业务个数, ID 从[0, num_of_operations-1]
vector<Path> paths;                                                     // num_of_paths 是路径个数, ID 从[0, num_of_paths-1]
unordered_map<int, unordered_map<int, int>> dists_min;                  // 记录两个节点之间直接连接的最短距离: 新生成的边的长度
unordered_map<int, unordered_map<int, vector<int>>> paths_of_nodes;     // 以节点序列记录两个节点之间的路径
vector<int> count_of_channel_used;                                      // 记录每个通道的使用次数
unordered_map<int, unordered_map<int, int>> jumps_min;                  // 记录两个节点之间直接连接的最短跳数的边: 新生成的边由此拷贝而来

/****************************************************以下是函数声明******************************************************/
vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id, unordered_set<int>& baned_edges);   // 返回源点与目标点之间的节点路径(channel_id 为 -1 时不考虑通道占用情况)
pair<int, vector<int>> find_cheapest_path_nodes(int nid_source,
                                                int nid_target,
                                                int cid,
                                                int cost_cheapest,
                                                unordered_set<int>& baned_edges);

pair<ChannelStatue, vector<int>> bfs_find_path_nodes_state_compression(int nid_source, int nid_target, unordered_set<int>& baned_edges);
ChannelStatue bfs_find_channel_id_unoccupied(int nid_source, int nid_target, unordered_set<int>& baned_edges);// 在不加边的情况下, 找出满足要求的通道状态

int add_edge(int nid0, int nid1);                        // 新增两个节点之间的边
int get_eid_channel_id_ok(int nid0, int nid1, int channel_id, unordered_set<int>& baned_edges);
void sort_data();
void place_amplifier_in_path();                          // 在路径中放置放大器
void place_amplifier(Path& ph);                          // 为每个路径放置放大器, 此函数可贪心得到最优解, 无需修改!
void optimization_select_shortest_edge() {};             // 在保证有解的情况下, 选择距离更短的边(只考虑两个节点之间的换边, 因赛题修改, 此处待重新实现)
void optimization_transfer_operation() {};               // 将使用率低的边的业务转移至其它边从而移除使用率低的边(因赛题修改, 此处待重新实现)
void optimization_path_rebuild();
bool pave();                                             // 铺设线路

void print_map();
void print_channel_used_total();
void print_info(stringstream& scout);

void path_edges_occupy(Path& ph);                        // 更新边的占用信息, 释放边
void path_edges_release(Path& ph);                       // 更新边的占用信息, 占用边

static inline bool timeout() {                           // 记录从程序起始到现在过去了多少时间
    using namespace std::chrono;
    std::chrono::steady_clock::time_point time_current = std::chrono::steady_clock::now();
    duration<double> delta = duration_cast<duration<double>>(time_current - time_begin);
    return delta.count() > process_time_max;
}

/****************************************************以下是主处理函数******************************************************/
void init_data() {
    // 输入
    cin >> num_of_nodes
        >> num_of_edges_initial
        >> num_of_operations
        >> num_of_paths
        >> num_of_channels
        >> distance_attenuation_max
        >> jump_attenuation_max;

    // 初始化基础变量
    nodes.resize(num_of_nodes);
    edges.resize(num_of_edges_initial);
    operations.resize(num_of_operations);
    paths.resize(num_of_paths);
    count_of_channel_used = vector<int>(num_of_channels, 0);

    // 把可用的通道编号位置为 1
    for (int i = 0; i < num_of_channels; ++i) {
        if (i < 64) {
            channel_statue_empty.low |= (uint64_t(1) << i);
        } else {
            channel_statue_empty.high |= (uint64_t(1) << (i - 64));
        }
    }

    // 初始化节点，节点数目不会变化
    for (int i = 0; i < num_of_nodes; ++i) {
        nodes[i].id = i;
    }

    // 初始化初始边，新增边必须节点间存在初始边
    for (int i = 0; i < num_of_edges_initial; ++i) {    // 读入 num_of_edges 行边
        cin >> edges[i].nid0 >> edges[i].nid1 >> edges[i].distance >> edges[i].jump;
        edges[i].id = i;
        edges[i].channel2operation = vector<int>(num_of_channels, -1);
        edges[i].channel2path = vector<int>(num_of_channels, -1);
        edges[i].channel_statue = channel_statue_empty;
        edges[i].priority = (1.0 * edges[i].distance / distance_attenuation_max) + (1.0 * edges[i].jump / jump_attenuation_max);
        nodes[edges[i].nid0].eids.push_back(edges[i].id);
        nodes[edges[i].nid1].eids.push_back(edges[i].id);
        // 记录节点边的最短距离
//        if (dists_min[edges[i].nid0].count(edges[i].nid1)) {
//            dists_min[edges[i].nid0][edges[i].nid1]
//                    = dists_min[edges[i].nid1][edges[i].nid0]
//                    = min(edges[i].distance, dists_min[edges[i].nid0][edges[i].nid1]);
//        } else {
//            dists_min[edges[i].nid0][edges[i].nid1]
//                    = dists_min[edges[i].nid1][edges[i].nid0]
//                    = edges[i].distance;
//        }
        // 记录节点之间最小跳数的那条边
        if (jumps_min[edges[i].nid0].count(edges[i].nid1)) {
            if (edges[i].priority < edges[jumps_min[edges[i].nid0][edges[i].nid1]].priority) {
                jumps_min[edges[i].nid0][edges[i].nid1] = jumps_min[edges[i].nid1][edges[i].nid0] = i;
            }
        } else {
            jumps_min[edges[i].nid0][edges[i].nid1] = jumps_min[edges[i].nid1][edges[i].nid0] = i;
        }

    }

    // 初始化业务
    int path_id = 0;
    for (int i = 0; i < num_of_operations; ++i) {   // 读入 num_of_operations 行业务
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
#ifdef TEST_HDU01
    print_map();
#endif
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
    auto count_channel_used_backup = count_of_channel_used;
    auto paths_best = paths;
    auto edges_best = edges;
    auto count_channel_used_best = count_of_channel_used;
    int edges_new_smallest = INF;

    while (true) {
        // 备份初始数据
        nodes = nodes_backup;
        edges = edges_backup;
        operations = operations_backup;
        paths = paths_backup;
        count_of_channel_used = count_channel_used_backup;

        sort_data();// 对业务进行排序

        if (!pave()) {// 铺设线路
            return;// 出现bug
        }

        // 多次计算获取随机到的最优解
        int edges_new = edges.size() - num_of_edges_initial;
        if (edges_new < edges_new_smallest) {
            edges_new_smallest = edges_new;
            paths_best = paths;
            edges_best = edges;
            count_channel_used_best = count_of_channel_used;
        }
        // 程序达到设定运行时间，退出
        if (timeout()) {
            break;
        }

    }

    // 获取最优解
    paths = paths_best;
    edges = edges_best;

    // 得到最优解后，在路径中央对放大器进行放置
    place_amplifier_in_path();

#ifdef TEST_HDU01
    print_channel_used_total();     // 打印通道使用
#endif

}

// 输出数据处理完的结果
void output_data() {
    // 输出边的信息，新增边个数等于 最终边 - 初始边
    cout << edges.size() - num_of_edges_initial << endl;
    for (int i = num_of_edges_initial; i < edges.size(); ++i) {
//        cout << edges[i].nid0 << " " << edges[i].nid1 << endl;
        cout << edges[i].base_eid << endl;
    }

    // 输出业务通道信息
    sort(paths.begin(), paths.end(), [&](Path& ph0, Path& ph1) {
        return ph0.id < ph1.id;
    });

    // 输出放大器的信息
    for (auto& path : paths) {
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
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        scout << "Node ID: " << node.id << ", ";
        for (auto& eid : node.eids) {
            auto& edge = edges[eid];
            scout << "(eid:" << edge.id << ",dist=" << edge.distance << ",nid_neighbor=" << edge.nid0 + edge.nid1 - node.id << ") ";
        }
        scout << "\n";
    }
    scout << "建图信息输出完毕\n";
    // 输出业务信息
    scout << "\n业务信息如下: \n";
    for (int pid = 0; pid < num_of_paths; ++pid) {
        auto& path = paths[pid];
        scout << "path.id = " << path.id << ", nid0 = " << path.nid0 << ", nid1 = " << path.nid1 << "\n";
    }
    scout << "业务信息输出完毕\n";

//    scout << "边状态压缩数据(重载<<=): ";
//    for (int cid = 128 - 1; cid >= 0; --cid) {
//        ChannelStatue tmp = 1;
//        tmp <<= uint(cid);
//        ChannelStatue statue = channel_statue_empty & tmp;
//        if (statue != ChannelStatue(0)) {
//            scout << "1";
//        } else {
//            scout << "0";
//        }
//    }
//    scout << "\n";
    scout << "边状态压缩数据: ";
    for (int cid = 128 - 1; cid >= 0; --cid) {
        scout << channel_statue_empty.get_bit(cid);
        if (cid != 0 && cid % 4 == 0) scout << "-";
    }
    scout << "\n";
    // 执行输出
    print_info(scout);
}

void print_channel_used_total() {
    scout << "每个通道的使用次数: ";
    for (int i = 0; i < num_of_channels; ++i) {
        scout << count_of_channel_used[i];
        if (i + 1 != num_of_channels) cout << ",";
        else scout << ".\n";
    }
    // 执行输出
    print_info(scout);
}

void print_info(stringstream& scout) {
#ifdef TEST_HDU01
    cout << scout.str() << '\n';
#endif
    // 清空字符串流
    scout.str("");
}

/****************************************************以下是排序相关函数的定义******************************************************/
// TODO: 排序如何优化
void sort_data() {

    // 对每个节点的邻接边进行排序, 把距离短的放在前面
    // TODO: 关键边? 或许应该把不怎么被使用的边放在前面, 把被很多业务使用的边放在后面
    for (int nid = 0; nid < num_of_nodes; ++nid) {
        auto& node = nodes[nid];
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
//            return edges[eid0].distance < edges[eid1].distance;
//            return edges[eid0].jump < edges[eid1].distance;
            return edges[eid0].priority < edges[eid1].priority;
        });
    }

    // 对业务进行排序

    //  随机排序
    random_device rd;
    mt19937 g(rd());

//    shuffle(operations.begin(), operations.end(), g);

    // TODO: 尝试两种shuffle, 一是2硬约束业务与普通业务 一起打乱; 二是找到分界线, 将2硬约束业务与普通业务 分别打乱
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

    return;

}

/****************************************************以下是主要搜索函数的定义******************************************************/
// TODO: 搜索如何优化

vector<int> bfs_find_path_nodes(int nid_source, int nid_target, int channel_id, unordered_set<int>& baned_edges) {
    // TODO: 使用双向搜索来加速 https://oi-wiki.org/search/bidirectional/
    // 当能不添加新边, 能通过当前通道连接这两个节点时, 返回经过的节点序列(从源点到目标点)
    // 否则返回一个空序列
    // channel_id == -1 表示不考虑通道编号
    unordered_map<int, int> prev;
    unordered_set<int> visited_edges;// 防止重复遍历同一条边
    queue<int> que;
    que.push(nid_source);
    prev[nid_source] = -1;
    while (!que.empty() && !prev.count(nid_target)) {
        int nid = que.front();
        que.pop();
        for (auto& eid : nodes[nid].eids) {
            auto& edge = edges[eid];
            if (baned_edges.count(eid) || visited_edges.count(eid)) continue;// 已被本业务其它路径使用
            visited_edges.insert(eid);
            if (channel_id != -1 && edge.channel2operation[channel_id] != -1) continue;// 考虑通道使用情况
            int nid_next = edge.nid0 + edge.nid1 - nid;
            if (prev.count(nid_next)) {
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

pair<ChannelStatue, vector<int>> bfs_find_path_nodes_state_compression(int nid_source, int nid_target, unordered_set<int>& baned_edges) {

    unordered_map<int, int> prev;
    queue<pair<int, ChannelStatue>> que;// get<0> 是节点ID, get<1> 是当前能通过的通道编号
    que.push(make_pair(nid_source, channel_statue_empty));// 初始设置为所有通道均可用

    ChannelStatue result_state(0);// 路径可用的通道编号
    vector<int> path_nodes;// 路径对应的节点序列

    while (!que.empty()) {
        auto [nid, state] = que.front();
        que.pop();
        if (nid == nid_target) {
            result_state = state;
            break;
        }
        for (auto& eid : nodes[nid].eids) {
            if (baned_edges.count(eid)) continue;
            ChannelStatue state_next = state & edges[eid].channel_statue;
            if (state_next.low != 0 || state_next.high != 0) {
                // TODO: 如何记忆化? 如何记录前驱节点?

            }
        }

    }

    if (!prev.count(nid_target)) {
        return make_pair(ChannelStatue(0), path_nodes);
    }

    int nid = nid_target;
    while (nid != -1) {
        path_nodes.push_back(nid);
        nid = prev[nid];
    }

    reverse(path_nodes.begin(), path_nodes.end());

    return make_pair(result_state, path_nodes);
}

ChannelStatue bfs_find_channel_id_unoccupied(int nid_source, int nid_target, unordered_set<int>& baned_edges) {

    // 在不加边的情况下, 找出满足要求的状态通道
    // 返回 0: 表示无可用的通道编号
    // 否则, 返回可用的通道编号

    queue<pair<int, ChannelStatue>> que;// get<0> 是节点ID, get<1> 是当前能通过的通道编号
    que.push(make_pair(nid_source, channel_statue_empty));// 初始设置为所有通道均可用
    unordered_set<int> visited_edges;
    vector<ChannelStatue> state_of_nodes(num_of_nodes, channel_statue_empty);// 记录节点的通道进入状态, 保证每个通道只会进入节点一次
    while (!que.empty()) {
        auto [nid, state] = que.front();
        que.pop();
        if (nid == nid_target) {
            return state;
        }
        for (auto& eid : nodes[nid].eids) {
            if (baned_edges.count(eid) || visited_edges.count(eid)) continue;
            visited_edges.insert(eid);
            int nid_next = edges[eid].nid0 + edges[eid].nid1 - nid;
            ChannelStatue state_next = state & edges[eid].channel_statue & state_of_nodes[nid_next];
            if (state_next.is_nonzero()) {
                // 更新节点 nid_next 的通道入度状态信息
                // 因为 state_of_nodes[nid_next] 在 state_next 上的1位置上均为1, 所以此举相当于将这些位置为0
                state_of_nodes[nid_next] -= state_next;
                que.push(make_pair(nid_next, state_next));
            }
        }

    }

    return 0;

}

// TODO 初步通过随机数进行剪枝，是否需要进一步优化？如何优化？
pair<int, vector<int>> find_cheapest_path_nodes(int nid_source,
                                                int nid_target,
                                                int cid,
                                                int cost_cheapest,
                                                unordered_set<int>& baned_edges) {

    // TODO: 在加边数量相同的情况下, 选择经过边的数量最小的路径
    // get<0> 表示代价
    unordered_map<int, int> prev_nid, prev_cost;
    unordered_set<int> visited_edges;// 防止重复遍历同一条边
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
            if (visited_edges.count(eid)) continue;
            visited_edges.insert(eid);
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

/****************************************************以下是路径铺设相关函数的定义******************************************************/
// TODO: 主要功能函数，如何进一步拆分
int add_edge(int nid0, int nid1) {
    Edge edge;
    edge.id = edges.size();
    edge.nid0 = nid0;
    edge.nid1 = nid1;
    edge.base_eid = jumps_min[nid0][nid1];
    edge.channel2operation = vector<int>(num_of_channels, -1);
    edge.channel2path = vector<int>(num_of_channels, -1);
    edge.channel_statue = channel_statue_empty;
    edge.distance = edges[edge.base_eid].distance;// 新边的长度与跳数继承拷贝边的数据
    edge.jump = edges[edge.base_eid].jump;
    edge.priority = edges[edge.base_eid].priority;
    nodes[edge.nid0].eids.push_back(edge.id);
    nodes[edge.nid1].eids.push_back(edge.id);

    edges.emplace_back(edge);

    return edge.id;
}

int get_eid_channel_id_ok(int nid0, int nid1, int channel_id, unordered_set<int>& baned_edges) {
    // 获取两个节点之间的边, 该边的 channel_id 未被使用
    // 返回 -1 表示无现成的可用边在此通道空闲
    int res = -1;
    for (auto& eid : nodes[nid0].eids) {
        auto& edge = edges[eid];
        if (baned_edges.count(eid)) continue;
        if (edge.channel2operation[channel_id] != -1) continue;
        int nid_neighbor = edge.nid0 + edge.nid1 - nid0;
        if (nid_neighbor != nid1) continue;
        res = eid;
        break;
    }
    return res;
}

void path_edges_occupy(Path& ph) {
    // 请保证 path_nodes, path_edges 及 path_channel_id 是可行的
    auto& op = operations[ph.operation_id];
    const int& path_channel_id = ph.channel_id;
    const auto& path_edges = ph.edges;
    const auto& path_nodes = ph.nodes;
    count_of_channel_used[path_channel_id] += path_edges.size();

#ifdef TEST_HDU01
    if (path_edges.empty()) {
        cerr << "Error: path_edges empty!\n";
        return;
    }
#endif

    for (auto& eid : path_edges) {
        auto& edge = edges[eid];
#ifdef TEST_HDU01
        if (edge.channel2operation[path_channel_id] != -1) {
            cerr << "Error: Reusing the same channel.(operation)\n";
        }
        if (edge.channel2path[path_channel_id] != -1) {
            cerr << "Error: Reusing the same channel.(path)\n";
        }
        if (edge.channel_statue.get_bit(path_channel_id) != 1) {
            cerr << "Error: Reusing the same channel.(state)\n";
        }
#endif
        edge.channel2operation[path_channel_id] = op.id;
        edge.channel2path[path_channel_id] = ph.id;
        edge.channel_statue.set_bit_zero(path_channel_id);
        ++edge.count_channel_used;
    }

#ifdef TEST_HDU01
    scout << "path.id = " << ph.id << ", channel_id = " << ph.channel_id << ": ";
    for (int i = 0; i < path_nodes.size(); ++i) {
        scout << path_nodes[i];
        if (i + 1 != path_nodes.size()) scout << "->";
    }
    scout << "\n";
//    scout << "op.channel_id: ";
//    for (auto& [cid, pids] : op.cid2pid) {
//        scout << cid << ", ";
//    }
//    scout << "\n";
    print_info(scout);
#endif

}

void path_edges_release(Path& ph) {
    auto& op = operations[ph.operation_id];
    const int& path_channel_id = ph.channel_id;
    const auto& path_edges = ph.edges;
    const auto& path_nodes = ph.nodes;
    count_of_channel_used[path_channel_id] -= path_edges.size();

    for (auto& eid : path_edges) {
        auto& edge = edges[eid];
#ifdef TEST_HDU01
        if (edge.channel2operation[path_channel_id] == -1) {
            cerr << "Error: Release an empty edge.(operation)\n";
        }
        if (edge.channel2path[path_channel_id] == -1) {
            cerr << "Error: Release an empty edge.(path)\n";
        }
        if (edge.channel_statue.get_bit(path_channel_id) != 0) {
            cerr << "Error: Release an empty edge.(state)\n";
        }
#endif
        edge.channel2operation[path_channel_id] = -1;
        edge.channel2path[path_channel_id] = -1;
        edge.channel_statue.set_bit_one(path_channel_id);
        --edge.count_channel_used;
    }

}

bool pave() {

    vector<int> pids;// 属于非2路径业务的路径

    for (auto& op : operations) {
        if (op.paths.size() != 2) {
            for (auto& pid : op.paths) {
                pids.push_back(pid);
            }
            continue;
        }
        for (int i = 0; i < op.paths.size(); ++i) {
            auto& path = paths[op.paths[i]];
            int nid_source = path.nid0, nid_target = path.nid1;
            int path_channel_id = -1;
            vector<int> path_nodes;
            vector<int> path_edges;
            if (i == 1 && op.paths.size() == 2) {   // 硬约束: 这条路径必须得与前一条路径使用的通道相同
                if (op.cid2pid.empty()) {           // 严重错误!!!
                    cerr << "Error: op.channel_id empty!\n";
                    return false;
                }
                path_channel_id = (*op.cid2pid.begin()).first;
                vector<int> path_nodes_tmp = bfs_find_path_nodes(nid_source, nid_target, path_channel_id, op.edges_used);
                if (path_nodes_tmp.empty()) {
                    // 需要加边来满足要求
#ifdef TEST_HDU01
                    scout << "未找到相同通道的可用线路\n";
                    print_info(scout);
#endif
                    auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, path_channel_id, INF, op.edges_used);
                    path_nodes = get<1>(pair_cost_nodes);
                } else {
#ifdef TEST_HDU01
                    scout << "找到了相同通道的可用线路\n";
                    print_info(scout);
#endif
                    path_nodes = path_nodes_tmp;
                }

            }
            else {

//                ChannelStatue channel_id_unoccupied = bfs_find_channel_id_unoccupied(nid_source, nid_target, op.edges_used);
//                if (channel_id_unoccupied.is_nonzero()) {
//#ifdef TEST_HDU01
//                    scout << "找到了可用线路(状态压缩)\n";
//                    print_info(scout);
//#endif
//                    int randomStart = rand() % num_of_channels;
//                    for (int i = 0; i < num_of_channels; ++i) {
//                        int cid = (randomStart + i) % num_of_channels;
//                        if (channel_id_unoccupied.get_bit(cid) == 1) {
//                            path_channel_id = cid;
//                            path_nodes = bfs_find_path_nodes(nid_source, nid_target, cid, op.edges_used);
//                            break;
//                        }
//                    }
//
//                }
//                else {
//#ifdef TEST_HDU01
//                    scout << "未找到可用线路(状态压缩)\n";
//                    print_info(scout);
//#endif
//                    unordered_set<int> baned_edges_emtpy;
//                    path_nodes = bfs_find_path_nodes(nid_source, nid_target, -1, baned_edges_emtpy);
//                    // 得到一个基础解, 用于剪枝
//                    path_channel_id = rand() % num_of_channels;
//                    int cost_cheapest = 0;
//                    for (int i = 1; i < path_nodes.size(); ++i) {
//                        // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
//                        int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id, op.edges_used);
//                        if (eid_old == -1) {
//                            ++cost_cheapest;
//                        }
//                    }
//                    // 随机一个通道编号, 求这个通道编号的最优路线
//                    // TODO: 基于整个图边通道的使用情况, 选择通道使用最少的那个通道编号
//                    int cid = rand() % num_of_channels;
//                    for (int i = 0; i < num_of_channels; ++i) {
//                        if (count_of_channel_used[i] < count_of_channel_used[cid]) {
//                            cid = i;
//                        }
//                    }
////#ifdef TEST_HDU01
////                    for (int i = 0; i < num_of_channels; ++i) {
////                        cout << count_of_channel_used[i];
////                        if (i + 1 != num_of_channels) cout << ",";
////                        else cout << ".";
////                    }
////#endif
//                    auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, cost_cheapest, op.edges_used);
//                    auto cost_tmp = get<0>(pair_cost_nodes);
//                    auto path_nodes_tmp = get<1>(pair_cost_nodes);
//                    // 将随机求得的最优路线与基础解做比较, 选择更优的
//                    if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
//                        cost_cheapest = cost_tmp;
//                        path_nodes = path_nodes_tmp;
//                        path_channel_id = cid;
//                    }
//
//                }

                // 以下是复赛使用代码
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
                // 随机一个通道编号, 求这个通道编号的最优路线
                int cid = rand() % num_of_channels;
                auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, cost_cheapest, op.edges_used);
                auto cost_tmp = get<0>(pair_cost_nodes);
                auto path_nodes_tmp = get<1>(pair_cost_nodes);
                // 将随机求得的最优路线与基础解做比较, 选择更优的
                if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
                    cost_cheapest = cost_tmp;
                    path_nodes = path_nodes_tmp;
                    path_channel_id = cid;
                }
                
            }


            // 用于快速提交解决方案, 得到一个初始解
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

            path.channel_id = path_channel_id;
            path.nodes = path_nodes;
            path.edges = path_edges;

            // 更新占用信息
            path_edges_occupy(path);
            op.edges_used.insert(path_edges.begin(), path_edges.end());
            op.cid2pid[path_channel_id].insert(path.id);

        }

    }


    random_device rd;
    mt19937 g(rd());

    shuffle(pids.begin(), pids.end(), g);

    for (auto& pid : pids) {

        auto& path = paths[pid];
        auto& op = operations[path.operation_id];
        int nid_source = path.nid0, nid_target = path.nid1;
        int path_channel_id = -1;
        vector<int> path_nodes;
        vector<int> path_edges;
        // 以下是复赛使用代码
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
        // 随机一个通道编号, 求这个通道编号的最优路线
        int cid = rand() % num_of_channels;
        auto pair_cost_nodes = find_cheapest_path_nodes(nid_source, nid_target, cid, cost_cheapest, op.edges_used);
        auto cost_tmp = get<0>(pair_cost_nodes);
        auto path_nodes_tmp = get<1>(pair_cost_nodes);
        // 将随机求得的最优路线与基础解做比较, 选择更优的
        if (!path_nodes_tmp.empty() && cost_tmp < cost_cheapest) {
            cost_cheapest = cost_tmp;
            path_nodes = path_nodes_tmp;
            path_channel_id = cid;
        }

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

        path.channel_id = path_channel_id;
        path.nodes = path_nodes;
        path.edges = path_edges;

        // 更新占用信息
        path_edges_occupy(path);
        op.edges_used.insert(path_edges.begin(), path_edges.end());
        op.cid2pid[path_channel_id].insert(path.id);

    }


    return true;

}

/****************************************************以下是放大器主要函数的定义******************************************************/
// TODO: 如何优化放大器防止成本

void place_amplifier(Path& path) {

    auto& path_nodes = path.nodes;
    auto& path_edges = path.edges;

    // 放置放大器
    vector<int> path_amplifiers;
    int dist_sum = edges[path_edges.front()].distance;
    int jump_sum = edges[path_edges.front()].jump;
    for (int i = 1; i < path_edges.size(); ++i) {
        // dist_sum/jump_sum  表示当前到达这个节点的线路长度/跳数次数
        auto& edge = edges[path_edges[i]];
        if (dist_sum + edge.distance >= distance_attenuation_max || jump_sum + edge.jump >= jump_attenuation_max) {
            path_amplifiers.push_back(path_nodes[i]);// 其实是这条边的左端点, 注意: 边个数与节点个数差 1
            dist_sum = edge.distance;
            jump_sum = edge.jump;
        } else {
            dist_sum += edge.distance;
            jump_sum += edge.jump;
        }
    }
    path.amplifiers = path_amplifiers;
}

void optimization_path_rebuild() {
    // 在不更换路径通道的情况下重新规划路线
    // TODO: 在更换路径通道的情况下, 如何不超时?
    // 把所有节点的边按长度升序排
    for (auto& node : nodes) {
        sort(node.eids.begin(), node.eids.end(), [&](int& eid0, int& eid1) {
            return edges[eid0].distance < edges[eid1].distance;
        });
    }

    for (auto& path : paths) {
        place_amplifier(path);
    }

    sort(paths.begin(), paths.end(), [&] (Path& ph0, Path& ph1) {
        return ph0.amplifiers.size() > ph1.amplifiers.size();
    });

    for (auto& path : paths) {
        auto& op = operations[path.operation_id];
        // 先释放已占用的边
        path_edges_release(path);
        for (auto& eid : path.edges) {
            op.edges_used.erase(eid);
        }
        op.cid2pid[path.channel_id].erase(path.id);
        if (op.cid2pid[path.channel_id].empty()) {
            op.cid2pid.erase(path.channel_id);
        }
        // 重新寻找一个可行解(必不需要新增边)
        int path_channel_id = path.channel_id;
        vector<int> path_nodes, path_edges;
        path_nodes = bfs_find_path_nodes(path.nid0, path.nid1, path_channel_id, op.edges_used);
        // 根据节点序列选择边
        for (int i = 1; i < path_nodes.size(); ++i) {
            // 检查是否存在 通道 path_channel_id 为空闲的边, 尽量少建新边
            int eid_old = get_eid_channel_id_ok(path_nodes[i - 1], path_nodes[i], path_channel_id, op.edges_used);
            if (eid_old != -1) {
                path_edges.emplace_back(eid_old);
                continue;
            }
            // 应是出现bug, 否则不可能到达此处
            // 需要加边, 直接恢复
            path_channel_id = path.channel_id;
            path_nodes = path.nodes;
            path_edges = path.edges;
            break;
        }

        path.channel_id = path_channel_id;
        path.nodes = path_nodes;
        path.edges = path_edges;
        // 占用边
        path_edges_occupy(path);
        op.edges_used.insert(path_edges.begin(), path_edges.end());// 记录这条业务使用的边
        op.cid2pid[path_channel_id].insert(path.id);

    }

}

void place_amplifier_in_path() {
//    // optimization_transfer_operation();// 将通道使用率低的边上的业务转移到其它边上从而移除掉该边, 但目前未奏效
//    optimization_path_rebuild();// 在保证可行性的前提下重新铺设线路, 这个能降低成本, 但偶尔会超时? 74847907 -> 74785404

    // 放大器的优化
    // 初始化: 放置放大器
    for (auto& op : operations) {
        for (auto& pid : op.paths) {
            place_amplifier(paths[pid]);
        }
    }

    // 后处理: 在只更换边的情况下, 检查是否存在路径长度更短的连接方案
    // optimization_select_shortest_edge();
}

inline int ChannelStatue::get_bit(int pos) {
    int res;
    if (pos < 64) {
        res = ((this->low >> pos) & 1);
    } else {
        res = ((this->high >> (pos - 64)) & 1);
    }
    return res;
}

inline bool ChannelStatue::set_bit_one(int pos) {
    if (this->get_bit(pos) == 1) return false;// 防止出现错误
    if (pos < 64) {
        this->low += (uint64_t(1) << pos);
    } else {
        this->high += (uint64_t(1) << (pos - 64));
    }
    return true;
}

inline bool ChannelStatue::set_bit_zero(int pos) {
    if (this->get_bit(pos) == 0) return false;// 防止出现错误
    if (pos < 64) {
        this->low -= (uint64_t(1) << pos);
    } else {
        this->high -= (uint64_t(1) << (pos - 64));
    }
    return true;
}