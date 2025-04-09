#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <random>
#include <cassert>
#include <chrono>
#include <numeric>
#include <unordered_map>

using namespace std;

// 乱数生成器
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// 都市の座標範囲
struct CityRange {
    int lx, rx, ly, ry;
};

// 都市間の距離情報
struct Edge {
    int u, v;
    double dist;
    bool is_real; // 実際の距離かどうか
    
    bool operator<(const Edge& other) const {
        if (is_real != other.is_real) return is_real > other.is_real; // 実際の距離を優先
        return dist < other.dist;
    }
};

// Union-Find
class UnionFind {
private:
    vector<int> parent, rank;
    
public:
    UnionFind(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }
    
    void unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) return;
        
        if (rank[x] < rank[y]) {
            parent[x] = y;
        } else {
            parent[y] = x;
            if (rank[x] == rank[y]) rank[x]++;
        }
    }
    
    bool same(int x, int y) {
        return find(x) == find(y);
    }
};

// 占いの結果から距離行列を更新する
void updateDistanceMatrix(vector<vector<double>>& dist_matrix, vector<vector<bool>>& is_real_dist, 
                          const vector<int>& cities, const vector<pair<int, int>>& mst_edges) {
    // 都市間の接続情報を構築
    int n = cities.size();
    vector<vector<pair<int, double>>> graph(n);
    
    // 辺の実際の距離を計算
    unordered_map<int, int> city_to_idx;
    for (int i = 0; i < n; i++) {
        city_to_idx[cities[i]] = i;
    }
    
    // 最小全域木の辺を追加
    for (const auto& edge : mst_edges) {
        int u = edge.first;
        int v = edge.second;
        
        int idx_u = city_to_idx[u];
        int idx_v = city_to_idx[v];
        
        // 都市間の距離を推定（中心座標間のユークリッド距離）
        double dx = dist_matrix[u][v];
        
        graph[idx_u].push_back({idx_v, dx});
        graph[idx_v].push_back({idx_u, dx});
        
        // 実際の距離として記録
        dist_matrix[u][v] = dist_matrix[v][u] = dx;
        is_real_dist[u][v] = is_real_dist[v][u] = true;
    }
    
    // 最短経路を計算して距離行列を更新
    for (int i = 0; i < n; i++) {
        int start_city = cities[i];
        
        // ダイクストラ法で最短経路を計算
        vector<double> distance(n, 1e9);
        vector<bool> visited(n, false);
        
        distance[i] = 0;
        
        for (int j = 0; j < n; j++) {
            int u = -1;
            for (int k = 0; k < n; k++) {
                if (!visited[k] && (u == -1 || distance[k] < distance[u])) {
                    u = k;
                }
            }
            
            if (distance[u] == 1e9) break;
            
            visited[u] = true;
            
            for (auto& [v, w] : graph[u]) {
                if (distance[u] + w < distance[v]) {
                    distance[v] = distance[u] + w;
                }
            }
        }
        
        // 距離行列を更新
        for (int j = 0; j < n; j++) {
            if (i != j) {
                int end_city = cities[j];
                
                // 既に実際の距離が分かっている場合はスキップ
                if (is_real_dist[start_city][end_city]) continue;
                
                // 最短経路の距離を使用
                if (distance[j] < 1e9) {
                    // 既存の推定値より小さい場合のみ更新
                    if (distance[j] < dist_matrix[start_city][end_city]) {
                        dist_matrix[start_city][end_city] = dist_matrix[end_city][start_city] = distance[j];
                    }
                }
            }
        }
    }
    
    // 三角不等式を利用して距離の上限と下限を更新
    for (int k = 0; k < cities.size(); k++) {
        int city_k = cities[k];
        for (int i = 0; i < cities.size(); i++) {
            if (i == k) continue;
            int city_i = cities[i];
            for (int j = 0; j < cities.size(); j++) {
                if (j == k || j == i) continue;
                int city_j = cities[j];
                
                // 三角不等式: dist(i,j) <= dist(i,k) + dist(k,j)
                double new_dist = dist_matrix[city_i][city_k] + dist_matrix[city_k][city_j];
                if (new_dist < dist_matrix[city_i][city_j]) {
                    dist_matrix[city_i][city_j] = dist_matrix[city_j][city_i] = new_dist;
                }
            }
        }
    }
}

// 都市の中心座標を推定
pair<double, double> estimateCityCenter(const CityRange& range) {
    return {(range.lx + range.rx) / 2.0, (range.ly + range.ry) / 2.0};
}

// 都市間の推定距離を計算
double estimateDistance(const pair<double, double>& city1, const pair<double, double>& city2) {
    double dx = city1.first - city2.first;
    double dy = city1.second - city2.second;
    return sqrt(dx * dx + dy * dy);
}

// クエリを実行して最小全域木を取得
vector<pair<int, int>> query(const vector<int>& cities) {
    int l = cities.size();
    
    cout << "? " << l;
    for (int city : cities) {
        cout << " " << city;
    }
    cout << endl;
    
    vector<pair<int, int>> edges;
    for (int i = 0; i < l - 1; i++) {
        int a, b;
        cin >> a >> b;
        edges.push_back({a, b});
    }
    
    return edges;
}

// 戦略的に都市を選んでクエリを実行
vector<int> selectCitiesForQuery(int N, int L, const vector<vector<double>>& dist_matrix, const vector<vector<bool>>& is_real_dist) {
    // できるだけ離れた都市を選ぶ
    vector<int> selected;
    vector<bool> used(N, false);
    
    // 未知の距離情報が多い都市を優先的に選ぶ
    vector<pair<int, int>> unknown_count;
    for (int i = 0; i < N; i++) {
        int count = 0;
        for (int j = 0; j < N; j++) {
            if (i != j && !is_real_dist[i][j]) {
                count++;
            }
        }
        unknown_count.push_back({count, i});
    }
    
    // 未知の距離情報が多い順にソート
    sort(unknown_count.begin(), unknown_count.end(), greater<pair<int, int>>());
    
    // 最初の都市を選択
    int first = unknown_count[0].second;
    selected.push_back(first);
    used[first] = true;
    
    // 残りの都市を選択（できるだけ離れた点を選ぶ）
    while (selected.size() < L) {
        double max_score = -1;
        int next = -1;
        
        for (int i = 0; i < N; i++) {
            if (used[i]) continue;
            
            double min_dist = 1e9;
            for (int j : selected) {
                min_dist = min(min_dist, dist_matrix[i][j]);
            }
            
            // 未知の距離情報の数を考慮
            int unknown = 0;
            for (int j = 0; j < N; j++) {
                if (i != j && !is_real_dist[i][j]) {
                    unknown++;
                }
            }
            
            // 距離と未知情報の両方を考慮したスコア
            double score = min_dist * (1.0 + unknown / (double)N);
            
            if (score > max_score) {
                max_score = score;
                next = i;
            }
        }
        
        selected.push_back(next);
        used[next] = true;
    }
    
    return selected;
}

// K-meansクラスタリングを行う
vector<vector<int>> kMeansClustering(const vector<vector<double>>& dist_matrix, int N, int M, const vector<int>& G) {
    // 初期クラスタ中心をランダムに選択
    vector<int> centers;
    vector<bool> used(N, false);
    
    // 最初のセンターをランダムに選択
    int first_center = rng() % N;
    centers.push_back(first_center);
    used[first_center] = true;
    
    // 残りのセンターを選択（できるだけ離れた点を選ぶ）
    for (int i = 1; i < M; i++) {
        double max_min_dist = -1;
        int next_center = -1;
        
        for (int j = 0; j < N; j++) {
            if (used[j]) continue;
            
            double min_dist = 1e9;
            for (int center : centers) {
                min_dist = min(min_dist, dist_matrix[j][center]);
            }
            
            if (min_dist > max_min_dist) {
                max_min_dist = min_dist;
                next_center = j;
            }
        }
        
        centers.push_back(next_center);
        used[next_center] = true;
    }
    
    // クラスタリングを実行
    vector<vector<int>> clusters(M);
    vector<int> assignments(N, -1);
    
    // 各クラスタに最初の都市を割り当て
    for (int i = 0; i < M; i++) {
        clusters[i].push_back(centers[i]);
        assignments[centers[i]] = i;
    }
    
    // 残りの都市を割り当て
    vector<int> remaining;
    for (int i = 0; i < N; i++) {
        if (assignments[i] == -1) {
            remaining.push_back(i);
        }
    }
    
    // 距離に基づいてソート
    sort(remaining.begin(), remaining.end(), [&](int a, int b) {
        double min_dist_a = 1e9, min_dist_b = 1e9;
        for (int center : centers) {
            min_dist_a = min(min_dist_a, dist_matrix[a][center]);
            min_dist_b = min(min_dist_b, dist_matrix[b][center]);
        }
        return min_dist_a < min_dist_b;
    });
    
    // 各クラスタのサイズ制約を満たすように割り当て
    for (int city : remaining) {
        // 最も近いクラスタを見つける
        int closest_cluster = -1;
        double min_dist = 1e9;
        
        for (int c = 0; c < M; c++) {
            if (clusters[c].size() >= G[c]) continue; // クラスタが既に満杯
            
            double avg_dist = 0;
            for (int member : clusters[c]) {
                avg_dist += dist_matrix[city][member];
            }
            avg_dist /= clusters[c].size();
            
            if (avg_dist < min_dist) {
                min_dist = avg_dist;
                closest_cluster = c;
            }
        }
        
        if (closest_cluster != -1) {
            clusters[closest_cluster].push_back(city);
            assignments[city] = closest_cluster;
        }
    }
    
    // クラスタサイズの調整
    bool changed = true;
    int iterations = 0;
    const int MAX_ITERATIONS = 100; // 最大反復回数を制限
    
    while (changed && iterations < MAX_ITERATIONS) {
        changed = false;
        iterations++;
        
        for (int c = 0; c < M; c++) {
            while (clusters[c].size() < G[c]) {
                // 他のクラスタから都市を移動
                int best_city = -1;
                int source_cluster = -1;
                double min_increase = 1e9;
                
                for (int other_c = 0; other_c < M; other_c++) {
                    if (other_c == c || clusters[other_c].size() <= G[other_c]) continue;
                    
                    for (int city : clusters[other_c]) {
                        double current_dist = 0;
                        for (int member : clusters[other_c]) {
                            if (member != city) {
                                current_dist += dist_matrix[city][member];
                            }
                        }
                        
                        double new_dist = 0;
                        for (int member : clusters[c]) {
                            new_dist += dist_matrix[city][member];
                        }
                        
                        double increase = new_dist - current_dist;
                        if (increase < min_increase) {
                            min_increase = increase;
                            best_city = city;
                            source_cluster = other_c;
                        }
                    }
                }
                
                if (best_city != -1) {
                    // 都市を移動
                    clusters[c].push_back(best_city);
                    clusters[source_cluster].erase(
                        remove(clusters[source_cluster].begin(), clusters[source_cluster].end(), best_city),
                        clusters[source_cluster].end()
                    );
                    assignments[best_city] = c;
                    changed = true;
                } else {
                    break;
                }
            }
        }
    }
    
    // 最終的なクラスタサイズを確認
    bool size_ok = true;
    for (int c = 0; c < M; c++) {
        if (clusters[c].size() != G[c]) {
            size_ok = false;
            break;
        }
    }
    
    // サイズ制約を満たせない場合は、ランダムに割り当て直す
    if (!size_ok) {
        // 複数回試行
        double best_dist = 1e18;
        vector<vector<int>> best_clusters;
        
        for (int attempt = 0; attempt < 5; attempt++) {
            vector<int> all_cities(N);
            iota(all_cities.begin(), all_cities.end(), 0);
            shuffle(all_cities.begin(), all_cities.end(), rng);
            
            vector<vector<int>> new_clusters(M);
            int city_idx = 0;
            for (int c = 0; c < M; c++) {
                for (int i = 0; i < G[c]; i++) {
                    new_clusters[c].push_back(all_cities[city_idx++]);
                }
            }
            
            // クラスタ内の総距離を計算
            double total_dist = 0;
            for (int c = 0; c < M; c++) {
                for (int i = 0; i < new_clusters[c].size(); i++) {
                    for (int j = i + 1; j < new_clusters[c].size(); j++) {
                        total_dist += dist_matrix[new_clusters[c][i]][new_clusters[c][j]];
                    }
                }
            }
            
            // より良い解を採用
            if (attempt == 0 || total_dist < best_dist) {
                best_clusters = new_clusters;
                best_dist = total_dist;
            }
        }
        
        clusters = best_clusters;
    }
    
    return clusters;
}

// 最小全域木を構築
vector<pair<int, int>> buildMST(const vector<int>& cities, const vector<vector<double>>& dist_matrix, const vector<vector<bool>>& is_real_dist) {
    int n = cities.size();
    vector<Edge> edges;
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int u = cities[i];
            int v = cities[j];
            edges.push_back({u, v, dist_matrix[u][v], is_real_dist[u][v]});
        }
    }
    
    sort(edges.begin(), edges.end());
    
    UnionFind uf(1000); // 十分大きな値
    vector<pair<int, int>> mst;
    
    for (const Edge& edge : edges) {
        if (!uf.same(edge.u, edge.v)) {
            uf.unite(edge.u, edge.v);
            mst.push_back({edge.u, edge.v});
            if (mst.size() == n - 1) break;
        }
    }
    
    return mst;
}

// グループ内の総距離を計算
double calculateGroupDistance(const vector<int>& group, const vector<vector<double>>& dist_matrix) {
    double total_dist = 0;
    for (int i = 0; i < group.size(); i++) {
        for (int j = i + 1; j < group.size(); j++) {
            total_dist += dist_matrix[group[i]][group[j]];
        }
    }
    return total_dist;
}

int main() {
    int N, M, Q, L, W;
    cin >> N >> M >> Q >> L >> W;
    
    vector<int> G(M);
    for (int i = 0; i < M; i++) {
        cin >> G[i];
    }
    
    vector<CityRange> city_ranges(N);
    for (int i = 0; i < N; i++) {
        cin >> city_ranges[i].lx >> city_ranges[i].rx >> city_ranges[i].ly >> city_ranges[i].ry;
    }
    
    // 都市の中心座標を推定
    vector<pair<double, double>> estimated_centers(N);
    for (int i = 0; i < N; i++) {
        estimated_centers[i] = estimateCityCenter(city_ranges[i]);
    }
    
    // 都市間の距離行列を初期化（推定値）
    vector<vector<double>> dist_matrix(N, vector<double>(N, 0));
    vector<vector<bool>> is_real_dist(N, vector<bool>(N, false));
    
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist = estimateDistance(estimated_centers[i], estimated_centers[j]);
            dist_matrix[i][j] = dist_matrix[j][i] = dist;
        }
    }
    
    // 占いを使って距離行列を更新
    int queries_used = 0;
    
    // 戦略的に都市を選んでクエリを実行
    while (queries_used < Q / 2) {
        vector<int> selected_cities = selectCitiesForQuery(N, L, dist_matrix, is_real_dist);
        
        // クエリを実行
        vector<pair<int, int>> mst_edges = query(selected_cities);
        queries_used++;
        
        // 距離行列を更新
        updateDistanceMatrix(dist_matrix, is_real_dist, selected_cities, mst_edges);
    }
    
    // K-meansクラスタリングでグループ分け
    vector<vector<int>> groups = kMeansClustering(dist_matrix, N, M, G);
    
    // 各グループ内で最小全域木を構築
    vector<vector<pair<int, int>>> group_edges(M);
    
    // グループのサイズが大きい順にクエリを使用
    vector<pair<int, int>> group_sizes;
    for (int g = 0; g < M; g++) {
        group_sizes.push_back({groups[g].size(), g});
    }
    sort(group_sizes.begin(), group_sizes.end(), greater<pair<int, int>>());
    
    for (const auto& p : group_sizes) {
        int size = p.first;
        int g = p.second;
        if (size <= 1) continue;
        
        // 残りのクエリを使って、各グループ内の最小全域木を取得
        if (queries_used < Q && size <= L) {
            vector<pair<int, int>> mst_edges = query(groups[g]);
            queries_used++;
            group_edges[g] = mst_edges;
        } else {
            // クエリを使い切った場合は、推定距離に基づいて最小全域木を構築
            group_edges[g] = buildMST(groups[g], dist_matrix, is_real_dist);
        }
    }
    
    // 結果を出力
    cout << "!" << endl;
    
    for (int g = 0; g < M; g++) {
        // グループに属する都市を出力
        for (int i = 0; i < groups[g].size(); i++) {
            cout << groups[g][i];
            if (i < groups[g].size() - 1) cout << " ";
        }
        cout << endl;
        
        // グループ内の道路を出力
        for (const auto& edge : group_edges[g]) {
            cout << edge.first << " " << edge.second << endl;
        }
    }
    
    return 0;
}
