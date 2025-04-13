#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <functional>
using namespace std;

// グローバル変数
const int N = 50;
int si, sj;
vector<vector<int>> tile(N, vector<int>(N));
vector<vector<int>> point(N, vector<int>(N));
vector<vector<bool>> visited(N, vector<bool>(N, false));
set<int> visited_tiles;

// 方向
const int dx[4] = {-1, 1, 0, 0};
const int dy[4] = {0, 0, -1, 1};
const char dir_char[4] = {'U', 'D', 'L', 'R'};

// 移動可能かチェック
bool can_move(int i, int j) {
    // 盤面外
    if (i < 0 || i >= N || j < 0 || j >= N) {
        return false;
    }
    
    // 既に訪れたタイル
    if (visited_tiles.count(tile[i][j])) {
        return false;
    }
    
    return true;
}

// 前方宣言
int calculate_score(const string& path);
string modify_path(const string& path);

// 経路の一部を破壊して再構築する
string reconstruct_path(const string& path) {
    if (path.size() < 10) return path;
    
    random_device rd;
    mt19937 gen(rd());
    
    // 破壊する経路の開始位置と長さを決定
    uniform_int_distribution<> start_dist(0, path.size() - 10);
    int start_pos = start_dist(gen);
    uniform_int_distribution<> length_dist(5, min(20, (int)path.size() - start_pos));
    int destroy_len = length_dist(gen);
    
    // 経路を分割
    string prefix = path.substr(0, start_pos);
    string suffix = (start_pos + destroy_len < path.size()) ? path.substr(start_pos + destroy_len) : "";
    
    // 分割点の状態を計算
    int i = si, j = sj;
    set<int> used_tiles = {tile[i][j]};
    
    for (int idx = 0; idx < start_pos; idx++) {
        char c = path[idx];
        int d;
        if (c == 'U') d = 0;
        else if (c == 'D') d = 1;
        else if (c == 'L') d = 2;
        else if (c == 'R') d = 3;
        else continue;
        
        i += dx[d];
        j += dy[d];
        used_tiles.insert(tile[i][j]);
    }
    
    // 次の分割点の状態を計算
    int ni = i, nj = j;
    for (int idx = start_pos + destroy_len; idx < path.size(); idx++) {
        char c = path[idx];
        int d;
        if (c == 'U') d = 0;
        else if (c == 'D') d = 1;
        else if (c == 'L') d = 2;
        else if (c == 'R') d = 3;
        else continue;
        
        ni += dx[d];
        nj += dy[d];
    }
    
    // 分割点間を繋ぐ経路を探索
    int max_attempts = 100;
    int best_score = -1;
    string best_middle = "";
    
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        string middle = "";
        int middle_score = 0;
        int curr_i = i, curr_j = j;
        set<int> curr_used_tiles = used_tiles;
        bool valid = true;
        
        // ランダムな経路で繋ぐ
        for (int step = 0; step < destroy_len * 2 && curr_i != ni && curr_j != nj; step++) {
            // ランダムな方向を選ぶ
            vector<int> valid_dirs;
            for (int d = 0; d < 4; d++) {
                int next_i = curr_i + dx[d];
                int next_j = curr_j + dy[d];
                
                // 盤面外
                if (next_i < 0 || next_i >= N || next_j < 0 || next_j >= N) {
                    continue;
                }
                
                // 既に訪れたタイル
                if (curr_used_tiles.count(tile[next_i][next_j])) {
                    continue;
                }
                
                valid_dirs.push_back(d);
            }
            
            if (valid_dirs.empty()) {
                valid = false;
                break;
            }
            
            uniform_int_distribution<> dir_dist(0, valid_dirs.size() - 1);
            int d = valid_dirs[dir_dist(gen)];
            
            curr_i += dx[d];
            curr_j += dy[d];
            middle += dir_char[d];
            curr_used_tiles.insert(tile[curr_i][curr_j]);
            middle_score += point[curr_i][curr_j];
        }
        
        // 有効な経路で最高スコアを更新
        if (valid && middle_score > best_score) {
            best_score = middle_score;
            best_middle = middle;
        }
    }
    
    // 再構築した経路を返す
    if (!best_middle.empty()) {
        string new_path = prefix + best_middle + suffix;
        int new_score = calculate_score(new_path);
        int old_score = calculate_score(path);
        
        if (new_score > old_score) {
            return new_path;
        }
    }
    
    return path;
}

// 貪欲法による解法
string greedy_solve() {
    string path = "";
    int i = si, j = sj;
    int total_score = point[i][j];
    
    // 初期位置のタイルを訪問済みに
    visited[i][j] = true;
    visited_tiles.insert(tile[i][j]);
    
    while (true) {
        int best_dir = -1;
        int best_score = -1;
        
        // 4方向を試す
        for (int d = 0; d < 4; d++) {
            int ni = i + dx[d];
            int nj = j + dy[d];
            
            if (can_move(ni, nj)) {
                // より高いスコアの方向を選択
                if (point[ni][nj] > best_score) {
                    best_score = point[ni][nj];
                    best_dir = d;
                }
            }
        }
        
        // 移動できる方向がなければ終了
        if (best_dir == -1) {
            break;
        }
        
        // 移動
        i += dx[best_dir];
        j += dy[best_dir];
        path += dir_char[best_dir];
        
        // 訪問済みに
        visited[i][j] = true;
        visited_tiles.insert(tile[i][j]);
        total_score += point[i][j];
    }
    
    return path;
}

// ビームサーチによる改良解法
string beam_search(int beam_width = 100, int max_steps = 1000) {
    // 状態を表す構造体
    struct State {
        int i, j;
        int score;
        string path;
        set<int> visited_tiles;
        
        // スコアによる降順ソート
        bool operator<(const State& other) const {
            return score > other.score;
        }
    };
    
    // 初期状態
    vector<State> beam = {
        {si, sj, point[si][sj], "", {tile[si][sj]}}
    };
    
    string best_path = "";
    int best_score = point[si][sj];
    
    for (int step = 0; step < max_steps; step++) {
        vector<State> next_beam;
        
        for (const auto& state : beam) {
            // 各方向に移動
            for (int d = 0; d < 4; d++) {
                int ni = state.i + dx[d];
                int nj = state.j + dy[d];
                
                // 盤面外
                if (ni < 0 || ni >= N || nj < 0 || nj >= N) {
                    continue;
                }
                
                // 既に訪れたタイル
                if (state.visited_tiles.count(tile[ni][nj])) {
                    continue;
                }
                
                // 新しい状態を作成
                State new_state = state;
                new_state.i = ni;
                new_state.j = nj;
                new_state.score += point[ni][nj];
                new_state.path += dir_char[d];
                new_state.visited_tiles.insert(tile[ni][nj]);
                
                next_beam.push_back(new_state);
                
                // 最良スコアを更新
                if (new_state.score > best_score) {
                    best_score = new_state.score;
                    best_path = new_state.path;
                }
            }
        }
        
        // 次のビームを選定
        if (next_beam.empty()) {
            break;
        }
        
        sort(next_beam.begin(), next_beam.end());
        beam = vector<State>(next_beam.begin(), next_beam.begin() + min(beam_width, (int)next_beam.size()));
    }
    
    return best_path;
}


// 前方宣言
string reconstruct_path(const string& path);

// 焼きなまし法による改良
string simulated_annealing(string initial_path, int max_iter = 10000, double temp_start = 100.0, double temp_end = 0.1) {
    // 乱数生成器
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> uniform_dist(0.0, 1.0);
    
    string current_path = initial_path;
    int current_score = calculate_score(current_path);
    
    string best_path = current_path;
    int best_score = current_score;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // 温度
        double temperature = temp_start * pow(temp_end / temp_start, (double)iter / max_iter);
        
        // パスの一部を変更またはパスの再構築
        string new_path;
        if (uniform_dist(gen) < 0.3) {  // 30%の確率で経路の再構築を行う
            new_path = reconstruct_path(current_path);
        } else {
            new_path = modify_path(current_path);
        }
        int new_score = calculate_score(new_path);
        
        // スコアの差分
        int delta = new_score - current_score;
        
        // 採用判定
        if (delta > 0 || uniform_dist(gen) < exp(delta / temperature)) {
            current_path = new_path;
            current_score = new_score;
            
            // 最良解の更新
            if (current_score > best_score) {
                best_path = current_path;
                best_score = current_score;
            }
        }
    }
    
    return best_path;
}

// パスのスコアを計算
int calculate_score(const string& path) {
    int i = si, j = sj;
    int score = point[i][j];
    
    set<int> used_tiles = {tile[i][j]};
    
    for (char c : path) {
        int d;
        if (c == 'U') d = 0;
        else if (c == 'D') d = 1;
        else if (c == 'L') d = 2;
        else if (c == 'R') d = 3;
        else continue;
        
        i += dx[d];
        j += dy[d];
        
        // 盤面外ならパスを途中で終了
        if (i < 0 || i >= N || j < 0 || j >= N) {
            break;
        }
        
        // 既に使用したタイルならパスを途中で終了
        if (used_tiles.count(tile[i][j])) {
            break;
        }
        
        used_tiles.insert(tile[i][j]);
        score += point[i][j];
    }
    
    return score;
}

// パスの一部を変更
string modify_path(const string& path) {
    if (path.empty()) {
        return path;
    }
    
    random_device rd;
    mt19937 gen(rd());
    
    // 変更方法をランダムに選択
    uniform_int_distribution<> method_dist(0, 2);
    int method = method_dist(gen);
    
    if (method == 0 && path.size() > 1) {
        // 末尾を削除
        uniform_int_distribution<> length_dist(1, min(5, (int)path.size()));
        int cut_length = length_dist(gen);
        return path.substr(0, path.size() - cut_length);
    } 
    else if (method == 1) {
        // パスの一部を変更
        if (path.size() < 3) return path;
        
        uniform_int_distribution<> pos_dist(0, path.size() - 3);
        int pos = pos_dist(gen);
        uniform_int_distribution<> length_dist(1, min(3, (int)path.size() - pos));
        int change_length = length_dist(gen);
        
        string new_path = path.substr(0, pos);
        
        // 変更部分をランダムに生成
        uniform_int_distribution<> dir_dist(0, 3);
        for (int i = 0; i < change_length; i++) {
            new_path += dir_char[dir_dist(gen)];
        }
        
        // 残りのパスをそのまま使う（有効性は calculate_score で判定）
        if (pos + change_length < path.size()) {
            new_path += path.substr(pos + change_length);
        }
        
        return new_path;
    }
    else {
        // 末尾に追加
        uniform_int_distribution<> length_dist(1, 5);
        int add_length = length_dist(gen);
        
        string new_path = path;
        uniform_int_distribution<> dir_dist(0, 3);
        
        for (int i = 0; i < add_length; i++) {
            new_path += dir_char[dir_dist(gen)];
        }
        
        return new_path;
    }
}

int main() {
    // 入力
    cin >> si >> sj;
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> tile[i][j];
        }
    }
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> point[i][j];
        }
    }
    
    // 時間測定開始
    auto start_time = chrono::high_resolution_clock::now();
    
    // 初期解（貪欲法）
    string greedy_path = greedy_solve();
    
    // ビームサーチで改良
    string beam_path = beam_search(50, 500);
    
    // 初期解として、貪欲法とビームサーチの良い方を採用
    string best_path;
    int greedy_score = calculate_score(greedy_path);
    int beam_score = 0;
    if (!beam_path.empty()) {
        beam_score = calculate_score(beam_path);
    }
    
    if (beam_score > greedy_score) {
        best_path = beam_path;
    } else {
        best_path = greedy_path;
    }
    
    // 焼きなまし法で改良を行う
    auto current_time = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count();
    
    // 残り時間で焼きなまし
    if (elapsed < 900) {  // 900ms以内で計算できるように調整
        int sa_iterations = min(20000, (int)(100000 * (900 - elapsed) / 900));
        double start_temp = 100.0;
        double end_temp = 0.01;
        best_path = simulated_annealing(best_path, sa_iterations, start_temp, end_temp);
    }
    
    // 余裕があれば経路の再構築を試みる
    current_time = chrono::high_resolution_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count();
    
    if (elapsed < 950) {
        for (int i = 0; i < 50; i++) {
            string new_path = reconstruct_path(best_path);
            int new_score = calculate_score(new_path);
            int best_score = calculate_score(best_path);
            if (new_score > best_score) {
                best_path = new_path;
            }
        }
    }
    
    // 出力
    cout << best_path << endl;
    
    return 0;
}
