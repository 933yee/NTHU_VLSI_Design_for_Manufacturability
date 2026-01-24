#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <climits>
#include <random>
#include <utility>
#include <cmath>
#include <chrono>
#include <stack>

using namespace std;
using namespace std::chrono;


struct Staple {
    int x, y;
};

struct CellType {
    CellType(int id, int width, int height) : id(id), width(width), height(height) {}
    CellType() {}
    int id;
    int width;
    int height;
    unordered_set<int> pins;
};

struct Cell {
    Cell(int id, int type_id, int x, int y, int max_displacement) 
        : id(id), type_id(type_id), x(x), y(y), max_displacement(max_displacement) {}
    Cell() {}
    int id;
    int type_id;
    int x;
    int y;
    int max_displacement;
};

struct CompactKey {
    int i;
    int a1, b1;
    int a2, b2;
    int a3, b3;

    bool operator==(const CompactKey& o) const {
        return tie(i, a1, b1, a2, b2, a3, b3) == tie(o.i, o.a1, o.b1, o.a2, o.b2, o.a3, o.b3);
    }

    string to_string() const {
        return 
               std::to_string(a1) + "_" + std::to_string(b1) + "_" +
               std::to_string(a2) + "_" + std::to_string(b2) + "_" +
               std::to_string(a3) + "_" + std::to_string(b3);
    }
};

namespace std {
    template<>
    struct hash<CompactKey> {
        size_t operator()(const CompactKey& k) const {
            return hash<string>()(k.to_string());
        }
    };
}

struct MATRONode {
    int i, s1, l1, s2, l2, s3, l3;
    vector<int> benefit_case; // benefit_case[5]
    vector<MATRONode*> prev_case; // prev_case[5]
    vector<int> prev_case_index; // prev_case_index[5]

    MATRONode(int i, int s1, int l1, int s2, int l2, int s3, int l3)
        : i(i), s1(s1), l1(l1), s2(s2), l2(l2), s3(s3), l3(l3),
          benefit_case(5, INT_MIN), prev_case(5, nullptr), prev_case_index(5, -1) {}
};

vector<string> split(string& str) {
    vector<string> result;
    istringstream iss(str);
    for (string s; iss >> s;) {
        result.push_back(s);
    }
    return result;
}

vector<string> input_preprocess(ifstream &input_file) {
    vector<string> objects;
    string line;
    while (getline(input_file, line)) 
        objects.push_back(line);
    return objects;
}

int chip_left_boundary = 0, chip_right_boundary = 0, chip_top_boundary = 0, chip_bottom_boundary = 0;
int number_of_rows = 0, row_height = 0, site_width = 0, number_of_sites_per_row = 0,
    number_of_cell_types = 0, number_of_cells = 0;

unordered_map<int, CellType*> cell_type_map;
unordered_map<int, Cell*> cell_map;
vector<vector<Cell*>> row_cells;
vector<vector<bool>> occupied, has_pin;
vector<Staple> staples;
vector<vector<int>> original_site_cell_map;
queue<MATRONode*> Q;
vector<MATRONode*> best_path;
vector<int> best_case_used; 
vector<vector<MATRONode*>> all_paths;
vector<vector<int>> all_cases;

unordered_map<CompactKey, MATRONode*> matro_node_table;
vector<MATRONode*> all_created_nodes;

// time
std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

MATRONode* create_or_get_node(int row_id, int i, int s1, int l1, int s2, int l2, int s3, int l3, queue<MATRONode*>* Q) {

    // 這個 cell 原本的 x 座標
    Cell* cell1 = (s1 >= 0 && row_id < number_of_rows) ? row_cells[row_id][s1] : nullptr;
    Cell* cell2 = (s2 >= 0 && row_id + 1 < number_of_rows) ? row_cells[row_id + 1][s2] : nullptr;
    Cell* cell3 = (s3 >= 0 && row_id + 2 < number_of_rows) ? row_cells[row_id + 2][s3] : nullptr;

    int origin_x1 = cell1 ? cell1->x / site_width : 0;
    int origin_x2 = cell2 ? cell2->x / site_width : 0;
    int origin_x3 = cell3 ? cell3->x / site_width : 0;


    int origin_site_cell_id_1 = row_id < number_of_rows ? original_site_cell_map[row_id][i] : -1;
    int origin_site_cell_id_2 = row_id + 1 < number_of_rows ? original_site_cell_map[row_id + 1][i] : -1;
    int origin_site_cell_id_3 = row_id + 2 < number_of_rows ? original_site_cell_map[row_id + 2][i] : -1;

    CompactKey key;
    key.i = i;
    key.a1 = cell1 ? (cell1->id - origin_site_cell_id_1) : -999;  // -999 just to distinguish
    key.a2 = cell2 ? (cell2->id - origin_site_cell_id_2) : -999;
    key.a3 = cell3 ? (cell3->id - origin_site_cell_id_3) : -999;
    key.b1 = cell1 ? (i - l1 - origin_x1) : -999;
    key.b2 = cell2 ? (i - l2 - origin_x2) : -999;
    key.b3 = cell3 ? (i - l3 - origin_x3) : -999;


    if (matro_node_table.count(key)) return matro_node_table[key];
    MATRONode* new_node = new MATRONode(i, s1, l1, s2, l2, s3, l3);
    matro_node_table[key] = new_node;
    Q->push(new_node); // Add to queue if provided
    all_created_nodes.push_back(new_node); // Store all created nodes for later use

    // cout key
    return new_node;
}

bool is_deferrable(int row_id, int cell_id, int i) {
    // Check if the site is deferrable
    if(row_id >= number_of_rows) return true;
    if (cell_id + 1 >= row_cells[row_id].size()) return true;

    Cell* cell = row_cells[row_id][cell_id + 1];
    
    return cell->max_displacement > i * site_width - cell->x;
}

bool is_placeable(int row_id, int cell_id, int i, int l) {
    // Check if the site is placeable
    if (row_id >= number_of_rows) return false;
    if (cell_id + 1 >= row_cells[row_id].size()) return false;
    if (cell_id >= 0) {
        // 看會不會與上一個 cell 撞到
        Cell* old_cell = row_cells[row_id][cell_id];
        int old_width_sites = (cell_type_map[old_cell->type_id]->width / site_width);
        if (l < old_width_sites) return false;
    }

    // 看能不能放下一個 cell
    Cell* new_cell = row_cells[row_id][cell_id + 1];

    if(new_cell->max_displacement < abs(i * site_width - new_cell->x))
        return false;

    int new_width_sites = (cell_type_map[new_cell->type_id]->width / site_width);
    return i + new_width_sites <= number_of_sites_per_row;
}

bool would_violate_antiparallel(int prev_case, int cur_case, int i, int row_id) {
    if ((prev_case == 1 && cur_case == 3) || (prev_case == 3 && cur_case == 1)) {
        return true;
    }
    if(cur_case == 1 || cur_case == 4){
        if(row_id - 2 >= 0 && i - 1 >= 0){
            int prev_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - (i - 1) - 1];
            int curr_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - i - 1];
            if (prev_case_prev_row == 2 && curr_case_prev_row != 2 && prev_case != 1) return true;
        }
        if(row_id - 2 >= 0 && i + 1 < number_of_sites_per_row){
            int prev_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - (i + 1) - 1];
            int curr_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - i - 1];
            if (prev_case_prev_row == 2 && curr_case_prev_row != 2) return true;
        }
    }
    if(cur_case == 2){
        if(row_id - 2 >= 0 && i - 1 >= 0){
            int prev_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - (i - 1) - 1];
            int curr_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - i - 1];
            if ((prev_case_prev_row == 3 || prev_case_prev_row == 4) && curr_case_prev_row != 3 && curr_case_prev_row != 4 && prev_case != 2) return true;
        }
        if(row_id - 2 >= 0 && i + 1 < number_of_sites_per_row){
            int prev_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - (i + 1) - 1];
            int curr_case_prev_row = all_cases[(row_id - 2) / 3][number_of_sites_per_row - i - 1];
            if ((prev_case_prev_row == 3 || prev_case_prev_row == 4) && curr_case_prev_row != 3 && curr_case_prev_row != 4) return true;
        }
    }
    
    return false;
}

bool can_insert_staple(int i, int row_id, int case_k, MATRONode* u) {
    int r1 = row_id;
    int r2 = row_id + 1;
    int r3 = row_id + 2;
    int r_prev3 = row_id - 1;

    auto no_cell_or_pin = [&](int row, int s, int l, bool is_prev3 = false) {
        if(row < 0 || row >= number_of_rows) return false;

        if(is_prev3) {
            int prev_case = all_cases[row / 3][number_of_sites_per_row - i - 1];
            if(prev_case == 3 || prev_case == 4) return false;

            MATRONode* prev_node = all_paths[row / 3][number_of_sites_per_row - i - 1];
            int s3 = prev_node->s3;
            int l3 = prev_node->l3;
            if(s3 < 0) return true;
            Cell* prev_cell = row_cells[row][s3];
            CellType* prev_type = cell_type_map[prev_cell->type_id];
            int start_site = i - l3;
            for (int pin_offset : prev_type->pins) {
                if (start_site + pin_offset == i - 1)
                    return false;
            }
            return true;
        }

        if (s < 0 || s >= row_cells[row].size()) return true;
        Cell* cell = row_cells[row][s];
        CellType* type = cell_type_map[cell->type_id];

        int width_sites = type->width / site_width;
        int start_site = i - l; // where this cell starts based on current site

        // 如果這 cell 橫跨了 site i 就會撞到
        // if (i >= start_site && i <= start_site + width_sites)
        //     return false;

        // 如果 cell 的某個 pin 剛好落在 site i，也不行
        for (int pin_offset : type->pins) {
            if (start_site + pin_offset == i - 1)
                return false;
        }

        return true;
    };

    switch (case_k) {
        case 0: return true;
        case 1: return no_cell_or_pin(r_prev3, u->s1, u->l1, true) && no_cell_or_pin(r1, u->s1, u->l1);
        case 2: return no_cell_or_pin(r1, u->s1, u->l1) && no_cell_or_pin(r2, u->s2, u->l2);
        case 3: return no_cell_or_pin(r2, u->s2, u->l2) && no_cell_or_pin(r3, u->s3, u->l3);
        case 4: return no_cell_or_pin(r_prev3, u->s1, u->l1, true) &&
                      no_cell_or_pin(r1,       u->s1, u->l1) &&
                      no_cell_or_pin(r2,       u->s2, u->l2) &&
                      no_cell_or_pin(r3,       u->s3, u->l3);
        default: return false;
    }
}



void update_benefit(MATRONode* u, MATRONode* v, int i, int row_id) {
    for (int cur_case = 0; cur_case < 5; ++cur_case) {
        if (!can_insert_staple(i, row_id, cur_case, v)) continue;

        for (int prev_case = 0; prev_case < 5; ++prev_case) {
            if (would_violate_antiparallel(prev_case, cur_case, i, row_id)) continue;

            int delta = (cur_case == 0 ? 0 : (cur_case == 4 ? 2 : 1));
            int candidate_benefit = u->benefit_case[prev_case] + delta;

            if (candidate_benefit > v->benefit_case[cur_case]) {
                v->benefit_case[cur_case] = candidate_benefit;
                v->prev_case[cur_case] = u;
                v->prev_case_index[cur_case] = prev_case;
            }
        }
    }
}


// Triple-row DP implementation for staple insertion.
// Insert this function after your draw_svg() and before main()

void reset() {
    unordered_set<MATRONode*> keep(best_path.begin(), best_path.end());
    for (MATRONode* node : all_created_nodes) {
        if (keep.find(node) == keep.end()) delete node;
    }
    all_created_nodes.clear();
    matro_node_table.clear();
    best_path.clear();
    best_case_used.clear();
    while (!Q.empty()) Q.pop();
}


void triple_row_dp(int row_id) {
    const int row1 = row_id;
    // cout << "row_id: " << row_id << endl;
    const int row2 = row_id + 1;
    const int row3 = row_id + 2;
    const int row_prev3 = row_id - 1; // new: previous bottom row to check anti-parallel violations


    // Initialize source node (no cells placed yet)
    MATRONode* source = create_or_get_node(row_id, 0, -1, 0, -1, 0, -1, 0, &Q);
    for (int i = 0; i < 5; ++i) source->benefit_case[i] = 0;

    // cout << "Source node: " << source->i << ", " << source->s1 << ", " << source->l1 << ", "
        //  << source->s2 << ", " << source->l2 << ", " << source->s3 << ", " << source->l3 << endl;
    // cout << "number_of_sites_per_row = " << number_of_sites_per_row << endl;
    // time
    auto start = high_resolution_clock::now();
    int t = 0;
    while (!Q.empty()) {
        if(t == number_of_sites_per_row) break;
        t++;
        // cout << "Processing node: " << t++ << endl;
        matro_node_table.clear();
        long long size = Q.size();
        for(long long j = 0; j < size; j++) {
            MATRONode* u = Q.front();
            Q.pop();
            // cout << "Processing node: " << u->i << ", " << u->s1 << ", " << u->l1 << ", "
            //     << u->s2 << ", " << u->l2 << ", " << u->s3 << ", " << u->l3 << endl;
            int i = u->i;
            // cout << "i: " << i << endl;
            if (i >= number_of_sites_per_row) break;
            // 啥都不做
            // cout << "Begin checking no cell" << endl;
            if(is_deferrable(row_id, u->s1, i) && is_deferrable(row_id+1, u->s2, i) && is_deferrable(row_id+2, u->s3, i)) {
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1, u->l1 + 1, u->s2, u->l2 + 1, u->s3, u->l3 + 1, &Q);
                // cout << "Cell created: " << v->i << ", " << v->s1 << ", " << v->l1 << ", "
                    //  << v->s2 << ", " << v->l2 << ", " << v->s3 << ", " << v->l3 << endl;
                update_benefit(u, v, i, row_id);
            }

            // cout << "End checking no cell" << endl;

            // 放 cell 到 R1
            // cout << "Begin checking cell in R1" << endl;
            if (is_placeable(row_id, u->s1, i, u->l1) && is_deferrable(row_id+1, u->s2, i) && is_deferrable(row_id+2, u->s3, i)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1 + 1, 1, u->s2, u->l2 + 1, u->s3, u->l3 + 1, &Q);
                update_benefit(u, v, i, row_id);
            }
            // cout << "End checking cell in R1" << endl;
            
            // 放 cell 到 R2
            if (is_deferrable(row_id, u->s1, i) && is_placeable(row_id + 1, u->s2, i, u->l2) && is_deferrable(row_id + 2, u->s3, i)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1, u->l1 + 1, u->s2 + 1, 1, u->s3, u->l3 + 1, &Q);
                update_benefit(u, v, i, row_id);
            }

            // 放 cell 到 R3
            if (is_deferrable(row_id, u->s1, i) && is_deferrable(row_id + 1, u->s2, i) && is_placeable(row_id + 2, u->s3, i, u->l3)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1, u->l1 + 1, u->s2, u->l2 + 1, u->s3 + 1, 1, &Q);
                update_benefit(u, v, i, row_id);
            }

            // 放 cell 到 R1, R2
            if (is_placeable(row_id, u->s1, i, u->l1) && is_placeable(row_id + 1, u->s2, i, u->l2) && is_deferrable(row_id + 2, u->s3, i)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1 + 1, 1, u->s2 + 1, 1, u->s3, u->l3 + 1, &Q);
                update_benefit(u, v, i, row_id);
            }

            // 放 cell 到 R1, R3
            if (is_placeable(row_id, u->s1, i, u->l1) && is_deferrable(row_id + 1, u->s2, i) && is_placeable(row_id + 2, u->s3, i, u->l3)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1 + 1, 1, u->s2, u->l2 + 1, u->s3 + 1, 1, &Q);
                update_benefit(u, v, i, row_id);
            }

            // 放 cell 到 R2, R3
            if (is_deferrable(row_id, u->s1, i) && is_placeable(row_id + 1, u->s2, i, u->l2) && is_placeable(row_id + 2, u->s3, i, u->l3)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1, u->l1 + 1, u->s2 + 1, 1, u->s3 + 1, 1, &Q);
                update_benefit(u, v, i, row_id);
            }

            // 放 cell 到 R1, R2, R3
            if (is_placeable(row_id, u->s1, i, u->l1) && is_placeable(row_id + 1, u->s2, i, u->l2) && is_placeable(row_id + 2, u->s3, i, u->l3)){
                MATRONode* v = create_or_get_node(row_id, i + 1, u->s1 + 1, 1, u->s2 + 1, 1, u->s3 + 1, 1, &Q);
                update_benefit(u, v, i, row_id);
            }
        }
        // cout << matro_node_table.size() << " nodes created." << endl;
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    // cout << "Time taken: " << duration.count() << " milliseconds" << endl;
    MATRONode* best_node = nullptr;
    int max_benefit = INT_MIN;
    for (auto& [key, node] : matro_node_table) {
        if (
            (row1 >= number_of_rows || node->s1 == row_cells[row1].size() - 1) &&
            (row2 >= number_of_rows || node->s2 == row_cells[row2].size() - 1) &&
            (row3 >= number_of_rows || node->s3 == row_cells[row3].size() - 1) 
        ) {
            for (int k = 0; k < 5; ++k) {
                if (node->benefit_case[k] > max_benefit) {
                    max_benefit = node->benefit_case[k];
                    best_node = node;
                }
            }
        }
    }
    // cout << matro_node_table.size() << " nodes created." << endl;
    // cout << "Best staple count = " << max_benefit << endl;

    // backtrack


    // 找 best case index
    int best_case = 0;
    for (int k = 0; k < 5; ++k) {
        if (best_node->benefit_case[k] == max_benefit) {
            best_case = k;
            break;
        }
    }

    best_path.push_back(best_node);
    best_case_used.push_back(best_case);

    while (best_node->prev_case[best_case] != nullptr) {
        MATRONode* prev = best_node->prev_case[best_case];

        int new_case = best_node->prev_case_index[best_case];
        // for (int k = 0; k < 5; ++k) {
        //     int delta = (best_case == 0) ? 0 : (best_case == 4 ? 2 : 1);
        //     if (prev->benefit_case[k] + delta == best_node->benefit_case[best_case]) {
        //         new_case = k;
        //         break;
        //     }
        // }

        best_node = prev;
        best_case = new_case;
        best_path.push_back(best_node);
        best_case_used.push_back(best_case);
    }
    // cout << "Best path size: " << best_path.size() << endl;
    // cout << "Best case used size: " << best_case_used.size() << endl;
}

void draw_svg_all() {

    // for (int g = 0; g < all_paths.size(); ++g) {
    //     const auto& path = all_paths[g];
    //     for (int p = path.size() - 1; p >= 0; --p) {
    //         MATRONode* node = path[p];
    //         int row_id = g * 3;
    //         auto update_cell_pos = [&](int row_offset, int site_index, int s_idx) {
    //             if (s_idx >= 0 && s_idx < row_cells[row_id + row_offset].size()) {
    //                 Cell* cell = row_cells[row_id + row_offset][s_idx];
    //                 int new_x = chip_left_boundary + site_index * site_width;
    //                 cell->x = new_x;
    //                 cell->y = chip_bottom_boundary + (row_id + row_offset) * row_height;
    //             }
    //         };

    //         update_cell_pos(0, node->i - node->l1, node->s1);
    //         update_cell_pos(1, node->i - node->l2, node->s2);
    //         update_cell_pos(2, node->i - node->l3, node->s3);
    //     }
    // }

    double scale = 0.1; 
    double svg_width = (chip_right_boundary - chip_left_boundary);
    double svg_height = (chip_top_boundary - chip_bottom_boundary);
    ofstream svg("layout.svg");
    svg << "<svg xmlns='http://www.w3.org/2000/svg' width='"
        << svg_width * scale << "' height='" << svg_height * scale << "'>\n";

    svg << "<g stroke='lightgray' stroke-width='0.5'>\n";
    for (int r = 0; r <= number_of_rows; ++r) {
        int y = chip_bottom_boundary + r * row_height;
        svg << "<line x1='" << chip_left_boundary * scale
            << "' y1='" << (chip_top_boundary - y) * scale
            << "' x2='" << chip_right_boundary * scale
            << "' y2='" << (chip_top_boundary - y) * scale
            << "' stroke='black' stroke-width='1'/>\n";
    }

    for (int s = 0; s <= number_of_sites_per_row; ++s) {
        int x = chip_left_boundary + s * site_width;
        svg << "<line x1='" << x * scale
            << "' y1='0' x2='" << x * scale
            << "' y2='" << (chip_top_boundary - chip_bottom_boundary) * scale
            << "' stroke='red' stroke-width='1'/>\n";
    }

    // 畫 cell
    for (auto& row : row_cells) {
        for (Cell* cell : row) {
            CellType* type = cell_type_map[cell->type_id];
            svg << "<rect x='" << cell->x * scale << "' y='" << (svg_height - cell->y - type->height) * scale
                << "' width='" << type->width * scale << "' height='" << type->height * scale
                << "' fill='rgba(173,216,230,0.8)' stroke='black'/>\n";

            int row_id = (cell->y - chip_bottom_boundary) / row_height;
            int start_site = (cell->x - chip_left_boundary) / site_width;

            for (int pin_site : type->pins) {
                int site_idx = start_site + pin_site;
                int pin_x = chip_left_boundary + site_idx * site_width + site_width / 4;
                int pin_y = chip_top_boundary - (row_id * row_height + 3 * row_height / 4);
                int pin_w = site_width / 2;
                int pin_h = row_height / 2;

                svg << "<rect x='" << pin_x * scale
                    << "' y='" << pin_y * scale
                    << "' width='" << pin_w * scale
                    << "' height='" << pin_h * scale
                    << "' fill='red'/>\n";
            }
            svg << "<text x='" << cell->x * scale
                << "' y='" << (svg_height - cell->y - type->height) * scale
                << "' font-size='6' fill='black'>" 
                << "id-" << cell->id << "-id" << "</text>\n";
        }
    }

    // 畫每個 site index
    // for (int r = 0; r < number_of_rows; ++r) {
    //     for (int s = 0; s < number_of_sites_per_row; ++s) {
    //         int x = chip_left_boundary + s * site_width;
    //         int y = chip_top_boundary - (r * row_height + row_height / 2);
    //         svg << "<text x='" << x * scale
    //             << "' y='" << y * scale
    //             << "' font-size='6' fill='black'>"
    //             << '(' << r << ", " << s << ')' << "</text>\n";
    //     }
    // }

    // staple
    for (int g = 0; g < all_paths.size(); ++g) {
        const auto& path = all_paths[g];
        const auto& cases = all_cases[g];
        int row_id = g * 3;

        for (int idx = path.size() - 1; idx >= 0; --idx) {
            MATRONode* node = path[idx];
            int case_id = cases[idx];
            int site = node->i - 1;

            if (site < 0) continue;

            int x = chip_left_boundary + site * site_width;
            auto draw_staple = [&](int row) {
                if (row < 0 || row >= number_of_rows) return;
                int y = chip_bottom_boundary + row * row_height;
                svg << "<rect x='" << x * scale
                    << "' y='" << (svg_height - y - row_height) * scale
                    << "' width='" << site_width * scale
                    << "' height='" << row_height * scale
                    << "' fill='rgba(0,255,0,0.3)' stroke='black'/>\n";
            };

            switch (case_id) {
                case 1: draw_staple(row_id); draw_staple(row_id - 1); break;
                case 2: draw_staple(row_id); draw_staple(row_id + 1); break;
                case 3: draw_staple(row_id + 1); draw_staple(row_id + 2); break;
                case 4: draw_staple(row_id); draw_staple(row_id - 1); draw_staple(row_id + 1); draw_staple(row_id + 2); break;
                default: break;
            }
        }
    }

    svg << "</g>\n";
    svg << "</svg>\n";
    svg.close();
    cout << "SVG file generated: layout.svg" << endl;
}

pair<int, int> cal_balance() {
    int staple_vdd = 0, staple_gnd = 0;
    int staple_count = 0;
    for (int g = 0; g < all_paths.size(); ++g) {
        const auto& path = all_paths[g];
        const auto& cases = all_cases[g];
        int row_id = g * 3;

        for (int idx = path.size() - 1; idx >= 0; --idx) {
            MATRONode* node = path[idx];
            int case_id = cases[idx];
            int site = node->i - 1;

            if (site < 0) continue;

            int x = chip_left_boundary + site * site_width;
            auto draw_staple = [&](int row) {
                if((row_id & 1) == 0){
                    if(row == row_id) 
                        staple_vdd++;
                    else if(row == row_id - 1 || row == row_id + 1) 
                        staple_gnd++;
                }else{
                    if(row == row_id) 
                        staple_gnd++;
                    else if(row == row_id - 1 || row == row_id + 1) 
                        staple_vdd++;
                }
                staple_count++;
            };

            switch (case_id) {
                case 1: draw_staple(row_id - 1); break;
                case 2: draw_staple(row_id);break;
                case 3: draw_staple(row_id + 1); break;
                case 4: draw_staple(row_id - 1); draw_staple(row_id + 1);  break;
                default: break;
            }
        }
    }
    // cout << "VDD staples: " << staple_vdd << endl;
    // cout << "GND staples: " << staple_gnd << endl;
    double ratio = (double)max(staple_vdd, staple_gnd) / (double)min(staple_vdd, staple_gnd);
    cout << "Balance ratio: " << ratio << endl;
    int removed_vdd = 0, removed_gnd = 0;
    if(ratio > 1.1){
        if(staple_vdd > staple_gnd){
            removed_vdd = ceil((double)staple_vdd - 1.1 * staple_gnd);
        }else{
            removed_gnd = ceil((double)staple_gnd - 1.1 * staple_vdd);
        }
    }
    return {removed_vdd, removed_gnd};
}

void write_out(ofstream &output_file) {
    for (auto& row : row_cells) {
        for (Cell* cell : row) {
            output_file << cell->id << " " << cell->x << " " << cell->y << " " << 0 << endl;
        }
    }
    pair<int, int> balance = cal_balance();
    int removed_vdd = balance.first;
    int removed_gnd = balance.second;
    // cout << "Removed VDD staples: " << removed_vdd << endl;
    // cout << "Removed GND staples: " << removed_gnd << endl;
    int staple_count = 0;
    for (int g = 0; g < all_paths.size(); ++g) {
        const auto& path = all_paths[g];
        const auto& cases = all_cases[g];
        int row_id = g * 3;

        for (int idx = path.size() - 1; idx >= 0; --idx) {
            MATRONode* node = path[idx];
            int case_id = cases[idx];
            int site = node->i - 1;

            if (site < 0) continue;

            int x = chip_left_boundary + site * site_width;
            auto draw_staple = [&](int row) {
                if((row_id & 1) == 0){
                    if(removed_vdd > 0 && row == row_id){ 
                        removed_vdd--;
                        return;
                    }else if(removed_gnd > 0 && (row == row_id - 1 || row == row_id + 1)){
                        removed_gnd--;
                        return;
                    }
                }else{
                    if(removed_gnd > 0 && row == row_id){
                        removed_gnd--;
                        return;
                    } else if(removed_vdd > 0 && (row == row_id - 1 || row == row_id + 1)){
                        removed_vdd--;
                        return;
                    }
                }

                int y = chip_bottom_boundary + row * row_height;
                output_file << x << " " << y << endl;
                staple_count++;
            };

            switch (case_id) {
                case 1: draw_staple(row_id - 1); break;
                case 2: draw_staple(row_id);break;
                case 3: draw_staple(row_id + 1); break;
                case 4: draw_staple(row_id - 1); draw_staple(row_id + 1);  break;
                default: break;
            }
        }
    }
    cout << "Total staples: " << staple_count << endl;
}


int main(int argc, char *argv[]) {
    // Open the input file
    ifstream input_file(argv[1]);
    if (!input_file.is_open()) {
        cout << "Error: could not open file " << argv[1] << endl;
        return 1;
    }

    ofstream output_file(argv[2]);
    if (!output_file.is_open()) {
        cout << "Error: could not open file " << argv[2] << endl;
        return 1;
    }

    vector<string> input_objects = input_preprocess(input_file);
    int object_size = input_objects.size();
    int input_indx = 0;

    vector<string> chip_boundaries = split(input_objects[input_indx++]);
    chip_left_boundary = stoi(chip_boundaries[0]);
    chip_bottom_boundary = stoi(chip_boundaries[1]);
    chip_right_boundary = stoi(chip_boundaries[2]);
    chip_top_boundary = stoi(chip_boundaries[3]);

    vector<string> row_info = split(input_objects[input_indx++]);
    number_of_rows = stoi(row_info[0]);
    row_height = stoi(row_info[1]);
    site_width = stoi(row_info[2]);
    number_of_sites_per_row = (chip_right_boundary - chip_left_boundary) / site_width;

    vector<string> cell_type_info = split(input_objects[input_indx++]);
    number_of_cell_types = stoi(cell_type_info[0]);

    vector<string> cell_info = split(input_objects[input_indx++]);
    number_of_cells = stoi(cell_info[0]);

    for (int index = 0; index < number_of_cell_types; index++) {
        vector<string> cell_type = split(input_objects[input_indx++]);
        int id = stoi(cell_type[0]);
        int width = stoi(cell_type[1]);
        int height = stoi(cell_type[2]);
        CellType* cell_type_obj = new CellType(id, width, height);
        for (int i = 3; i < cell_type.size(); i++) {
            int pin_id = stoi(cell_type[i]);
            cell_type_obj->pins.insert(pin_id);
        }
        cell_type_map[id] = cell_type_obj;
    }

    for (int index = 0; index < number_of_cells; index++) {
        vector<string> cell = split(input_objects[input_indx++]);
        int id = stoi(cell[0]);
        int type_id = stoi(cell[1]);
        int x = stoi(cell[2]);
        int y = stoi(cell[3]);
        int max_displacement = stoi(cell[4]);

        double ratio;
        if (number_of_cells < 10000) ratio = 1.0;
        else if (number_of_cells < 60000) ratio = 0.9;
        else if (number_of_cells < 100000) ratio = 0.85;
        else if (number_of_cells < 300000) ratio = 0.6;
        else ratio = 0.5;

        Cell* cell_obj = new Cell(id, type_id, x, y, max_displacement * ratio);
        cell_map[id] = cell_obj;
    }

    row_cells.resize(number_of_rows);
    original_site_cell_map.resize(number_of_rows, vector<int>(number_of_sites_per_row, -1));
    for (auto& pair  : cell_map) {
        Cell* cell = pair.second;
        int row_id = (cell->y - chip_bottom_boundary) / row_height;
        row_cells[row_id].push_back(cell);
        
        // for every site this cell occupies
        int start_site = (cell->x - chip_left_boundary) / site_width;
        int site_span = cell_type_map[cell->type_id]->width / site_width;
        for (int i = 0; i < site_span; ++i) {
            original_site_cell_map[row_id][start_site + i] = cell->id;
        }
    }
    // Sort each row's cells by their x coordinate
    for (int r = 0; r < number_of_rows; ++r) {
        sort(row_cells[r].begin(), row_cells[r].end(), [](Cell* a, Cell* b) {
            return a->x < b->x;
        });
    }

    for (int r = 0; r < number_of_rows; ++r) {
        for (int s = 0; s < number_of_sites_per_row; ++s) {
            if (original_site_cell_map[r][s] == -1 && s > 0) {
                original_site_cell_map[r][s] = original_site_cell_map[r][s - 1];
            }
        }
    }

    // occupied.resize(number_of_rows, vector<bool>(number_of_sites_per_row, false));
    // has_pin.resize(number_of_rows, vector<bool>(number_of_sites_per_row, false));

    // for (int r = 0; r < number_of_rows; ++r) {
    //     for (Cell* cell : row_cells[r]) {
    //         int start_site = (cell->x - chip_left_boundary) / site_width;
    //         int site_span = cell_type_map[cell->type_id]->width / site_width;
    //         for (int i = 0; i < site_span; ++i) {
    //             occupied[r][start_site + i] = true;
    //         }
    //         for (int pin_site : cell_type_map[cell->type_id]->pins) {
    //             has_pin[r][start_site + pin_site] = true;
    //         }
    //     }
    // }
    // cout << "number_of_rows = " << number_of_rows << endl;
    for (int i = 0; i < number_of_rows; i += 3) {
        triple_row_dp(i);
        all_paths.push_back(best_path);
        all_cases.push_back(best_case_used);

        const auto& path = all_paths.back();
        for (int p = path.size() - 1; p >= 0; --p) {
            MATRONode* node = path[p];
            int row_id = i;
            auto update_cell_pos = [&](int row_offset, int site_index, int s_idx) {
                if (s_idx >= 0 && s_idx < row_cells[row_id + row_offset].size()) {
                    Cell* cell = row_cells[row_id + row_offset][s_idx];
                    int new_x = chip_left_boundary + site_index * site_width;
                    cell->x = new_x;
                    // cell->y = chip_bottom_boundary + (row_id + row_offset) * row_height;
                }
            };

            update_cell_pos(0, node->i - node->l1, node->s1);
            update_cell_pos(1, node->i - node->l2, node->s2);
            update_cell_pos(2, node->i - node->l3, node->s3);
        }

        reset();
    }
    write_out(output_file);
    // draw_svg_all();
    input_file.close();
    output_file.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    cout << "Execution time: "
         << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count()
         << " s" << endl;
    return 0;
}