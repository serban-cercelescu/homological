#include <cassert>

#include <algorithm>
#include <execution>
#include <iostream>
#include <map>
#include <set>
#include <thread>
#include <utility>
#include <vector>
using namespace std;

using pii = pair<int, int>;
using i64 = long long;


class Sparse {
    int n, m;
    map<int, map<int, i64>> mx;

public:
    i64 get(int i, int j) {
        auto rowit = mx.find(i);
        if (rowit == mx.end()) {
            return 0;
        }
        else {
            auto colit = rowit->second.find(j);
            if (colit == rowit->second.end()) {
                return 0;
            }
            else {
                return colit->second;
            }
        }
    }

    void put(int i, int j, i64 val) {
        if (val == 0) {
            auto rowit = mx.find(i);
            if (rowit != mx.end()) {
                auto &row = rowit->second;
                row.erase(j);
                if (row.empty())
                    mx.erase(i);
            }
        }
        else {
            mx[i][j] = val;
        }
    }

    pii size() {
        return pii(n, m);
    }

    map<int, i64>& operator [] (int i) {
        return mx[i];
    }

    void swap_row(int src, int trg) {
        swap(mx[src], mx[trg]);
        if (mx[src].empty())
            mx.erase(src);
        if (mx[trg].empty())
            mx.erase(trg);
    }

    set<int> column_entries(int col) {
        set<int> rows;
        for (auto &row: mx) {
            int row_idx = row.first;
            auto &row_content = row.second;

            if (row_content.find(col) != row_content.end())
                rows.insert(row_idx);
        }
        return rows;
    }

    set<int> row_entries(int row) {
        set<int> cols;
        auto it = mx.find(row);
        if (it == mx.end()) {
            return cols;
        }
        else {
            for (auto &col: it->second)
                cols.insert(col.first);
            return cols;
        }
    }

    void mul_row(int src, i64 qt) {
        auto it = mx.find(src);
        if (it == mx.end()) {
            return;
        }
        else {
            vector<pair<int, i64>> row = get_row(src);
            for_each(
                execution::unseq,
                begin(row),
                end(row),
                [&](const pair<int, long long>& entry) {
                    put(src, entry.first, entry.second * qt);
                }
            );
        }
    }

    void mul_col(int src, i64 qt) {
        set<int> rows = column_entries(src);

        for_each(
            execution::unseq,
            begin(rows),
            end(rows),
            [&](const int &row) {
                auto itrow = mx.find(row);
                if (itrow != mx.end()) {
                    auto itcol = itrow->second.find(src);
                    if (itcol != itrow->second.end()) {
                        itcol->second *= qt;
                    }
                }
            }
        );
    }

    void swap_col(int src, int trg) {
        set<int> rows = column_entries(src);
        rows.merge(column_entries(trg));

        for_each(
            execution::unseq,
            begin(rows),
            end(rows),
            [&](const int &row) {
                int a = get(row, src);
                int b = get(row, trg);
                put(row, src, b);
                put(row, trg, a);
            }
        );
    }

    void add_row(int src, int trg, i64 qt) {
        auto it = mx.find(src);
        if (it != mx.end()) {
            for_each(
                execution::unseq,
                begin(it->second),
                end(it->second),
                [&](const auto& row_entry) {
                    int col_pos = row_entry.first;
                    i64 entry_val = row_entry.second;
                    i64 target_val = get(trg, col_pos);
                    put(trg, col_pos, target_val + entry_val * qt);
                }
            );
        }
    }

    void add_col(int src, int trg, i64 qt) {
        set<int> rows = column_entries(src);
        rows.merge(column_entries(trg));

        for_each(
            execution::unseq,
            begin(rows),
            end(rows),
            [&](const int &row) {
                i64 a = get(row, src);
                i64 b = get(row, trg);
                put(row, trg, b + a * qt);
            }
        );
    }

    vector<pair<int, i64>> get_row(int i) {
        vector<pair<int, i64>> row;
        if (mx.find(i) == mx.end())
            return row;
        else {
            for (auto &entry: mx[i])
                row.push_back(entry);
        }
        return row;
    }

    vector<pair<int, i64>> get_col(int j) {
        vector<pair<int, i64>> col;
        for (auto &row: mx)
            if (row.second.find(j) != row.second.end())
                col.push_back(make_pair(row.first, row.second[j]));
        return col;
    }

    void set_row(int i, vector<pair<int, i64>> row) {
        mx[i].clear();
        mx[i] = map<int, i64>(row.begin(), row.end());
    }

    void set_column(int j, vector<pair<int, i64>> col) {
        for_each(
            execution::unseq,
            begin(col),
            end(col),
            [&](const pair<int, i64>& entry) {
                put(entry.first, j, entry.second);
            }
        );
        
        for_each(
            execution::unseq,
            begin(col),
            end(col),
            [&](const pair<int, i64>& entry) {
                put(entry.first, j, entry.second);
            }
        );
    }

    Sparse(int _n, int _m) {
        n = _n;
        m = _m;

        if (n >= (1 << 21) || m >= (1 << 21))
            throw runtime_error("Sparse matrix too large, size limit is 2097151");
    }
};

class Matrix {
    int n, m;
    vector<vector<i64>> mx;

public:
    pii size() {
        return pii(n, m);
    }

    void swap_row(int src, int trg) {
        swap(mx[src], mx[trg]);
    }

    void swap_col(int src, int trg) {
        for (int i = 0; i < n; ++i)
            swap(mx[i][src], mx[i][trg]);
    }
    
    void mul_col(int src, i64 qt) {
        for (int i = 0; i < n; ++i)
            mx[i][src]*= qt;
    }
    
    void mul_row(int src, i64 qt) {
        for (int j = 0; j < m; ++j)
            mx[src][j]*= qt;
    }

    void add_row(int src, int trg, i64 qt) {
        for (int j = 0; j < m; ++j)
            mx[trg][j]+= mx[src][j] * qt;
    }

    void add_col(int src, int trg, i64 qt) {
        for (int i = 0; i < n; ++i)
            mx[i][trg]+= mx[i][src] * qt;
    }

    vector<i64> &operator [] (int i) {
        return mx[i];
    }    

    Matrix(int _n, int _m) {
        n = _n;
        m = _m;
        mx = vector<vector<i64>>(n, vector<i64>(m, 0));
    }
};

Matrix operator * (Matrix &a, Matrix &b) {
    int n, m, p;
    n = a.size().first;
    m = a.size().second;
    p = b.size().second;

    Matrix c(n, p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            for (int k = 0; k < m; ++k)
                c[i][j]+= a[i][k] * b[k][j];
    return c;
}

Matrix identity(int n) {
    Matrix matrix(n, n);
    for (int i = 0; i < n; ++i)
        matrix[i][i] = 1;
    return matrix;
}


// returns a matrix such that matrix * [a, b] = [gcd(a, b), 0]
Matrix extended_euclid(i64 a, i64 b) { // I'm very skeptical that this will work
    if (b == 0)
        return identity(2);
    Matrix rec = extended_euclid(b, a % b);
    Matrix matrix(2, 2);
    matrix[0][0] = 0;
    matrix[0][1] = 1;
    matrix[1][0] = 1;
    matrix[1][1] = -a / b;
    return rec * matrix;
}

bool place_nonzero(Sparse &matrix, int t) {
    int n, m;
    tie(n, m) = matrix.size();
    
    for (int i = t; i <= n; ++i) {
        for (int j = t; j <= m; ++j) {
            if (matrix.get(i, j) != 0) {
                matrix.swap_row(t, i);
                matrix.swap_col(t, j);
                return true;
            }
        }
    }
    return false;
}

void optimize_pivot(Sparse &matrix, int t) { // ensures that matrix[t][t] is the gcd of its row and column
    int n = matrix.size().first;
    int m = matrix.size().second;    
    
    if (matrix.get(t, t) < 0) {
        matrix.mul_row(t, -1);
    }

    set<int> rows = matrix.column_entries(t);
    for (int i: rows) if (i > t) {
        i64 a = matrix.get(t, t);
        i64 b = matrix.get(i, t);

        if (b % a == 0)
            continue;
        assert(a >= 0);
        if (b < 0) {
            b = -b;
            matrix.mul_row(i, -1);
        }

        Matrix euclid_matrix = extended_euclid(a, b); // euclid_matrix * [a, b] = [gcd(a, b), 0]
        assert(euclid_matrix[0][0] * a + euclid_matrix[0][1] * b == __gcd(a, b));
        assert(euclid_matrix[1][0] * a + euclid_matrix[1][1] * b == 0);

        auto row_i = matrix.get_row(i);
        auto row_t = matrix.get_row(t);

        matrix.set_row(t, {});
        matrix.set_row(i, {});

        for (auto &entry: row_t) {
            int j = entry.first;
            i64 val_t = entry.second;
            i64 val_i = [&]() -> i64 {
                int lb = -1;
                for (int msk = 1 << 20; msk > 0; msk/= 2)
                    if (lb + msk < row_i.size() && row_i[lb + msk].first <= j)
                        lb+= msk;
                if (lb == -1 || row_i[lb].first != j)
                    return 0;
                else
                    return row_i[lb].second;
            }();

            matrix.put(t, j, euclid_matrix[0][0] * val_t + euclid_matrix[0][1] * val_i);
            matrix.put(i, j, euclid_matrix[1][0] * val_t + euclid_matrix[1][1] * val_i);
        }

        for (auto &entry: row_i) {
            int j = entry.first;
            i64 val_i = entry.second;
            i64 val_t = [&]() -> i64 {
                int lb = -1;
                for (int msk = 1 << 20; msk > 0; msk/= 2)
                    if (lb + msk < row_t.size() && row_t[lb + msk].first <= j)
                        lb+= msk;
                if (lb == -1 || row_t[lb].first != j)
                    return 0;
                else
                    return row_t[lb].second;
            }();

            matrix.put(t, j, euclid_matrix[0][0] * val_t + euclid_matrix[0][1] * val_i);
            matrix.put(i, j, euclid_matrix[1][0] * val_t + euclid_matrix[1][1] * val_i);
        }
    }

    // we make all other entries in the row zero
    set<int> cols = matrix.row_entries(t);
    for (int j: cols) if (j > t) {
        i64 a = matrix.get(t, t); 
        i64 b = matrix.get(t, j);

        if (b % a == 0)
            continue;
        assert(a >= 0);
        if (b < 0) {
            b = -b;
            matrix.mul_col(j, -1);
        }

        Matrix euclid_matrix = extended_euclid(a, b); // euclid_matrix * [a, b] = [gcd(a, b), 0]
        assert(euclid_matrix[0][0] * a + euclid_matrix[0][1] * b == __gcd(a, b));
        assert(euclid_matrix[1][0] * a + euclid_matrix[1][1] * b == 0);

        auto col_t = matrix.get_col(t);
        auto col_j = matrix.get_col(j);

        matrix.set_column(t, {});
        matrix.set_column(j, {});

        for (auto &entry: col_t) {
            int i = entry.first;
            i64 val_t = entry.second;
            i64 val_j = [&]() -> i64 {
                int lb = -1;
                for (int msk = 1 << 20; msk > 0; msk/= 2)
                    if (lb + msk < col_j.size() && col_j[lb + msk].first <= i)
                        lb+= msk;
                if (lb == -1 || col_j[lb].first != i)
                    return 0;
                else
                    return col_j[lb].second;
            }();

            matrix.put(i, t, euclid_matrix[0][0] * val_t + euclid_matrix[0][1] * val_j);
            matrix.put(i, j, euclid_matrix[1][0] * val_t + euclid_matrix[1][1] * val_j);
        }

        for (auto &entry: col_j) {
            int i = entry.first;
            i64 val_j = entry.second;
            i64 val_t = [&]() -> i64 {
                int lb = -1;
                for (int msk = 1 << 20; msk > 0; msk/= 2)
                    if (lb + msk < col_t.size() && col_t[lb + msk].first <= i)
                        lb+= msk;
                if (lb == -1 || col_t[lb].first != i)
                    return 0;
                else
                    return col_t[lb].second;
            }();

            matrix.put(i, t, euclid_matrix[0][0] * val_t + euclid_matrix[0][1] * val_j);
            matrix.put(i, j, euclid_matrix[1][0] * val_t + euclid_matrix[1][1] * val_j);
        }
    }

}

void clear_columns(Sparse &matrix, int t) {
    set<int> columns = matrix.row_entries(t);
    i64 a = matrix.get(t, t);
    for (auto j: columns) if (j > t) {
        i64 b = matrix.get(t, j);
        if (b == 0)
            continue;
        matrix.add_col(t, j, -b / a);
    }
}

void clear_rows(Sparse &matrix, int t) {
    set<int> rows = matrix.column_entries(t);
    i64 a = matrix.get(t, t);
    for (auto i: rows) if (i > t) {
        i64 b = matrix.get(i, t);
        if (b == 0)
            continue;
        matrix.add_row(t, i, -b / a);
    }
}

vector<int> smith_normal_form(Sparse matrix) {
    int n, m;
    tie(n, m) = matrix.size();

    for (int t = 1; t <= min(n, m); ++t) {
        if (!place_nonzero(matrix, t))
            break;
        // matrix[t][t] will now be the pivot
        optimize_pivot(matrix, t);
        clear_columns(matrix, t);
        clear_rows(matrix, t);

    }

    bool sorted = false;
    vector<int> diagf;
    diagf.reserve(min(n, m));
    for (int i = 1; i <= min(n, m); ++i)
        diagf.push_back(matrix[i][i]);

    for (int i = 0; i < diagf.size(); ++i) {
        for (int j = i + 1; j < diagf.size(); ++j) {
            if (diagf[j] > diagf[i]) {
                swap(diagf[i], diagf[j]);
            }
        }
        if (diagf[i] == 0) {
            sort(diagf.begin(), diagf.begin() + i + 1);
            sorted = true;
            break;
        }

        for (int j = i + 1; j < diagf.size(); ++j) if (diagf[j] != 0) {
            int a = diagf[i];
            int b = diagf[j];
            while (a % b != 0) {
                a%= b;
                swap(a, b);
            }
        }
    }

    if (!sorted)
        sort(begin(diagf), end(diagf));

    return diagf;
}

extern "C" {
    i64 *csmith_normal_form(i64[]);

    i64 *csmith_normal_form(i64 _mx[]) {
        int n = _mx[0], m = _mx[1];
        Sparse mx(n, m);
        
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= m; ++j)
                mx.put(i, j, _mx[(i - 1) * m + j - 1 + 2]);
        
        vector<int> snf = smith_normal_form(mx);
        i64 *ans = new i64[snf.size() + 1];

        ans[0] = snf.size();
        for (int i = 0; i < snf.size(); ++i) {
            ans[i + 1] = snf[i];
        }

        return ans;
    }
}
