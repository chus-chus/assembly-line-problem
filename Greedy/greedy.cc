/*
   Algorithmics and Programming III (Data Science & Engineering).
   Technical University of Catalonia.

   Project: Greedy approach.

                Aguilera González, Cristina
                            &
                Antoñanzas Acero, Jesús M.                      */

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

/* A struct representing each car class. */
struct carclass {
    int quantity;
    vector<bool> upgrades;
};

/* Save the current solution: its penalties, time to obtain it and
 the solution itself. */
struct final_chain {
    int penalties;
    double time;
    vector<int> chain;
};

const int MAX_INT = 999999;

/* Set it as global as it enables quick access. */
final_chain result {MAX_INT, 0.0, {}};

/* Start time on startup. */
auto timestart = chrono::steady_clock::now();

/* For reading and writing, we need to access these files. It is more
   more convenient if they are global variables. */
char* output_file;
char* input_file;

/* Writes current solution found in "result" on the previously specified file. */
void record_sol () {
    ofstream out(output_file);
    out << result.penalties << " " << setprecision(1) << fixed << result.time;
    out << endl;
    int n = result.chain.size();
    for (int i = 0; i < n; ++i) {
        if (i == 0) out << result.chain[i];
        else out << " " << result.chain[i];
    }
    out.close();
}

/* Returns elapsed time since program execution up to the function call
(in seconds). */
double get_time () {
    auto timeend = chrono::steady_clock::now();
    auto time = timeend - timestart;
    return chrono::duration <double, milli> (time).count()/1000;
}

/* Records found solution. */
void new_sol (int npenal, const vector<int>& chain) {
    result.penalties = npenal;
    result.time = get_time();
    result.chain = chain;
    record_sol ();
}

/* Computes penalties of a full chain of cars. For each upgrade station and for
every car, count the penalties associated with it in its corresponding window.
Notice that we start and finish on non-existing cars in order to tackle
the so called "incomplete" penalties. */
int penalties (const vector<int>& chain, const vector<carclass>& classes,
               const vector<pair<int,int>>& capacities) {
    int M = classes[0].upgrades.size();
    int n = chain.size();
    int npenal = 0;
    for (int i = 0; i < M; ++i) {
        int ne = capacities[i].second;
        int ce = capacities[i].first;
        for (int j = 2 - ne; j <= n - 2; ++j) {
            int counter = 0;
            for (int y = j; y < j + ne; ++y) {
                if (y >= 0 and y < n) {
                    if (classes[chain[y]].upgrades[i] == true) ++counter;
                }
            }
            if (counter - ce > 0) npenal += counter - ce;
        }
    }
    return npenal;
}

/* Decides which kind of car fits best in the chain after a car of class
'prev_car'.
We say a car is contextually better if its coefficient associated with the
previous car is smaller than for any other car.
The coefficient is:

    "correlation_coef with previous car" / ("# cars left of the class")^2

So, when we have lots of cars left of a certain class, we need to release them,
as it can later become a penalty bottleneck. Also, the bigger "correlation_coef"
is, the less likely that car is to be put in that position.
Pre: chain contains at least one car. */
int contextual_best (int prev_car, const vector<vector<double>>& correlation_coef,
                     const vector<carclass>& classes) {
    int best_car = 0;
    double best_coef = double (MAX_INT);
    int K = correlation_coef.size();
    for (int i = 0; i < K; ++i) {
        double coef = correlation_coef[prev_car][i]/pow(classes[i].quantity, 2);
        if (classes[i].quantity != 0 and coef < best_coef) {
            best_coef = coef;
            best_car = i;
        }
    }
    return best_car;
}

/* Returns a chain of C cars, the first one being of class "first_car".
"classes" is a copy parameter because it's what we use to keep track of how many
cars we have used for each iteration of the function "solve". */
vector<int> find_chain (int C, const vector<vector<double>>& correlation_coef,
                        vector<carclass> classes, int first_car) {
    vector<int> chain(C, 0);
    chain[0] = first_car;
    --classes[first_car].quantity;
    for (int i = 1; i < C; ++i) {
        int to_insert = contextual_best(chain[i-1], correlation_coef, classes);
        --classes[to_insert].quantity;
        chain[i] = to_insert;
    }
    return chain;
}

/* Returns the best chain it can compute after trying all types of cars
as starting cars. */
void solve (int C, const vector<carclass>& classes,
            const vector<pair<int,int>>& capacities,
            const vector<vector<double>>& correlation_coef) {
    int K = classes.size();
    for (int first_car = 0; first_car < K; ++first_car) {
        vector<int> chain;
        chain = find_chain (C, correlation_coef, classes, first_car);
        int penal = penalties (chain, classes, capacities);
        if (penal < result.penalties) new_sol (penal, chain);
    }
}

/* Returns a vector containing the relationship between each pair of car types
(the sum of the division 'n_e/c_e' for each upgrade in common). */
vector<vector<double>> find_coinc (const vector<carclass>& classes,
                                   const vector<pair<int,int>>& capacities) {
    int K = classes.size();
    int M = capacities.size();
    double coef;
    vector<vector<double>> correlation_coef(K, vector<double>(K, MAX_INT));

    /* Compare the i'th class with all [i, K] others. We can do this and still
    compute all coefficients because the relationship between classes
    (0, 1) is the same as (1, 0). */
    for (int i = 0; i < K; ++i) {
        for (int j = i; j < K; ++j) {
            coef = 0;
            for (int u = 0; u < M; ++u) {
                if (classes[i].upgrades[u] and classes[j].upgrades[u]) {
                    coef += capacities[u].second/double(capacities[u].first);
                }
            }
            correlation_coef[i][j] = coef;
            if (i != j) correlation_coef[j][i] = coef;
        }
    }
    return correlation_coef;
}

void read (vector<carclass>& classes, vector<pair<int,int>>& capacities) {
    ifstream input(input_file);

    int C, M, K;
    input >> C >> M >> K;

    /* Read capacities of each station from input file. A station can upgrade
    c_e out of n_e cars at a time without penalty.
                    First -> c_e | Second -> n_e */
    for (int i = 0; i < M; ++i) input >> capacities[i].first;
    for (int i = 0; i < M; ++i) input >> capacities[i].second;

    /* Read last portion of the input. For each upgrade kind, first read the
    index, then the quantity of cars of that kind and finally which of the M
    upgrades the class needs. */
    int index;

    for (int i = 0; i < K; ++i) {
        input >> index;
        input >> classes[index].quantity;
        for (int j = 0; j < M; ++j) {
            bool upgrade; input >> upgrade;
            classes[index].upgrades.push_back(upgrade);
        }
    }
    input.close();
}

int main (int argc, char** argv) {
    input_file = argv[1];
    output_file = argv[2];

    ifstream input(input_file);

    /* C = Cars, M = upgrades, K = classes. */
    int C, M, K;
    input >> C >> M >> K;
    input.close();

    vector<pair<int,int>> capacities(M);
    vector<carclass> classes(K);
    read (classes, capacities);

    /* "correlation_coef": in each position 'i' (between 0 and K - 1), contains
    the relationship between car class 'i' and 'j' in form of a coefficient. */
    vector<vector<double>> correlation_coef;
    correlation_coef = find_coinc (classes, capacities);

    solve (C, classes, capacities, correlation_coef);
}
