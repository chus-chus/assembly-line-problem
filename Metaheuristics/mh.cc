/*
   Algorithmics and Programming III (Data Science & Engineering).
   Technical University of Catalonia.

   Project: Metaheuristics approach (Simulated Annealing).

                Aguilera González, Cristina
                            &
                Antoñanzas Acero, Jesús M.                      */

#include <algorithm>
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

/* Returns elapsed time since program execution up to the function call
(in seconds). */
double get_time () {
    auto timeend = chrono::steady_clock::now();
    auto time = timeend - timestart;
    return chrono::duration <double, milli> (time).count()/1000;
}

/* Writes current solution found in "result" on the previously specified file. */
void record_sol () {
    ofstream out (output_file);
    out << result.penalties << " " << setprecision (1) << fixed << result.time;
    out << endl;
    int n = result.chain.size();
    for (int i = 0; i < n; ++i) {
        if (i == 0) out << result.chain[i];
        else out << " " << result.chain[i];
    }
    out.close();
}

/* Records found solution. */
void new_sol (int npenal, const vector<int>& chain) {
    result.penalties = npenal;
    result.time = get_time();
    result.chain = chain;
    record_sol();
}

/* Computes penalties of a full chain of cars. For each upgrade station and for
every car, count the penalties associated with it in its corresponding window.
Notice that we start and finish on non-existing cars in order to tackle
the so called "incomplete" penalties. */
int whole_penalties (const vector<int>& chain, const vector<carclass>& classes,
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

/* Adjust temperature following a geometric law. */
double temperature (double T) {
    return 0.9999*T;
}

/* Returns a double between 0 and 1 (both included), depending on the
temperature "T" and the increment in penalties. */
double get_probability (double T, int neigh_penal, int curr_penal) {
    return exp(-(neigh_penal - curr_penal)/T);
}

/* Compute penalties associated to the cars in positions indexs[0] and
indexs[1] of the chain. Notice that for each car and upgrade there is more
than one interval covering it, so we have to take into account those as well
("incomplete" penalties are also computed if needed). */
int index_penalties (const vector<int>& chain, const vector<carclass>& classes,
                     const vector<pair<int,int>>& capacities,
                     const vector<int>& indexs) {
    int M = classes[0].upgrades.size();
    int n = chain.size();
    int npenal = 0;
    for (int i = 0; i < M; ++i) {
        int ne = capacities[i].second;
        int ce = capacities[i].first;
        for (int u = 0; u < 2; ++u) {
            int index = indexs[u];
            /* For each index, these are the intervals starting in "j" that
            cover it. */
            for (int j = index + 1 - ne; j <= index; ++j) {
                int counter = 0;
                /* Count upgrades of each interval of length "ne". */
                for (int y = j; y < j + ne; ++y) {
                    if (y >= 0 and y < n) {
                        if (classes[chain[y]].upgrades[i] == true) ++counter;
                    }
                }
                /* If we have more upgrades than the station can handle. */
                if (counter - ce > 0) npenal += counter - ce;
            }
        }
    }
    return npenal;
}

/* Randomly pick two different integers between 0 and n (both included). */
vector<int> pick_indexs (const vector<int>& chain) {
    int n = chain.size();
    int i = rand() % n;
    int j = rand() % n;
    while (i == j) i = rand() % n;
    return {i, j};
}

/* Returns a neighbour of chain "curr_chain" and computes its penalties,
writing them in "neigh_penal". We consider a neighbour a chain with two cars
in different positions. Instead of computing the penalties for the whole
neighbour chain, the procedure we follow is:

    1- Compute penalties of initial chain (those corresponding to the two
       cars we already know are going to swap).
    2- Swap cars.
    3- Compute penalties of neighbour chain, but only those associated with
       the cars we've just swapped (in their new position).
    4- Penalties of the new chain will be the penalties of the whole
       initial chain minus the ones computed in step 1 plus the ones
       computed in step 3.

The cost of this technique may be similar to calculating all penalties if
the chain is small, but for large chains there is a noticeable advantage. */
vector<int> compute_neigh (const vector<int>& curr_chain, int curr_penal,
                           const vector<carclass>& classes, int& neigh_penal,
                           const vector<pair<int,int>>& capacities) {
    vector<int> indexs = pick_indexs (curr_chain);
    vector<int> neighbour = curr_chain;

    int beforeswap = index_penalties (curr_chain, classes, capacities, indexs);
    swap (neighbour[indexs[0]], neighbour[indexs[1]]);
    int afterswap = index_penalties (neighbour, classes, capacities, indexs);

    neigh_penal = (curr_penal - beforeswap) + afterswap;
    return neighbour;
}

/* Performs Simulated Annealing. */
void sim_annealing (const vector<carclass>& classes,
                    const vector<pair<int,int>>& capacities) {
    /* Set the seed for random generation. */
    unsigned rseed = chrono::system_clock::now().time_since_epoch().count();
    srand(rseed);
    double T = 1000;

    /* Set initial chain as the one computed before "sim_annealing" was called. */
    vector<int> curr_chain = result.chain;
    int curr_penal = result.penalties;

    /* We hope the user will stop the program execution after some running time. */
    while (true) {
        int neigh_penal;
        vector<int> neighbour = compute_neigh (curr_chain, curr_penal, classes,
                                               neigh_penal, capacities);
        /* Found a better solution. Set it as the current one. */
        if (neigh_penal < curr_penal) {
            new_sol (neigh_penal, neighbour);
            curr_chain = neighbour;
            curr_penal = neigh_penal;
        }

        /* Even if the new solution is worse (or the same), set it as the
        current one with a certain probability in order to avoid getting stuck
        in local minima. */
        else {
            double prob = get_probability (T, neigh_penal, curr_penal);
            if ((double) rand() / (RAND_MAX) < prob) {
                curr_chain = neighbour;
                curr_penal = neigh_penal;
            }
        }
        T = temperature(T);
    }
}

/* Generate a chain containing all cars available. */
vector<int> generate_chain (const vector<carclass>& classes, int C) {
    int K = classes.size();
    int position = 0;
    vector<int> chain(C);
    for (int i = 0; i < K; ++i) {
        int n = classes[i].quantity;
        for (int j = 0; j < n; ++j) {
            chain[position] = i;
            ++position;
        }
    }
    return chain;
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

/* Important note: Execution shall to be stopped whenever one wants in order
to get the desired level of approximation. When done so, the solution will be
found in the output file. */
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

    /* Starting chain for simulated annealing. */
    vector<int> initial_chain = generate_chain (classes, C);
    vector<int> used(C);

    new_sol (whole_penalties (initial_chain, classes, capacities), initial_chain);
    sim_annealing (classes, capacities);
}
