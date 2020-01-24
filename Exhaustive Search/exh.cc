/*
   Algorithmics and Programming III (Data Science & Engineering).
   Technical University of Catalonia.

   Project: Exhaustive Search approach.

                Aguilera González, Cristina
                            &
                Antoñanzas Acero, Jesús M.                      */

#include <chrono>
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

/* We'll update the solution while the program is running, so we're better
 off using a global variable. */
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

/* Records found solution. */
void new_sol (int npenal, const vector<int>& chain) {
    result.penalties = npenal;
    result.time = get_time ();
    result.chain = chain;
    record_sol ();
}

/* Given an integer "partialpen", a partial solution "chain" and an integer
"index", updates the chain's penalties. Before function call, "partialpen"
contained the chain's penalties up to index - 1. "Penalties" will compute
the chain's penalties contained within the window of corresponding
width [u, index] for each upgrade station. */
void penalties (const vector<int>& chain, const vector<carclass>& classes,
                const vector<pair<int,int>>& capacities, int index, int& partialpen) {
    int M = classes[0].upgrades.size();  /* # of upgrades. */
    int npenal = 0;                      /* Total penalties for each upgrade. */

    for (int i = 0; i < M; ++i) {        /* Each upgrade station. */
        /* Reminder: i'th station can handle "ce" out of "ne" cars at a time. */
        int ne = capacities[i].second;
        int ce = capacities[i].first;

        /* Counter saves how many upgrades we're doing in a window. */
        int counter = 0;

        /* Count upgrades in the window (index - ne, index] and "incomplete"
        penalties. "classes[chain[y]].upgrades[i]" tells if the y'th car from
        the chain (chain[y]) needs the i'th upgrade. */
        for (int y = index; y > index - ne; --y) {
            if (y >= 0 and classes[chain[y]].upgrades[i]) ++counter;
        }

        /* If we have more upgrades than the station can handle. */
        if (counter - ce > 0) npenal += counter - ce;

        /* If at the end of the chain, check for incomplete penalties. */
        if (index == chain.size() - 1) {
            for (int j = index + 2 - ne; j < index; ++j) {
                counter = 0;
                for (int y = j; y < j + ne; ++y) {
                    if (y >= 0 and y <= index) {
                        if (classes[chain[y]].upgrades[i]) ++counter;
                    }
                }
                if (counter - ce > 0) npenal += counter - ce;
            }
        }
    }
    partialpen += npenal;

    /* We have a better solution if the below condition checks. */
    if (partialpen < result.penalties and index == chain.size() - 1) {
        new_sol (partialpen, chain);
    }
}

/* Generate permutations of size C containing all available cars.
Each car appears "classes[i].quantity" times. */
void permutations (int C, int K, const vector<carclass>& classes,
                   vector<int>& chain, vector<int>& used, int partialpen,
                   const vector<pair<int,int>>& capacities, int index) {
    if (index == C) return;
    else {
        for (int i = 0; i < K; ++i) {
            if (used[i] < classes[i].quantity) {
                chain[index] = i;
                ++used[i];
                int temporal_pen = partialpen;
                penalties (chain, classes, capacities, index, partialpen);

                /* Only continue with the current chain if it has less penalties
                than the previous chain generated. */
                if (partialpen < result.penalties) {
                    permutations (C, K, classes, chain, used, partialpen,
                                  capacities, index + 1);
                }
                partialpen = temporal_pen;
                --used[i];
            }
        }
    }
}

void solve (int C, const vector<carclass>& classes,
            const vector<pair<int,int>>& capacities) {
    int K = classes.size();

    /* "partialpen" will contain penalties of a chain up to some index. */
    int partialpen = 0;
    vector<int> chain(C, 0);
    vector<int> used(K, 0);
    permutations (C, K, classes, chain, used, partialpen, capacities, 0);
}

void read (vector<carclass>& classes, vector<pair<int,int>>& capacities) {
    ifstream input (input_file);

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
    /* One should remember these are global variables. */
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

    solve (C, classes, capacities);
}
