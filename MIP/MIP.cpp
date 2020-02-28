#include <string>
#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>
#include "gurobi_c++.h"

using namespace std;

int main() {
    string line;
    getline(cin, line);
    while(line[0] == 'c')
        getline(cin, line);

    istringstream iss(line);

    int n, m;
    string d1, d2;

    iss >> d1 >> d2 >> n >> m;

    vector<vector<int>> adj = vector<vector<int>>(n, vector<int>());
    cout << "Graph:" << endl;
    for(int i = 0; i < n; i++) {
        cout << i << ": ";
        for(int j : adj[i])
            cout << j << " ";
        cout << endl;
    }

    for(int i = 0; i < m; i++) {
        getline(cin, line);
        while(line[0] == 'c') {
            cout << "Skipping comment." << endl;
            getline(cin, line);
        }
        istringstream iss(line);

        int a, b; 
        iss >> a >> b;
        a--; b--;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    cout << "Graph:" << endl;
    for(int i = 0; i < n; i++) {
        cout << i << ": ";
        for(int j : adj[i])
            cout << j << " ";
        cout << endl;
    }

    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "mip1.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        cout << "Adding root vars." << endl;

        auto is_root_vars = vector<GRBVar>(n);
        GRBLinExpr one_root = 0;
        for(int i = 0; i < n; i++) {
            is_root_vars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            one_root += is_root_vars[i];
        }

        model.addConstr(one_root == 1);

        cout << "Adding parent_of vars." << endl;

        auto parent_of_vars = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
        for(int i = 0; i < n; i++) {
            parent_of_vars[i][i] = model.addVar(0.0, 0.0, 0.0, GRB_BINARY);
            GRBLinExpr one_parent = is_root_vars[i];
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    parent_of_vars[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    one_parent += parent_of_vars[i][j];
                }
            }
            model.addConstr(one_parent == 1);
        }

        cout << "Adding ascendant_of vars." << endl;

        auto asc_of_vars = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
        for(int i = 0; i < n; i++) {
            asc_of_vars[i][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            model.addConstr(asc_of_vars[i][i] == 1);

            for(int j = 0; j < n; j++) {
                if(i != j) {
                    asc_of_vars[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                }
            }
        }

        // Ascendant constraints with quadratic constraints:
        /*for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i == j)
                    continue;
                GRBQuadExpr asc_i_j = 0;
                for(int k = 0; k < n; k++) {
                    if(k == i)
                        continue;
                    asc_i_j += parent_of_vars[i][k] * asc_of_vars[k][j];
                }

                model.addQConstr(asc_i_j == asc_of_vars[i][j]);
            }
        }*/

        // Ascendant constraints with binary variables:
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i == j)
                    continue;
                GRBLinExpr asc_i_j = 0;
                

                for(int k = 0; k < n; k++) {
                    if(k == i)
                        continue;
                    auto bit = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    auto aux = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                    model.addConstr(aux <= parent_of_vars[i][k]);
                    model.addConstr(aux <= asc_of_vars[k][j]);
                    model.addConstr(aux >= parent_of_vars[i][k] - bit);
                    model.addConstr(aux >= asc_of_vars[k][j] + bit - 1);

                    asc_i_j += aux;
                }
                model.addConstr(asc_i_j == asc_of_vars[i][j]);
            }
        }

        cout << "Adding edge constraints." << endl;
        for(int i = 0; i < n; i++) {
            for(int j : adj[i]) {
                model.addConstr(asc_of_vars[i][j] + asc_of_vars[j][i] >= 1);
            }
        }

        cout << "Adding depth vars." << endl;

        auto depth_vars = vector<GRBVar>(n);
        for(int i = 0; i < n; i++) {
            depth_vars[i] = model.addVar(0.0, n - 1, 0.0, GRB_CONTINUOUS);
        }

        cout << "Adding depth var constraints." << endl;

        // With quadratic constraints.
        /* for(int i = 0; i < n; i++) {
            GRBQuadExpr depth_of_i = -n*is_root_vars[i];
            for(int j = 0; j < n; j++) {
                if(i == j)
                    continue;

                depth_of_i += parent_of_vars[i][j] * (1 + depth_vars[j]);
            }
            model.addQConstr(depth_vars[i] >= depth_of_i);
        }*/
        // With linear constraints.
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i == j)
                    continue;
                model.addConstr(depth_vars[i] >= 1 + depth_vars[j]
                                                 - n*(1 - parent_of_vars[i][j]));
            }
        }

        cout << "Adding full depth var." << endl;

        GRBVar full_depth = model.addVar(0, n - 1, 0.0, GRB_CONTINUOUS);
        for(int i = 0; i < n; i++)
            model.addConstr(full_depth >= depth_vars[i]);

        cout << "Setting objective." << endl;
        GRBLinExpr obj = full_depth;
        model.setObjective(obj, GRB_MINIMIZE);

        cout << "Optimizing." << endl;
        model.optimize();

        cout << "'Parent of' table" << endl;
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                cout << (parent_of_vars[i][j].get(GRB_DoubleAttr_X) > 0.5) << " ";
            }
            cout << endl;
        }

        cout << "'Ascendant of' table" << endl;
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                cout << (asc_of_vars[i][j].get(GRB_DoubleAttr_X) > 0.5) << " ";
            }
            cout << endl;
        }

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

}

