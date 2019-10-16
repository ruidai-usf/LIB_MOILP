/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: david
 *
 * Created on January 30, 2017, 1:19 PM
 */

#include <cstdlib>
#include <iostream> /*This is a package for reading/writing from/in screen in C++ */
#include <fstream> /*This is a package for reading/writing files C++ */
#include <ilcplex/ilocplex.h>
#include <vector> /*This is a package for creating dynamic lists */
#include <ctime>

using namespace std;

#define epsilon 0.001
#define Infinity 1000000.0

/*We Now Declare Our Global C++ Variables*/

double start_time(0);
double finish_time(0);

int Constraints_num;
int Is_one_ND_point(0);
int Var_num;
int N_Objectives(2);

double* Obj1_coef;
double* Obj2_coef;

double** Var_coef;

double* RHS_value;
int* Type_value;

/*We Now Declare Our CPLEX Variables*/

ILOSTLBEGIN
IloEnv env;
IloModel WS_IP_model(env); /* This Is Your Model Name */
IloCplex WS_IP_cplex(WS_IP_model); /* This is Your Solver Name */
IloNumVarArray var_x;
IloExpr Math_Expression(env);
IloObjective Objective_Func(env);
IloRangeArray Constraints(env);
IloRangeArray Extra_Constraints(env);

void Reading_The_Input_File(char* Input_file_name) {
    /* For reading from a file we should use the ifstream function */

    ifstream InputFile;
    InputFile.open(Input_file_name);

    InputFile >> Constraints_num;
    InputFile >> Var_num;

    /* You can now create the coefficients of objective function */

    Obj1_coef = new double [Var_num];

    for (int i = 0; i < Var_num; i++) {
        InputFile >> Obj1_coef[i];
    }

    Obj2_coef = new double [Var_num];

    for (int i = 0; i < Var_num; i++) {
        InputFile >> Obj2_coef[i];
    }

    /* You can now create the coefficient of A matrix and b */

    Var_coef = new double* [Constraints_num];

    for (int i = 0; i < Constraints_num; i++) {

        Var_coef[i] = new double [Var_num];
        for (int j = 0; j < Var_num; j++) {
            InputFile >> Var_coef[i][j];
        }
    }

    RHS_value = new double[Constraints_num];
    for (int i = 0; i < Constraints_num; i++) {
        InputFile >> RHS_value[i];
    }

    Type_value = new int[Constraints_num];
    for (int i = 0; i < Constraints_num; i++) {
        InputFile >> Type_value[i];
    }

    /* We are done with Input File, so we close it.  */

    InputFile.close();
}

void Create_The_CPLEX_Model() {

    /* The following will be used to make the name of variables compatible with our model */

    char VarName[25];

    /*Creating continuous var_x */

    var_x = IloNumVarArray(env, Var_num);

    for (int i = 0; i < Var_num; i++) {
        sprintf(VarName, "x(%d)", i + 1);
        var_x[i] = IloNumVar(env, 0.0, Infinity, ILOINT, VarName);
    }

    /* Creating the Constraints */

    for (int i = 0; i < Constraints_num; i++) {

        Math_Expression.clear();

        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Var_coef[i][j] * var_x[j];
        }

        Math_Expression += (-1) * RHS_value[i];

        //        cout << Type_value[i] << endl;

        if (Type_value[i] == 1) {
            Constraints.add(Math_Expression <= 0);
        } else {
            Constraints.add(Math_Expression == 0);
        }

    }

    /* Creating the Objective Function */

    Math_Expression.clear();

    /* Adding the Objective Function and Constraints to the Model */

    WS_IP_model.add(Constraints);
}

class Box {
public:

    /*The top and bottom vertexes of the box*/

    double* vertex_T;
    double* vertex_B;

    /*The nondominated points returned*/

    double* point_found;
    int* point_solution;
    int is_ND_point;

    /*The following is a constructor of a node that initialize a value for variables of the node */
    Box() {
        vertex_T = new double [N_Objectives];
        vertex_B = new double [N_Objectives];
        point_found = new double [N_Objectives];
        point_solution = new int [Var_num];

        is_ND_point = 0;

        for (int i = 0; i < N_Objectives; i++) {
            vertex_T[i] = 0;
            vertex_B[i] = 0;
            point_found[i] = 0;
        }

        for (int i = 0; i < Var_num; i++) {
            point_solution[i] = 0;
        }

    }

    /*We now update the information about this node based on its parent's information */

    void Generate_New_Nondominated_Points() {

        Math_Expression.clear();

        for (int j = 0; j < Var_num; j++) {
            Math_Expression += ((vertex_T[1] - vertex_B[1]) * Obj1_coef[j]+(vertex_B[0]
                    - vertex_T[0]) * Obj2_coef[j]) * var_x[j];
        }

        Objective_Func = IloMinimize(env, Math_Expression);

        WS_IP_model.add(Objective_Func);

        WS_IP_cplex.extract(WS_IP_model);

        WS_IP_cplex.setOut(env.getNullStream());
        if (WS_IP_cplex.solve() == 0) {
            cout << "Infeasible" << endl;
        }

        for (int i = 0; i < Var_num; i++) {
            point_solution[i] = WS_IP_cplex.getValue(var_x[i]);
        }

        for (int j = 0; j < Var_num; j++) {
            point_found[0] += Obj1_coef[j] * point_solution[j];
            point_found[1] += Obj2_coef[j] * point_solution[j];
        }


        double weighted_obj_val_vertex;
        double weighted_obj_val_point;

        weighted_obj_val_vertex = (vertex_T[1] - vertex_B[1]) * vertex_T[0]+(vertex_B[0]
                - vertex_T[0]) * vertex_T[1];

        for (int j = 0; j < Var_num; j++) {
            weighted_obj_val_point += ((vertex_T[1] - vertex_B[1]) * Obj1_coef[j]+(vertex_B[0]
                    - vertex_T[0]) * Obj2_coef[j]) * point_solution[j];
        }
        //        cout<<weighted_obj_val_vertex<<endl;
        //        cout<<weighted_obj_val_point<<endl;

        if (weighted_obj_val_point <= weighted_obj_val_vertex - epsilon) {
            is_ND_point = 1;
        }

        WS_IP_cplex.clear();
        WS_IP_model.remove(Objective_Func);

    }

    virtual ~Box() {
        delete [] vertex_T;
        delete [] vertex_B;
        delete [] point_found;
        delete [] point_solution;
    }
private:

};


vector <Box*>Tree_of_Boxes;

void Create_The_Box_ZERO() {

    /* Create a New Node and call it New_Box */
    Box* Initial_Box = new Box;

    /*get the value of Vertex_B[1]*/

    Math_Expression.clear();

    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj2_coef[j] * var_x[j];
    }

    Objective_Func = IloMinimize(env, Math_Expression);

    /*add objective function into the model*/
    WS_IP_model.add(Objective_Func);

    /*generate present cplex model*/
    WS_IP_cplex.extract(WS_IP_model);

//    WS_IP_cplex.exportModel("Model1.lp");

    WS_IP_cplex.setOut(env.getNullStream());
    if (WS_IP_cplex.solve() == 0) {
        cout << "Infeasible" << endl;
    }

    Initial_Box->vertex_B[1] = WS_IP_cplex.getObjValue();

    //    cout << Initial_Box->vertex_B[1] << endl;

    WS_IP_cplex.clear();
    WS_IP_model.remove(Objective_Func);

    /*get the value of vertex_B[0]*/

    Math_Expression.clear();

    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj1_coef[j] * var_x[j];
    }

    Objective_Func = IloMinimize(env, Math_Expression);

    WS_IP_model.add(Objective_Func);

    Math_Expression.clear();

    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj2_coef[j] * var_x[j];
    }

    Extra_Constraints.add(Math_Expression - Initial_Box->vertex_B[1] <= epsilon);

    WS_IP_model.add(Extra_Constraints);

    WS_IP_cplex.extract(WS_IP_model);

//    WS_IP_cplex.exportModel("Model2.lp");

    WS_IP_cplex.setOut(env.getNullStream());
    if (WS_IP_cplex.solve() == 0) {
        cout << "Infeasible" << endl;
    }

    Initial_Box->vertex_B[0] = WS_IP_cplex.getObjValue();

    //    cout << Initial_Box->vertex_B[0] << endl;

    WS_IP_cplex.clear();
    WS_IP_model.remove(Extra_Constraints);
    Extra_Constraints.clear();

    /*get the value of vertex_T[0]*/

    WS_IP_cplex.extract(WS_IP_model);

//    WS_IP_cplex.exportModel("Model3.lp");

    WS_IP_cplex.setOut(env.getNullStream());
    if (WS_IP_cplex.solve() == 0) {
        cout << "Infeasible" << endl;
    }

    Initial_Box->vertex_T[0] = WS_IP_cplex.getObjValue();

    //    cout << Initial_Box->vertex_T[0] << endl;

    WS_IP_cplex.clear();
    WS_IP_model.remove(Objective_Func);

    if (Initial_Box->vertex_B[0] == Initial_Box->vertex_T[0]) {
        Is_one_ND_point = 1;
    }

    /*get the value of vertex_T[1]*/
    if (Is_one_ND_point == 0) {

        Math_Expression.clear();

        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj2_coef[j] * var_x[j];
        }

        Objective_Func = IloMinimize(env, Math_Expression);

        WS_IP_model.add(Objective_Func);

        Math_Expression.clear();

        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj1_coef[j] * var_x[j];
        }

        Extra_Constraints.add(Math_Expression - Initial_Box->vertex_T[0] <= epsilon);

        WS_IP_model.add(Extra_Constraints);

        WS_IP_cplex.extract(WS_IP_model);

//        WS_IP_cplex.exportModel("Model4.lp");

        WS_IP_cplex.setOut(env.getNullStream());
        if (WS_IP_cplex.solve() == 0) {
            cout << "Infeasible" << endl;
        }

        Initial_Box->vertex_T[1] = WS_IP_cplex.getObjValue();

        //        cout << Initial_Box->vertex_T[1] << endl;

        WS_IP_cplex.clear();

        WS_IP_model.remove(Objective_Func);
        WS_IP_model.remove(Extra_Constraints);
        Extra_Constraints.clear();
    }

    Tree_of_Boxes.push_back(Initial_Box);

}

vector <double*> ND_points_set;

void Add_the_new_nondominated_point_to_the_set(double* New_nondominated_point) {
    bool It_is_Added(0);

    for (int i = 0; i < ND_points_set.size(); i++) {
        if (New_nondominated_point[0] <= ND_points_set.at(i)[0] - epsilon) {
            ND_points_set.insert(ND_points_set.begin() + i, New_nondominated_point);
            It_is_Added = 1;
            break;
        } else {
            if (New_nondominated_point[0] == ND_points_set.at(i)[0]) {
                cout << "Warning: a nondominated point is generated which we did not add to the set" << endl;
                It_is_Added = 1;
                break;
            }
        }
    }

    if (It_is_Added == 0) {
        ND_points_set.push_back(New_nondominated_point);
    }

}

void Writing_The_Output_File(char* Report_file, char* ND_set_file, char* Case_name) {
    finish_time = clock();

    ofstream OutputFile;
    OutputFile.open(ND_set_file);
//    cout << ND_points_set.size() << endl;
    for (int i = 0; i < ND_points_set.size(); i++) {
        OutputFile << ND_points_set.at(i)[0] << "    " << ND_points_set.at(i)[1] << endl;
    }
    OutputFile << endl;
    OutputFile.close();

    OutputFile.open(Report_file, ios::app);
    OutputFile << Case_name << " ";
    OutputFile << (finish_time - start_time) / CLOCKS_PER_SEC << " ";
    OutputFile << ND_points_set.size() << " ";

    OutputFile << endl;

    OutputFile.close();
}

int main(int argc, char** argv) {

    start_time = clock();
    Reading_The_Input_File(argv[1]);
    Create_The_CPLEX_Model();
    Create_The_Box_ZERO();

    double* Initial1_point = new double [N_Objectives];

    for (int i = 0; i < N_Objectives; i++) {
        Initial1_point[i] = Tree_of_Boxes.front()->vertex_B[i];
    }

    ND_points_set.push_back(Initial1_point);

    if (Is_one_ND_point == 0) {

        double* Initial2_point = new double [N_Objectives];

        for (int i = 0; i < N_Objectives; i++) {
            Initial2_point[i] = Tree_of_Boxes.front()->vertex_T[i];
        }

        ND_points_set.insert(ND_points_set.begin(), Initial2_point);

        while (Tree_of_Boxes.size() > 0) {
            Box* temp1_box = new Box;
            Box* temp2_box = new Box;

            Tree_of_Boxes.front()->Generate_New_Nondominated_Points();
            //            cout<<Tree_of_Boxes.front()->point_found[0]<<"; "<<Tree_of_Boxes.front()->point_found[1]<<endl;

            if (Tree_of_Boxes.front()->is_ND_point == 1) {
                
                double* New_point = new double [N_Objectives];
                for (int i = 0; i < N_Objectives; i++) {
                    New_point[i] = Tree_of_Boxes.front()->point_found[i];
                    //                cout << New1_point[i] << "   ";
                }
                Add_the_new_nondominated_point_to_the_set(New_point);
                //                cout<<Tree_of_Boxes.front()->point_found[0]<<"; "<<Tree_of_Boxes.front()->point_found[1]<<endl;

                temp1_box->vertex_B[0] = Tree_of_Boxes.front()->point_found[0];
                temp1_box->vertex_B[1] = Tree_of_Boxes.front()->point_found[1];
                temp1_box->vertex_T[0] = Tree_of_Boxes.front()->vertex_T[0];
                temp1_box->vertex_T[1] = Tree_of_Boxes.front()->vertex_T[1];

                temp2_box->vertex_B[0] = Tree_of_Boxes.front()->vertex_B[0];
                temp2_box->vertex_B[1] = Tree_of_Boxes.front()->vertex_B[1];
                temp2_box->vertex_T[0] = Tree_of_Boxes.front()->point_found[0];
                temp2_box->vertex_T[1] = Tree_of_Boxes.front()->point_found[1];

                Tree_of_Boxes.push_back(temp1_box);
                Tree_of_Boxes.push_back(temp2_box);
            } else {
                temp1_box->~Box();
                temp2_box->~Box();
            }

            Tree_of_Boxes.front()->~Box();
            Tree_of_Boxes.erase(Tree_of_Boxes.begin());
        }
    }

    Writing_The_Output_File(argv[2], argv[3], argv[4]);

    return 0;
}

