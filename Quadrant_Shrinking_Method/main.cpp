/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: david
 *
 * Created on March 23, 2017, 4:23 PM
 */

#include <cstdlib>
#include <iostream> /*This is a package for reading/writing from/in screen in C++ */
#include <fstream> /*This is a package for reading/writing files C++ */
#include <ilcplex/ilocplex.h>
#include <vector> /*This is a package for creating dynamic lists */
#include <ctime>

using namespace std;

#define epsilon 0.9
#define small_value 0.01
#define Infinity 1000000.0

/*We Now Declare Our Global C++ Variables*/

double start_time(0);
double finish_time(0);

int Constraints_num;
int Var_num;
int N_Objectives(3);

double* Obj1_coef;
double* Obj2_coef;
double* Obj3_coef;

double** Var_coef;

double* RHS_value;
int* Type_value;

/*We Now Declare Our CPLEX Variables*/

ILOSTLBEGIN
IloEnv env;
IloModel QSM_model(env); /* This Is Your Model Name */
IloCplex QSM_cplex(QSM_model); /* This is Your Solver Name */
IloNumVarArray var_x;

IloNumVarArray startVar(env);
IloNumArray startVal(env);

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

    Obj3_coef = new double [Var_num];

    for (int i = 0; i < Var_num; i++) {
        InputFile >> Obj3_coef[i];
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

    QSM_model.add(Constraints);

    /*import the decision variables into the start variable*/

    for (int j = 0; j < Var_num; j++) {

        startVar.add(var_x[j]);
    }
}

/*define class to save information of every Quadrant*/
class Quadrant {
public:

    /*define Quadrant(U)*/
    double* quadrant_U;
    double* quadrant_U_candidate1;
    double* quadrant_U_candidate2;

    /*indicate whether a unknown ND point is found*/
    int is_empty;

    /*ND_point information*/
    double* point_found;
    int* point_solution;

    /*The following is a constructor of a node that initialize a value for variables of the node */
    Quadrant() {
        quadrant_U = new double [N_Objectives - 1];

        quadrant_U_candidate1 = new double [N_Objectives - 1]; //always add the candidate1 first then add the candidate2
        quadrant_U_candidate2 = new double [N_Objectives - 1];

        point_found = new double [N_Objectives];
        point_solution = new int [Var_num];

        is_empty = 0;


        for (int i = 0; i < (N_Objectives - 1); i++) {
            quadrant_U[i] = 0;
            quadrant_U_candidate1[i] = 0;
            quadrant_U_candidate2[i] = 0;
        }

        for (int i = 0; i < N_Objectives; i++) {
            point_found[i] = 0;
        }

        for (int i = 0; i < Var_num; i++) {
            point_solution[i] = 0;
        }

    }

    virtual ~Quadrant() {
        delete [] quadrant_U;
        delete [] point_found;
        delete [] point_solution;
        delete [] quadrant_U_candidate1;
        delete [] quadrant_U_candidate2;
    }
private:

};

void Searching_Quadrant_For_Right_Boundary(Quadrant* Target_quadrant) {

    /*add objective function for first stage*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj3_coef[j] * var_x[j];
    }
    Objective_Func = IloMinimize(env, Math_Expression);

    QSM_model.add(Objective_Func);

    /*add constraint for Quadrant(U_1)*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj1_coef[j] * var_x[j];
    }
    Extra_Constraints.add(Math_Expression - Target_quadrant->quadrant_U[0] <= 0);

    /*add constraint for Quadrant(U_2)*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj2_coef[j] * var_x[j];
    }
    Extra_Constraints.add(Math_Expression - Target_quadrant->quadrant_U[1] <= 0);

    QSM_model.add(Extra_Constraints);

    /*extract the model*/
    QSM_cplex.extract(QSM_model);

    //    QSM_cplex.exportModel("1.lp");

    QSM_cplex.setOut(env.getNullStream());
    if (QSM_cplex.solve() == 0) {
        Target_quadrant->is_empty = 1;
        cout << "Infeasible" << endl;
    }

    if (Target_quadrant->is_empty == 0) {
        int* temp_solution;
        double* temp_ND_point;

        temp_solution = new int [Var_num];
        for (int i = 0; i < Var_num; i++) {
            startVal.add(QSM_cplex.getValue(var_x[i]));
            //            cout << QSM_cplex.getValue(var_x[i]) << "; ";
            temp_solution[i] = QSM_cplex.getValue(var_x[i]);
        }
        //        cout << endl;

        temp_ND_point = new double [N_Objectives];
        for (int j = 0; j < Var_num; j++) {
            temp_ND_point[0] += Obj1_coef[j] * temp_solution[j];
            temp_ND_point[1] += Obj2_coef[j] * temp_solution[j];
            temp_ND_point[2] += Obj3_coef[j] * temp_solution[j];
        }

        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();

        /*add objective function for second stage*/
        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += (Obj1_coef[j] + Obj2_coef[j] + Obj3_coef[j]) * var_x[j];
        }
        Objective_Func = IloMinimize(env, Math_Expression);

        QSM_model.add(Objective_Func);

        /*add extra constraints for second stage*/
        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj1_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[0] <= 0);

        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj2_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[1] <= 0);

        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj3_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[2] <= 0);

        QSM_model.add(Extra_Constraints);

        /*extract the model*/
        QSM_cplex.extract(QSM_model);

        //        QSM_cplex.exportModel("2.lp");

        /*feed feasible solution to cplex solver*/
        QSM_cplex.addMIPStart(startVar, startVal);

        QSM_cplex.setOut(env.getNullStream());

        if (QSM_cplex.solve() == 0) {
            cout << "The second stage: Infeasible" << endl;
        }

        for (int i = 0; i < Var_num; i++) {
            Target_quadrant->point_solution[i] = QSM_cplex.getValue(var_x[i]);
        }

        for (int j = 0; j < Var_num; j++) {
            Target_quadrant->point_found[0] += Obj1_coef[j] * Target_quadrant->point_solution[j];
            Target_quadrant->point_found[1] += Obj2_coef[j] * Target_quadrant->point_solution[j];
            Target_quadrant->point_found[2] += Obj3_coef[j] * Target_quadrant->point_solution[j];
        }

        Target_quadrant->quadrant_U_candidate1[0] = Target_quadrant->point_found[0] - epsilon;
        Target_quadrant->quadrant_U_candidate1[1] = Target_quadrant->quadrant_U[1];

        Target_quadrant->quadrant_U_candidate2[0] = Target_quadrant->quadrant_U[0];
        Target_quadrant->quadrant_U_candidate2[1] = Target_quadrant->point_found[1] - epsilon;

        startVal.clear();

        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();

        delete [] temp_solution;
        delete [] temp_ND_point;
    } else {
        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();
    }
}

void Searching_Quadrant_For_Top_Boundary(Quadrant* Target_quadrant) {

    /*add objective function for first stage*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj3_coef[j] * var_x[j];
    }
    Objective_Func = IloMinimize(env, Math_Expression);

    QSM_model.add(Objective_Func);

    /*add constraint for Quadrant(U_1)*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj1_coef[j] * var_x[j];
    }
    Extra_Constraints.add(Math_Expression - Target_quadrant->quadrant_U[0] <= 0);

    /*add constraint for Quadrant(U_2)*/
    Math_Expression.clear();
    for (int j = 0; j < Var_num; j++) {
        Math_Expression += Obj2_coef[j] * var_x[j];
    }
    Extra_Constraints.add(Math_Expression - Target_quadrant->quadrant_U[1] <= 0);

    QSM_model.add(Extra_Constraints);

    /*extract the model*/
    QSM_cplex.extract(QSM_model);

    //    QSM_cplex.exportModel("T1.lp");

    QSM_cplex.setOut(env.getNullStream());
    if (QSM_cplex.solve() == 0) {
        Target_quadrant->is_empty = 1;
        cout << "Infeasible" << endl;
    }

    if (Target_quadrant->is_empty == 0) {
        int* temp_solution;
        double* temp_ND_point;

        temp_solution = new int [Var_num];
        for (int i = 0; i < Var_num; i++) {
            startVal.add(QSM_cplex.getValue(var_x[i]));
            temp_solution[i] = QSM_cplex.getValue(var_x[i]);
        }

        temp_ND_point = new double [N_Objectives];
        for (int j = 0; j < Var_num; j++) {
            temp_ND_point[0] += Obj1_coef[j] * temp_solution[j];
            temp_ND_point[1] += Obj2_coef[j] * temp_solution[j];
            temp_ND_point[2] += Obj3_coef[j] * temp_solution[j];
        }

        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();

        /*add objective function for second stage*/
        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += (Obj1_coef[j] + Obj2_coef[j] + Obj3_coef[j]) * var_x[j];
        }
        Objective_Func = IloMinimize(env, Math_Expression);

        QSM_model.add(Objective_Func);

        /*add extra constraints for second stage*/
        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj1_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[0] <= 0);

        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj2_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[1] <= 0);

        Math_Expression.clear();
        for (int j = 0; j < Var_num; j++) {
            Math_Expression += Obj3_coef[j] * var_x[j];
        }
        Extra_Constraints.add(Math_Expression - temp_ND_point[2] <= 0);

        QSM_model.add(Extra_Constraints);

        /*extract the model*/
        QSM_cplex.extract(QSM_model);

        //        QSM_cplex.exportModel("T2.lp");

        /*feed feasible solution to cplex solver*/
        QSM_cplex.addMIPStart(startVar, startVal);

        QSM_cplex.setOut(env.getNullStream());

        if (QSM_cplex.solve() == 0) {
            cout << "The second stage: Infeasible" << endl;
        }

        for (int i = 0; i < Var_num; i++) {
            Target_quadrant->point_solution[i] = QSM_cplex.getValue(var_x[i]);
        }

        for (int j = 0; j < Var_num; j++) {
            Target_quadrant->point_found[0] += Obj1_coef[j] * Target_quadrant->point_solution[j];
            Target_quadrant->point_found[1] += Obj2_coef[j] * Target_quadrant->point_solution[j];
            Target_quadrant->point_found[2] += Obj3_coef[j] * Target_quadrant->point_solution[j];
        }

        Target_quadrant->quadrant_U_candidate1[0] = Target_quadrant->quadrant_U[0];
        Target_quadrant->quadrant_U_candidate1[1] = Target_quadrant->point_found[1] - epsilon;

        Target_quadrant->quadrant_U_candidate2[0] = Target_quadrant->point_found[0] - epsilon;
        Target_quadrant->quadrant_U_candidate2[1] = Target_quadrant->quadrant_U[1];

        startVal.clear();

        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();

        delete [] temp_solution;
        delete [] temp_ND_point;
    } else {
        QSM_cplex.clear();
        QSM_model.remove(Objective_Func);
        QSM_model.remove(Extra_Constraints);
        Extra_Constraints.clear();
    }
}


vector <Quadrant*> Quadrant_List;

vector <double*> List_D;
vector <double*> ND_points_set;

void Add_the_new_nondominated_point_to_the_set(double* New_nondominated_point) {
    bool It_is_Repeated(0);
    bool It_is_Added(0);

    for (int i = 0; i < ND_points_set.size(); i++) {
        if (New_nondominated_point[0] >= ND_points_set.at(i)[0] - small_value &&
                New_nondominated_point[0] <= ND_points_set.at(i)[0] + small_value &&
                New_nondominated_point[1] >= ND_points_set.at(i)[1] - small_value &&
                New_nondominated_point[1] <= ND_points_set.at(i)[1] + small_value &&
                New_nondominated_point[2] >= ND_points_set.at(i)[2] - small_value &&
                New_nondominated_point[2] <= ND_points_set.at(i)[2] + small_value) {
            It_is_Repeated = 1;
            break;
        }
    }
    //    cout<<It_is_Repeated<<endl;

    if (It_is_Repeated == 0) {
        for (int i = 0; i < ND_points_set.size(); i++) {
            if (New_nondominated_point[2] <= ND_points_set.at(i)[2] - small_value) {
                ND_points_set.insert(ND_points_set.begin() + i, New_nondominated_point);
                It_is_Added = 1;
                break;
            } else {
                if (New_nondominated_point[2] == ND_points_set.at(i)[2]) {
                    if (New_nondominated_point[1] <= ND_points_set.at(i)[1] - small_value) {
                        ND_points_set.insert(ND_points_set.begin() + i, New_nondominated_point);
                        It_is_Added = 1;
                        break;
                    }
                }
            }
        }

        if (It_is_Added == 0) {
            ND_points_set.push_back(New_nondominated_point);
        }
    }

}

void Writing_The_Output_File(char* Report_file, char* ND_set_file, char* Case_name) {
    finish_time = clock();

    ofstream OutputFile;
    OutputFile.open(ND_set_file);
    //    cout << ND_points_set.size() << endl;

    //    for (int i = 0; i < Quadrant_List.size(); i++) {
    //        OutputFile << Quadrant_List.at(i)->quadrant_U[0] << "; " << Quadrant_List.at(i)->quadrant_U[1] << "       "
    //                << Quadrant_List.at(i)->point_found[0] << "  " << Quadrant_List.at(i)->point_found[1]
    //                << "  " << Quadrant_List.at(i)->point_found[2] << "   |||   " << Quadrant_List.at(i)->quadrant_U_candidate1[0]
    //                << "  " << Quadrant_List.at(i)->quadrant_U_candidate1[1] << " ; " << Quadrant_List.at(i)->quadrant_U_candidate2[0]
    //                << "  " << Quadrant_List.at(i)->quadrant_U_candidate2[1] << endl;
    //        for (int j = 0; j < Var_num; j++) {
    //            OutputFile << Quadrant_List.at(i)->point_solution[j] << "; ";
    //        }
    //        OutputFile << endl;
    //    }

    for (int i = 0; i < ND_points_set.size(); i++) {
        OutputFile << ND_points_set.at(i)[0] << "  " << ND_points_set.at(i)[1] << "  " << ND_points_set.at(i)[2] << endl;
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

    double* Initial_quadrant_U = new double [N_Objectives - 1];

    for (int i = 0; i < (N_Objectives - 1); i++) {
        Initial_quadrant_U[i] = Infinity;
    }

    List_D.push_back(Initial_quadrant_U);

    while (List_D.size() > 0) {

        /*right boundary searching*/
        int right_boundary_not_done = 1;
        while (right_boundary_not_done == 1) {
            //            cout << "Right: " << List_D.front()[0] << "; " << List_D.front()[1] << endl;

            Quadrant* temp_right_quadrant = new Quadrant;

            for (int i = 0; i < (N_Objectives - 1); i++) {
                temp_right_quadrant->quadrant_U[i] = List_D.front()[i];
            }

            delete [] List_D.front();
            List_D.erase(List_D.begin());

            Searching_Quadrant_For_Right_Boundary(temp_right_quadrant);
            //            cout<<"SQRight!!!"<<endl;

            Quadrant_List.push_back(temp_right_quadrant);

            if (temp_right_quadrant->is_empty == 1) {
                right_boundary_not_done = 0;
            } else {
                double* new_found_ND_point = new double [N_Objectives];

                for (int i = 0; i < N_Objectives; i++) {
                    new_found_ND_point[i] = temp_right_quadrant->point_found[i];
                }
                Add_the_new_nondominated_point_to_the_set(new_found_ND_point);

                //                cout << List_D.size() << endl;
                if (temp_right_quadrant->quadrant_U_candidate1[0] >= List_D.front()[0]
                        + small_value || List_D.size() <= 0) {
                    double* quadrant_U_candidate_temp1 = new double [N_Objectives - 1];
                    for (int i = 0; i < (N_Objectives - 1); i++) {
                        quadrant_U_candidate_temp1[i] = temp_right_quadrant->quadrant_U_candidate1[i];
                    }
                    List_D.insert(List_D.begin(), quadrant_U_candidate_temp1);
                }

                double* quadrant_U_candidate_temp2 = new double [N_Objectives - 1];
                for (int i = 0; i < (N_Objectives - 1); i++) {
                    quadrant_U_candidate_temp2[i] = temp_right_quadrant->quadrant_U_candidate2[i];
                }
                List_D.insert(List_D.begin(), quadrant_U_candidate_temp2);
            }
        }

        /*top boundary searching*/
        int top_boundary_not_done = 1;
        if (List_D.size() <= 0) {
            break;
        } else {
            while (top_boundary_not_done == 1) {

                Quadrant* temp_top_quadrant = new Quadrant;

                for (int i = 0; i < (N_Objectives - 1); i++) {
                    temp_top_quadrant->quadrant_U[i] = List_D.back()[i];
                }

                delete [] List_D.back();
                List_D.pop_back();

                Searching_Quadrant_For_Top_Boundary(temp_top_quadrant);

                Quadrant_List.push_back(temp_top_quadrant);

                if (temp_top_quadrant->is_empty == 1) {
                    top_boundary_not_done = 0;
                } else {
                    double* new_found_ND_point = new double [N_Objectives];

                    for (int i = 0; i < N_Objectives; i++) {
                        new_found_ND_point[i] = temp_top_quadrant->point_found[i];
                    }
                    Add_the_new_nondominated_point_to_the_set(new_found_ND_point);

                    //                cout << List_D.size() << endl;
                    int first_push_back = 0;

                    if (List_D.size() <= 0) {
                        first_push_back = 1;
                    } else {
                        if (temp_top_quadrant->quadrant_U_candidate1[1] >= List_D.back()[1] + small_value) {
                            first_push_back = 1;
                        }
                    }

                    if (first_push_back == 1) {
                        double* quadrant_U_candidate_temp3 = new double [N_Objectives - 1];
                        for (int i = 0; i < (N_Objectives - 1); i++) {
                            quadrant_U_candidate_temp3[i] = temp_top_quadrant->quadrant_U_candidate1[i];
                        }
                        List_D.push_back(quadrant_U_candidate_temp3);
                    }

                    double* quadrant_U_candidate_temp4 = new double [N_Objectives - 1];
                    for (int i = 0; i < (N_Objectives - 1); i++) {
                        quadrant_U_candidate_temp4[i] = temp_top_quadrant->quadrant_U_candidate2[i];
                    }
                    List_D.push_back(quadrant_U_candidate_temp4);
                }
            }
        }
    }

    Writing_The_Output_File(argv[2], argv[3], argv[4]);

    return 0;
}

