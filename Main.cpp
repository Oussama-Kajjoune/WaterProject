
#include "libxl.h"
#include <sstream>
#include "ilcplex/ilocplex.h";
#include <iostream>
#include <windows.h>
#include <chrono>;
#include <random>;

using namespace libxl;

typedef IloArray<IloNumVarArray> NumVar2D;  // enables us to define 2-D decision variables
typedef IloArray<NumVar2D> NumVar3D;        // enables us to define 3-D decision variables
typedef IloArray<IloExprArray> Expr2D;      // enables us to define 2-D expression


Sheet* getSheetByName(Book* book, const wchar_t* name)
{
    for (int i = 0; i < book->sheetCount(); ++i)
    {
        if (wcscmp(book->getSheet(i)->name(), name) == 0)
        {
            return book->getSheet(i);
        }
    }
    return 0;
}

int main()
{
#pragma region Open_Excel_File
    Book* book = xlCreateXMLBookW();
    book->setKey(L"kajjo", L"windows-2a21240508c7ed0261ba6f61adj7n2e5");
    book->load(L"ProblemData.xlsx");
    Sheet* sheet_Input = getSheetByName(book, L"Parameters");

    auto start = std::chrono::high_resolution_clock::now();
#pragma endregion

#pragma region Read_Problem_Data_From_Excel
     
    int nT = sheet_Input->readNum(1, 2);
    int nR = sheet_Input->readNum(2, 2);
    int nS = sheet_Input->readNum(3, 2);
    int nN = sheet_Input->readNum(4, 2);
    int nL = sheet_Input->readNum(5, 2);
    int nK = sheet_Input->readNum(6, 2);
    int nP = sheet_Input->readNum(7, 2);
    int nO = sheet_Input->readNum(8, 2);
    int** F = new int* [nO];
    for (int n = 0; n < nO; n++)
    {
        F[n] = new int[nO];
        for (int m = 0; m < nO; m++)
        {
            F[n][m] = sheet_Input->readNum(10 + n, 3 + m);
        }
    }

    int t_s = sheet_Input->readNum(23, 2);
    int year = sheet_Input->readNum(24, 2);
    double rho = sheet_Input->readNum(25, 2);
    double beta = (pow(1 + rho, year) - 1) / (rho * pow(1 + rho, year));
    double g = sheet_Input->readNum(26, 2);
    double rho_w = sheet_Input->readNum(27, 2);
    double ec = sheet_Input->readNum(28, 2);
    double v_max = sheet_Input->readNum(29, 2);
    double v_min = sheet_Input->readNum(30, 2);
    double eta = sheet_Input->readNum(31, 2);
    double pi = sheet_Input->readNum(32, 2);

    double* tau = new double[nP];
    for (int n = 0; n < nP; n++)
    {
        tau[n] = sheet_Input->readNum(34, 2 + n);
    }
    double* sigma = new double[nP];
    for (int n = 0; n < nP; n++)
    {
        sigma[n] = sheet_Input->readNum(35, 2 + n);
    }
    double* epsilon = new double[nR];
    for (int n = 0; n < nR; n++)
    {
        epsilon[n] = sheet_Input->readNum(36, 2 + n);
    }
    double* gamma = new double[nR];
    for (int n = 0; n < nR; n++)
    {
        gamma[n] = sheet_Input->readNum(37, 2 + n);
    }
    double* theta = new double[nR + nS];
    for (int n = nR; n < nR + nS; n++)
    {
        theta[n] = sheet_Input->readNum(38, 2 + n - nR);
    }
    double* fix = new double[nR];
    for (int n = 0; n < nR; n++)
    {
        fix[n] = sheet_Input->readNum(39, 2 + n);
    }
    double* k = new double[nR + nS + nN];
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        k[n] = sheet_Input->readNum(40, 2 + n - nR - nS);
    }
    double** delta_Z = new double* [nO];
    for (int n = 0; n < nO; n++)
    {
        delta_Z[n] = new double[nO];
        for (int m = 0; m < nO; m++)
        {
            delta_Z[n][m] = sheet_Input->readNum(43 + n, 3 + m);
        }
    }
    double** e = new double* [nO];
    for (int n = 0; n < nO; n++)
    {
        e[n] = new double[nO];
        for (int m = 0; m < nO; m++)
        {
            e[n][m] = sheet_Input->readNum(57 + n, 3 + m);
        }
    }
    double** d_flat = new double* [nL * nK];
    for (int n = 0; n < nL * nK; n++)
    {
        d_flat[n] = new double[nT];
        for (int t = 0; t < nT; t++)
        {
            d_flat[n][t] = sheet_Input->readNum(71 + n, 4 + t);
        }
    }
    double*** d = new double** [nR + nS + nN + nL];
    for (int l = nR + nS + nN; l < nR + nS + nN + nL; l++)
    {
        d[l] = new double* [nK];
        for (int k = 0; k < nK; k++)
        {
            d[l][k] = new double[nT];
            for (int t = 0; t < nT; t++)
            {
                d[l][k][t] = d_flat[k + nK * (l - (nR + nS + nN))][t];
            }
        }
    }

#pragma endregion

#pragma region Create_Env_Model_Cplex
    IloEnv env;
    IloModel Model(env);
#pragma endregion

#pragma region Decision_variables
    NumVar2D Q_r(env, nR);
    for (int r = 0; r < nR; r++)
    { 
        Q_r[r] = IloNumVarArray(env, nT, 0.0, IloInfinity, ILOFLOAT);
    }
    NumVar2D I_s(env, nS);
    for (int s = nR; s < nR + nS; s++)
    {
        I_s[s] = IloNumVarArray(env, nT + 1, 0.0, IloInfinity, ILOFLOAT);
    }
    IloNumVarArray C_r_max(env, nR, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray I_s_max(env, nR + nS, 0.0, IloInfinity, ILOFLOAT);
    NumVar3D Q(env, nO);                                
    for (int i = 0; i < nO; i++)
    {
        Q[i] = NumVar2D(env, nO);
        for (int j = 0; j < nO; j++)
        {
            Q[i][j] = IloNumVarArray(env, nT, 0.0, IloInfinity, ILOFLOAT);
        }
    }
    NumVar2D Q_l(env, nR + nS + nN + nL);
    for (int l = nR + nS + nN; l < nR + nS + nN + nL; l++)
    {
        Q_l[l] = IloNumVarArray(env, nT, 0.0, IloInfinity, ILOFLOAT);
    }
    NumVar2D D(env, nR + nS + nN);
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        D[n] = IloNumVarArray(env, nR + nS + nN + nL, 0.0, IloInfinity, ILOFLOAT);
    }
    NumVar3D Y_prime(env, nR + nS + nN);
    for (int i = nR + nS; i < nR + nS + nN; i++)
    {
        Y_prime[i] = NumVar2D(env, nR + nS + nN + nL);
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            Y_prime[i][j] = IloNumVarArray(env, nP, 0.0, 1.0, ILOBOOL);
        }
    }
    NumVar2D Y(env, nR + nS + nN);
    for (int i = nR + nS; i < nR + nS + nN; i++)
    {
        Y[i] = IloNumVarArray(env, nR + nS + nN + nL, 0.0, 1.0, ILOBOOL);
    }
    IloNumVarArray W(env, nR + nS + nN, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray CECRO(env, nR, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray CECN(env, nR + nS + nN, 0.0, IloInfinity, ILOFLOAT);
    NumVar2D CECP(env, nR + nS + nN);
    for (int j = nR + nS; j < nR + nS + nN; j++)
    {
        CECP[j] = IloNumVarArray(env, nR + nS + nN + nL, 0.0, IloInfinity, ILOFLOAT);
    }
    IloNumVarArray CECS(env, nR + nS, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray ELC_r(env, nR, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray ELC_n(env, nR + nS + nN, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray OECRO(env, nR, 0.0, IloInfinity, ILOFLOAT);
    IloNumVarArray OECN(env, nR + nS + nN, 0.0, IloInfinity, ILOFLOAT);
    IloNumVar TC(env, 0.0, IloInfinity, ILOFLOAT);
#pragma endregion

#pragma region Expressions
    // To Constraints -->
#pragma endregion

#pragma region Objective_Function
    Model.add(IloMinimize(env, TC));
#pragma endregion

#pragma region Constraints
    // Flow_conservation_at_the_desalination_plant:
    for (int r = 0; r < nR; r++)
    {
        for (int t = 0; t < nT; t++)
        {
            IloExpr exp1(env);
            for (int j = 0; j < nO; j++)
            {
                if (F[r][j] == 1)
                {
                    exp1 += Q[r][j][t];
                }
            }
            Model.add(Q_r[r][t] == exp1);
        }
    }
    // Flow_conservation_at_pumping_stations:
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int t = 0; t < nT; t++)
        {
            IloExpr exp2(env);
            IloExpr exp3(env);
            for (int i = 0; i < nO; i++)
            {
                if (F[i][n] == 1)
                {
                    exp2 += Q[i][n][t];
                }
            }
            for (int j = 0; j < nO; j++)
            {
                if (F[n][j] == 1)
                {
                    exp3 += Q[n][j][t];
                }
            }
            Model.add(exp2 == exp3);
        }
    }
    // Flow_conservation_at_demand_locations:
    for (int l = nR + nS + nN; l < nR + nS + nN + nL; l++)
    {
        for (int t = 0; t < nT; t++)
        {
            IloExpr exp4(env);
            IloExpr exp5(env);
            for (int i = 0; i < nO; i++)
            {
                if (F[i][l] == 1)
                {
                    exp4 += Q[i][l][t];
                }
            }
            for (int j = 0; j < nO; j++)
            {
                if (F[l][j] == 1)
                {
                    exp5 += Q[l][j][t];
                }
            }
            Model.add(Q_l[l][t] == exp4 - exp5);
        }
    }
    // Satisfy_water_requirements:
    for (int l = nR + nS + nN; l < nR + nS + nN + nL; l++)
    {
        for (int t = 0; t < nT; t++)
        {
            IloExpr exp6(env);
            for (int k = 0; k < nK; k++)
            {
                exp6 += d[l][k][t];
            }
            Model.add(Q_l[l][t] * t_s >= exp6);
        }
    }
    // Initial_storage:
    for (int s = nR; s < nR + nS; s++)
    {
        Model.add(I_s[s][0] == 0);
    }
    // Flow_conservation_at_storage:
    for (int s = nR; s < nR + nS; s++)
    {
        for (int t = 0; t < nT; t++)
        {
            IloExpr exp7(env);
            for (int j = 0; j < nO; j++)
            {
                if (F[s][j] == 1)
                {
                    exp7 += t_s * Q[s][j][t];
                }
            }
            Model.add(I_s[s][t + 1] == I_s[s][t] + t_s * Q_r[0][t] - exp7);
        }
    }
    // Maximum_production:
    for (int r = 0; r < nR; r++)
    {
        for (int t = 0; t < nT; t++)
        {
            Model.add(Q_r[r][t] <= C_r_max[r]);
        }
    }
    // Maximum_storage:
    for (int s = nR; s < nR + nS; s++)
    {
        for (int t = 0; t < nT; t++)
        {
            Model.add(I_s[s][t] <= I_s_max[s]);
        }
    }
    // Selection_of_pipeline_diameter:
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            if (F[n][j] == 1)
            {
                IloExpr exp8(env);
                IloExpr exp9(env);
                for (int p = 0; p < nP; p++)
                {
                    exp8 += Y_prime[n][j][p];
                    exp9 += tau[p] * Y_prime[n][j][p];
                }
                Model.add(exp8 == Y[n][j]);
                Model.add(D[n][j] == exp9);
            }
        }
    }
    // Maximum_velocity_in_pipelines:
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int t = 0; t < nT; t++)
        {
            for (int j = nR + nS; j < nR + nS + nN + nL; j++)
            {
                if (F[n][j] == 1)
                {
                    IloExpr exp10(env);
                    for (int p = 0; p < nP; p++)
                    {
                        exp10 += tau[p] * tau[p] * Y_prime[n][j][p];
                    }
                    Model.add(Q[n][j][t] <= (pi * exp10 * 2 * v_max) / 4);
                }
            }
        }
    }
    // W
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        IloExpr exp11(env);
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            for (int t = 0; t < nT; t++)
            {
                if (F[n][j] == 1)
                {
                    exp11 += Q[n][j][t] * delta_Z[n][j];
                }
            }
        }
        Model.add(W[n] == ((rho_w * g) / eta) * exp11);
    }
    // CECRO
    for (int r = 0; r < nR; r++)
    {
        IloExpr exp12(env);
        for (int t = 0; t < nT; t++)
        {
            exp12 += Q_r[r][t];
        }
        Model.add(CECRO[r] == (gamma[r] * exp12) / beta);
    }
    // CECN
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        Model.add(CECN[n] == (k[n] * W[n]) / beta);
    }
    // CECP
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            IloExpr exp13(env);
            for (int p = 0; p < nP; p++)
            {
                if (F[n][j] == 1)
                {
                    exp13 += Y_prime[n][j][p] * sigma[p];
                }
            }
            Model.add(CECP[n][j] == (e[n][j] * exp13) / beta);
        }
    }
    // CECS
    for (int s = nR; s < nR + nS; s++)
    {
        Model.add(CECS[s] == (I_s_max[s] * theta[s]) / beta);
    }
    // ELC_r
    for (int r = 0; r < nR; r++)
    {
        IloExpr exp14(env);
        for (int t = 0; t < nT; t++)
        {
            exp14 += Q_r[r][t];
        }
        Model.add(ELC_r[r] == epsilon[r] * exp14);
    }
    // ELC_n
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        Model.add(ELC_n[n] == (t_s * nT) / 3600 * W[n]);
    }
    // OECRO
    for (int r = 0; r < nR; r++)
    {
        IloExpr exp15(env);
        for (int t = 0; t < nT; t++)
        {
            exp15 += Q_r[r][t];
        }
        Model.add(OECRO[r] == ec * ELC_r[r] + fix[r] * exp15);
    }
    // OECN
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        Model.add(OECN[n] == ec * ELC_n[n]);
    }
    // TC
    IloExpr exp16(env);
    IloExpr exp17(env);
    IloExpr exp18(env);
    IloExpr exp19(env);
    for (int r = 0; r < nR; r++)
    {
        exp16 += CECRO[r] + OECRO[r];
    }
    for (int s = nR; s < nR + nS; s++)
    {
        exp17 += CECS[s];
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        exp18 += CECN[n] + OECN[n];
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            if (F[n][j] == 1)
            {
                exp19 += CECP[n][j];
            }
        }
    }
    Model.add(TC == exp16 + exp17 + exp18 + exp19);
            
#pragma endregion

#pragma region Solve_Problem
    IloCplex cplex(Model);
    cplex.setOut(env.getNullStream());
    if (!cplex.solve()) {
        env.error() << "Failed to optimize the Master Problem!!!" << std::endl;
        throw(-1);
    }

    double obj = cplex.getObjValue();

    auto end = std::chrono::high_resolution_clock::now();
    auto Elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "\t Elapsed Time(ms): " << Elapsed.count() << std::endl;

    std::cout << "\n\n\t The objective value is: " << obj << std::endl;
#pragma endregion

#pragma region Write_Results_In_Excel
    Sheet* sheet_Output = getSheetByName(book, L"Results");
    sheet_Output->writeNum(0, 3, Elapsed.count());
    sheet_Output->writeNum(1, 3, obj);

    for (int r = 0; r < nR; r++)
    {
        for (int t = 0; t < nT; t++)
            sheet_Output->writeNum(3, 3 + t, cplex.getValue(Q_r[r][t]));
    }
    for (int s = nR; s < nR + nS; s++)
    {
        for (int t = 0; t < nT + 1; t++)
            sheet_Output->writeNum(4, 3 + t, cplex.getValue(I_s[s][t]));
    }
    for (int r = 0; r < nR; r++)
    {
        sheet_Output->writeNum(5, 3 + r, cplex.getValue(C_r_max[r]));
    }
    for (int s = nR; s < nR + nS; s++)
    {
        sheet_Output->writeNum(6, 3 + s - nR, cplex.getValue(I_s_max[s]));
    }
    for (int l = nR + nS + nN; l < nR + nS + nN + nL; l++)
    {
        for (int t = 0; t < nT; t++)
            sheet_Output->writeNum(7 + l - (nR + nS + nN), 3 + t, cplex.getValue(Q_l[l][t]));
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            if (F[n][j] == 1)
            {
                sheet_Output->writeNum(12 + n - (nR + nS), 3 + j, cplex.getValue(D[n][j]));
                sheet_Output->writeNum(17 + n - (nR + nS), 3 + j, cplex.getValue(Y[n][j]));
            }
        }
    }
    sheet_Output->writeNum(22, 3, beta);
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        sheet_Output->writeNum(23, 3 + n - (nR + nS), cplex.getValue(W[n]));
    }
    for (int r = 0; r < nR; r++)
    {
        sheet_Output->writeNum(24, 3 + r, cplex.getValue(CECRO[r]));
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        sheet_Output->writeNum(25, 3 + n - (nR + nS), cplex.getValue(CECN[n]));
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        for (int j = nR + nS; j < nR + nS + nN + nL; j++)
        {
            sheet_Output->writeNum(26 + n - (nR + nS), 3 + j, cplex.getValue(CECP[n][j]));
        }
    }
    for (int s = nR; s < nR + nS; s++)
    {
        sheet_Output->writeNum(31, 3 + s - nR, cplex.getValue(CECS[s]));
    }
    for (int r = 0; r < nR; r++)
    {
        sheet_Output->writeNum(32, 3 + r, cplex.getValue(ELC_r[r]));
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        sheet_Output->writeNum(33, 3 + n - (nR + nS), cplex.getValue(ELC_n[n]));
    }
    for (int r = 0; r < nR; r++)
    {
        sheet_Output->writeNum(34, 3 + r, cplex.getValue(OECRO[r]));
    }
    for (int n = nR + nS; n < nR + nS + nN; n++)
    {
        sheet_Output->writeNum(35, 3 + n - (nR + nS), cplex.getValue(OECN[n]));
    }
    sheet_Output->writeNum(36, 3, obj);
    for (int i = 0; i < nO; i++)
    {
        for (int j = 0; j < nO; j++)
        {
            for (int t = 0; t < nT; t++)
            {
                if (F[i][j] == 1)
                {
                    sheet_Output->writeNum(39 + j + i * nO, 3 + t, cplex.getValue(Q[i][j][t]));
                }
            }
        }
    }
#pragma endregion

#pragma region Save_Excel_File_and_Open_It
    if (book->save(L"ProblemData.xlsx"))
    {
        ::ShellExecute(NULL, L"open", L"ProblemData.xlsx", NULL, NULL, SW_SHOW);
    }
    else
    {
        std::cout << book->errorMessage() << std::endl;
    }

    book->release();
#pragma endregion

    return 0;
}