#include <iostream>
#include <fstream>
#include <conio.h>
#include <mtl/matrix.h>
#include <itl/preconditioner/ilu.h>
#include <itl/interface/mtl.h>
#include <itl/krylov/bicgstab.h>
#include "input_data.h"

using namespace std;
using namespace mtl;
using namespace itl;

//исходные данные
double ll;
double tt;
double a;
double alpha;
double kappaa;
double b;
double beta;
double thetab;
double thetainit;
double umin;
double umax;
double thetad;

//функционал качества
cost_func_t cost_func;

//начальное приближение
init_guess_t init_guess_type;

//параметры сетки
int N;
int tnum;

int ns; //размер рабочей области

//параметры метода Ньютона
const int numiter_newton = 10;
const double eps_newton = 1e-8;

//параметры итерационного метода решения СЛАУ
const int max_iter = 1000;
const double iter_rtol = 1e-10;


int neq;
int neq_phi;

typedef matrix<double, rectangle<>, array<compressed<> >,
    row_major>::type Matrix;

Matrix *A_phi, *A;
dense1D<double> *rhs_phi, *rhs;

typedef enum { THETA, PHI } var_t;

double ***u;
double ***phiinit;

double ****theta, ****phi;
double ***thetaold, ***thetanew, ***phinew;

double ***p1, ***p2;
double ***p1old, ***p2old;
double ***gr;

bool **vis_u;


int get_index(var_t v, int i, int j, int k)
{
    if (i < 0 || i > N || j < 0 || j > N || k < 0 || k > N) {
        cout << "\nBAD INDICES!!!\n";
        exit(1);
    }
    int r = i*(N + 1)*(N + 1) + j*(N + 1) + k;
    if (v == PHI)
        r += (N + 1)*(N + 1)*(N + 1);
    return r;
}

void fill_mat_phi(int i1, int j1, int k1, int i2, int j2, int k2, double num)
{
    int eq = get_index(THETA, i1, j1, k1);
    int c = get_index(THETA, i2, j2, k2);
    (*A_phi)(eq, c) = num;
}

void fill_rhs_phi(int i, int j, int k, double num)
{
    int eq = get_index(THETA, i, j, k);
    (*rhs_phi)[eq] = num;
}

void fill_mat(var_t v1, int i1, int j1, int k1,
    var_t v2, int i2, int j2, int k2, double num)
{
    int eq = get_index(v1, i1, j1, k1);
    int c = get_index(v2, i2, j2, k2);
    (*A)(eq, c) = num;
}

void fill_rhs(var_t v, int i, int j, int k, double num)
{
    int eq = get_index(v, i, j, k);
    (*rhs)[eq] = num;
}

int count_borders(int i, int j, int k)
{
    bool bi = i == 0 || i == N;
    bool bj = j == 0 || j == N;
    bool bk = k == 0 || k == N;
    return (int)bi + (int)bj + (int)bk;
}

void get_ij_border(int i, int j, int k,
    int &i1, int &j1, int &k1, int &i2, int &j2, int &k2,
    int &i3, int &j3, int &k3, int &i4, int &j4, int &k4,
    int &i5, int &j5, int &k5)
    //(i1, j1, k1) -- внутренний узел, остальные -- на той же грани
{
    int ind[3] = {i, j, k};
    int res[5][3];
    int cnt = 0;
    for (int s = 0; s < 3; ++s) {
        if (ind[s] == 0 || ind[s] == N) {
            res[0][s] = (ind[s] == 0 ? 1 : N-1);
            for (int q = 0; q < 3; ++q) {
                if (q != s)
                    res[0][q] = ind[q];
            }
        }
        else {
            for (int dind = -1; dind <= 1; dind += 2) {
                ++cnt;
                res[cnt][s] = ind[s] + dind;
                for (int q = 0; q < 3; ++q) {
                    if (q != s)
                        res[cnt][q] = ind[q];
                }
            }
        }
    }

    int *p[5][3] = {{&i1, &j1, &k1}, {&i2, &j2, &k2},
        {&i3, &j3, &k3}, {&i4, &j4, &k4}, {&i5, &j5, &k5}};
    for (int t = 0; t < 5; ++t) {
        for (int s = 0; s < 3; ++s)
            *p[t][s] = res[t][s];
    }
}

void get_ij_edge(int i, int j, int k,
    int &i1, int &j1, int &k1, int &i2, int &j2, int &k2,
    int &i3, int &j3, int &k3, int &i4, int &j4, int &k4)
    //(i1, j1, k1), (i2, j2, k2) -- узлы на смежных гранях,
    //(i3, j3, k3), (i4, j4, j4) -- на том же ребре
{
    int ind[3] = {i, j, k};
    int res[4][3];
    for (int s = 0; s < 3; ++s) {
        if (ind[s] > 0 && ind[s] < N) {
            res[2][s] = ind[s] - 1;
            res[3][s] = ind[s] + 1;
            for (int q = 0; q < 3; ++q) {
                if (q != s) {
                    res[2][q] = ind[q];
                    res[3][q] = ind[q];
                }
            }
            int cnt = 0;
            for (int q = 0; q < 3; ++q) {
                if (q != s) {
                    res[cnt][q] = (ind[q] == 0 ? 1 : N-1);
                    for (int t = 0; t < 3; ++t) {
                        if (t != q)
                            res[cnt][t] = ind[t];
                    }
                    ++cnt;
                }
            }
        }
    }

    int *p[4][3] = {{&i1, &j1, &k1}, {&i2, &j2, &k2},
        {&i3, &j3, &k3}, {&i4, &j4, &k4}};
    for (int t = 0; t < 4; ++t) {
        for (int s = 0; s < 3; ++s)
            *p[t][s] = res[t][s];
    }
}

void get_ij_corner(int i, int j, int k,
    int &i1, int &j1, int &k1, int &i2, int &j2, int &k2,
    int &i3, int &j3, int &k3)
    //(i1, j1, k1), (i2, j2, k2), (i3, j3, k3) -- узлы на смежных рёбрах
{
    int ind[3] = {i, j, k};
    int res[3][3];
    for (int s = 0; s < 3; ++s) {
        res[s][s] = (ind[s] == 0 ? 1 : N-1);
        for (int q = 0; q < 3; ++q) {
            if (q != s)
                res[s][q] = ind[q];
        }
    }

    int *p[3][3] = {{&i1, &j1, &k1}, {&i2, &j2, &k2},
        {&i3, &j3, &k3}};
    for (int t = 0; t < 3; ++t) {
        for (int s = 0; s < 3; ++s)
            *p[t][s] = res[t][s];
    }
}

double adj_src(int m, int i, int j, int k)
{
    if (cost_func == COST_FUNC_J1)
        return theta[m][i][j][k] - thetad;
    else
        return 0;
}

bool file_exists(string fname)
{
    return ifstream(fname.c_str()) != NULL;
}

string file_number(int num)
{
    ostringstream out;
    out << num;
    string s = out.str();
    if (num < 10)
        s = "0" + s;
    return s;
}

string control_file_name(int num)
{
    return file_number(num) + "_control.txt";
}

string k_file_name(int num)
{
    return file_number(num) + "_k.txt";
}

string jj_file_name(int num)
{
    return file_number(num) + "_jj.txt";
}

string crit_file_name(int num)
{
    return file_number(num) + "_crit.txt";
}

string theta_phi_file_name(int num)
{
    return file_number(num) + "_theta_phi.txt";
}

string log_file_name(int num)
{
    return file_number(num) + "_log.txt";
}

void fill_control()
/* Заполнение управления с рабочей области на всю границу.
   Рабочая область: 0 <= i <= ns, 0 <= j <= i.
*/
{
    //заполняем левый верхний квадрат
    for (int i = 0; i <= ns; ++i) {
        for (int j = 0; j <= i; ++j)
            u[j][i][0] = u[i][j][0];
    }
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            int i1 = (i <= ns ? i : N - i);
            int j1 = (j <= ns ? j : N - j);
            u[i][j][0] = u[i1][j1][0];
        }
    }
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            u[i][j][N] = u[i][j][0];
            u[0][i][j] = u[N][i][j] = u[i][j][0];
            u[i][0][j] = u[i][N][j] = u[i][j][0];
        }
    }
}

void correct_u(double &uval)
//поправка на точность
{
    if (fabs(uval - umin) < fabs(uval - umax))
        uval = umin;
    else
        uval = umax;
}

void read_control(int iter_num)
//читаем управление из файла и распространяем на всю границу
{
    ifstream fin(control_file_name(iter_num).c_str());
    for (int i = 0; i <= ns; ++i) {
        for (int j = 0; j <= i; ++j) {
            double dummy;
            fin >> u[i][j][0] >> dummy;
            correct_u(u[i][j][0]);
        }
    }
    fill_control();
}

void read_grad(int iter_num)
//читаем градиент из файла (только в рабочей области)
{
    ifstream fin(control_file_name(iter_num).c_str());
    for (int i = 0; i <= ns; ++i) {
        for (int j = 0; j <= i; ++j) {
            double dummy;
            fin >> dummy >> gr[i][j][0];
        }
    }
}

int get_prev_iterations()
{
    int prev_iterations = 0;
    while (file_exists(control_file_name(prev_iterations + 1)))
        ++prev_iterations;
    return prev_iterations;
}

bool the_same_control(int iter_num)
{
    //если на итерации iter_num было k != 1, то игнорируем
    int k;
    ifstream fin_k(k_file_name(iter_num).c_str());
    fin_k >> k;
    if (k != 1)
        return false;

    ifstream fin(control_file_name(iter_num).c_str());
    for (int i = 0; i <= ns; ++i) {
        for (int j = 0; j <= i; ++j) {
            double uval, dummy;
            fin >> uval >> dummy;
            correct_u(uval);
            if (uval != u[i][j][0])
                return false;
        }
    }
    return true;
}

bool **new_bool_array2d(int dim)
//выделяет память под двумерный массив dim x dim
{
    bool **arr = new bool*[dim];
    arr[0] = new bool[dim * dim];
    for (int i = 1; i < dim; ++i)
        arr[i] = arr[i-1] + dim;
    return arr;
}

void delete_bool_array2d(bool **arr, int dim)
//освобождает память под двумерный массив dim x dim
{
    delete[] arr[0];
    delete[] arr;
}

double ***new_double_array3d(int dim)
//выделяет память под трехмерный массив dim x dim x dim
{
    double ***arr = new double**[dim];
    arr[0] = new double*[dim * dim];
    for (int i = 1; i < dim; ++i)
        arr[i] = arr[i-1] + dim;

    arr[0][0] = new double[dim * dim * dim];
    for (int i = 0; i < dim; ++i) {
        if (i > 0)
            arr[i][0] = arr[i-1][0] + dim * dim;
        for (int j = 1; j < dim; ++j) {
            arr[i][j] = arr[i][j-1] + dim;
        }
    }

    return arr;
}

void delete_double_array3d(double ***arr, int dim)
//освобождает память под трехмерный массив dim x dim x dim
{
    delete[] arr[0][0];
    delete[] arr[0];
    delete[] arr;
}

double ****new_double_array4d(int dim0, int dim)
//выделяет память под четырехмерный массив dim0 x dim x dim x dim
{
    double ****arr = new double***[dim0];
    for (int i = 0; i < dim0; ++i)
        arr[i] = new_double_array3d(dim);
    return arr;
}

void delete_double_array4d(double ****arr, int dim0, int dim)
//освобождает память под четырехмерный массив dim0 x dim x dim x dim
{
    for (int i = 0; i < dim0; ++i)
        delete_double_array3d(arr[i], dim);
    delete[] arr;
}

int main(int argc, char *argv[])
{
    int iterations_number;

    if (argc != 2) {
        cout << "Usage: opt_control3d [num_iterations | -all]\n";
        exit(1);
    }
    else {
        if (string(argv[1]) == "-all") {
            //создание файлов all
            cout << "Print all\n";
            int prev_iterations = get_prev_iterations();
            {
                ofstream fout_jj("all_jj.txt");
                fout_jj.precision(17);
                for (int iter = 1; iter <= prev_iterations; ++iter) {
                    double jj;
                    ifstream fin_jj(jj_file_name(iter).c_str());
                    fin_jj >> jj;
                    fout_jj << iter << "   " << jj << "\n";
                }

                ofstream fout_crit("all_crit.txt");
                for (int iter = 1; iter <= prev_iterations; ++iter) {
                    int crit;
                    ifstream fin_crit(crit_file_name(iter).c_str());
                    fin_crit >> crit;
                    fout_crit << iter << "   " << crit << "\n";
                }

                ofstream fout_k("all_k.txt");
                for (int iter = 1; iter <= prev_iterations; ++iter) {
                    int k;
                    ifstream fin_k(k_file_name(iter).c_str());
                    fin_k >> k;
                    fout_k << iter << "   " << k << "\n";
                }

                ofstream fout_log("all_log.txt");
                for (int iter = 1; iter <= prev_iterations; ++iter) {
                    ifstream fin_log(log_file_name(iter).c_str());
                    while (!fin_log.eof()) {
                        string s;
                        getline(fin_log, s);
                        fout_log << s << endl;
                    }
                }

            }
            exit(0);
        }
        else {
            char *end_ptr;
            iterations_number = strtol(argv[1], &end_ptr, 10);
            if (*end_ptr || errno == ERANGE) {
                cout << "Incorrect parameter\n";
                exit(1);
            }
            cout << iterations_number << " iterations\n\n";
        }
    }

    get_input_data();

    double u_init_guess;
    if (init_guess_type == INIT_GUESS_UMIN)
        u_init_guess = umin;
    else //init_guess_type == INIT_GUESS_UMAX
        u_init_guess = umax;
    ns = N/2;
    neq = 2*(N+1)*(N+1)*(N+1);
    neq_phi = (N+1)*(N+1)*(N+1);

    Matrix A_phi(neq_phi, neq_phi), A(neq, neq);
    dense1D<double> rhs_phi(neq_phi), rhs(neq);
    dense1D<double> x_phi(neq_phi), x(neq);

    ::A_phi = &A_phi;
    ::A = &A;
    ::rhs_phi = &rhs_phi;
    ::rhs = &rhs;

    theta = new_double_array4d(tnum+1, N+1);
    phi = new_double_array4d(tnum+1, N+1);

    u = new_double_array3d(N+1);
    phiinit = new_double_array3d(N+1);
    thetaold = new_double_array3d(N+1);
    thetanew = new_double_array3d(N+1);
    phinew = new_double_array3d(N+1);
    p1 = new_double_array3d(N+1);
    p2 = new_double_array3d(N+1);
    p1old = new_double_array3d(N+1);
    p2old = new_double_array3d(N+1);
    gr = new_double_array3d(N+1);

    vis_u = new_bool_array2d(N+1);

    double h = ll/N;
    double tau = tt/tnum;
    double phib = pow(thetab, 4);

    //считаем количество файлов
    int prev_iterations = get_prev_iterations();

    int last_iteration = prev_iterations + iterations_number;
    for (int iteration_index = prev_iterations + 1;
        iteration_index <= last_iteration; ++iteration_index)
    {

        int kn; //в скольких узлах изменяем управление

        if (iteration_index == 1) {
            for (int i = 0; i <= ns; ++i) {
                for (int j = 0; j <= i; ++j)
                    u[i][j][0] = u_init_guess;
            }
            fill_control();
            kn = (ns + 1)*(ns + 2)/2;
        }
        else {
            read_control(iteration_index - 1);
            ifstream fin_k(k_file_name(iteration_index - 1).c_str());
            fin_k >> kn;

            //проверяем условие окончания итераций
            int crit;
            ifstream fin_crit(crit_file_name(iteration_index - 1).c_str());
            fin_crit >> crit;
            if (crit == 0) {
                cout << "No critical nodes\n";
                getch();
                exit(0);
            }

            //проверяем зацикливание
            int the_same_iter = 0;
            for (int iter = iteration_index - 2; iter >= 1; --iter) {
                if (the_same_control(iter)) {
                    the_same_iter = iter;
                    break;
                }
            }
            if (the_same_iter != 0) {
                cout << "The same control at iterations " << the_same_iter <<
                    " and " << iteration_index - 1 << "\n";
                getch();
                exit(0);
            }
        }

        //вычисление начального условия для интенсивности излучения
        cout << "(" << iteration_index << "/" << last_iteration << ") phi0\n";

        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    int borders = count_borders(i, j, k);
                    int i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, i5, j5, k5;
                    switch (borders) {
                    case 0:
                        fill_mat_phi(i, j, k, i-1, j, k,
                            -alpha/(h*h));
                        fill_mat_phi(i, j, k, i, j-1, k,
                            -alpha/(h*h));
                        fill_mat_phi(i, j, k, i, j, k-1,
                            -alpha/(h*h));
                        fill_mat_phi(i, j, k, i, j, k,
                            2*alpha/(h*h) + 2*alpha/(h*h) + 2*alpha/(h*h) + kappaa);
                        fill_mat_phi(i, j, k, i, j, k+1,
                            -alpha/(h*h));
                        fill_mat_phi(i, j, k, i, j+1, k,
                            -alpha/(h*h));
                        fill_mat_phi(i, j, k, i+1, j, k,
                            -alpha/(h*h));
                        fill_rhs_phi(i, j, k,
                            kappaa*pow(thetainit, 4));
                        break;
                    case 1:
                        get_ij_border(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                            i4, j4, k4, i5, j5, k5);
                        fill_mat_phi(i, j, k, i1, j1, k1,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i2, j2, k2,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i3, j3, k3,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i4, j4, k4,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i5, j5, k5,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i, j, k,
                            alpha/h + u[i][j][k] + alpha/h + alpha/h
                            + kappaa*h/2);
                        fill_rhs_phi(i, j, k,
                            u[i][j][k]*phib + kappaa*h/2*pow(thetainit, 4));
                        break;
                    case 2:
                        get_ij_edge(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                            i4, j4, k4);
                        fill_mat_phi(i, j, k, i1, j1, k1,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i2, j2, k2,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i3, j3, k3,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i4, j4, k4,
                            -alpha/(2*h));
                        fill_mat_phi(i, j, k, i, j, k,
                            alpha/h + u[i][j][k] + alpha/h + u[i][j][k] + alpha/h
                            + kappaa*h/2);
                        fill_rhs_phi(i, j, k,
                            u[i][j][k]*phib + u[i][j][k]*phib
                            + kappaa*h/2*pow(thetainit, 4));
                        break;
                    case 3:
                        get_ij_corner(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3);
                        fill_mat_phi(i, j, k, i1, j1, k1,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i2, j2, k2,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i3, j3, k3,
                            -alpha/h);
                        fill_mat_phi(i, j, k, i, j, k,
                            alpha/h + u[i][j][k] + alpha/h + u[i][j][k]
                            + alpha/h + u[i][j][k] + kappaa*h/2);
                        fill_rhs_phi(i, j, k,
                            u[i][j][k]*phib + u[i][j][k]*phib + u[i][j][k]*phib
                            + kappaa*h/2*pow(thetainit, 4));
                        break;
                    }
                }
            }
        }

        {
            //решение СЛАУ
            for (int s = 0; s < neq_phi; ++s)
                x_phi[s] = 0;
            ILU<Matrix> precond(A_phi);
            basic_iteration<double> iter(rhs_phi, max_iter, iter_rtol);
            bicgstab(A_phi, x_phi, rhs_phi, precond(), iter);
            if (iter.error_code()) {
                cout << "\n not converge!\n";
                exit(1);
            }
        }

        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    int s = get_index(THETA, i, j, k);
                    phiinit[i][j][k] = x_phi[s];
                }
            }
        }

        //начальные условия
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    theta[0][i][j][k] = thetainit;
                    phi[0][i][j][k] = phiinit[i][j][k];
                }
            }
        }

        //начальное приближение для итерационного метода решения СЛАУ
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    int s = get_index(THETA, i, j, k);
                    x[s] = theta[0][i][j][k];
                    s = get_index(PHI, i, j, k);
                    x[s] = phi[0][i][j][k];
                }
            }
        }

        for (int m = 1; m <= tnum; ++m) {
            cout << "(" << iteration_index << "/" << last_iteration << ") m = " << m << "  :  ";

            //начальное приближение для метода Ньютона
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k)
                        thetaold[i][j][k] = theta[m-1][i][j][k];
                }
            }

            for (int iter = 1; iter <= numiter_newton; ++iter) {
                cout << iter << " ";

                //заполнение матрицы СЛАУ и вектора правых частей

                for (int i = 0; i <= N; ++i) {
                    for (int j = 0; j <= N; ++j) {
                        for (int k = 0; k <= N; ++k) {
                            int borders = count_borders(i, j, k);
                            int i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, i5, j5, k5;
                            switch (borders) {
                            case 0:
                                fill_mat(THETA, i, j, k, THETA, i-1, j, k,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, THETA, i, j-1, k,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, THETA, i, j, k-1,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, THETA, i, j, k,
                                    1/tau + 0.5*( 2*a/(h*h) + 2*a/(h*h) + 2*a/(h*h)
                                        + b*kappaa*4*pow(thetaold[i][j][k], 3) ));
                                fill_mat(THETA, i, j, k, THETA, i, j, k+1,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, THETA, i, j+1, k,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, THETA, i+1, j, k,
                                    -0.5*a/(h*h));
                                fill_mat(THETA, i, j, k, PHI, i, j, k,
                                    -0.5*b*kappaa);
                                fill_rhs(THETA, i, j, k,
                                    theta[m-1][i][j][k]/tau + 0.5*b*kappaa*3*pow(thetaold[i][j][k], 4)
                                    - 0.5*( -a/(h*h)*(theta[m-1][i-1][j][k]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i+1][j][k])
                                        - a/(h*h)*(theta[m-1][i][j-1][k]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i][j+1][k])
                                        - a/(h*h)*(theta[m-1][i][j][k-1]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i][j][k+1])
                                        + b*kappaa*(pow(theta[m-1][i][j][k], 4) - phi[m-1][i][j][k]) ));

                                fill_mat(PHI, i, j, k, THETA, i, j, k,
                                    -kappaa*4*pow(thetaold[i][j][k], 3));
                                fill_mat(PHI, i, j, k, PHI, i-1, j, k,
                                    -alpha/(h*h));
                                fill_mat(PHI, i, j, k, PHI, i, j-1, k,
                                    -alpha/(h*h));
                                fill_mat(PHI, i, j, k, PHI, i, j, k-1,
                                    -alpha/(h*h));
                                fill_mat(PHI, i, j, k, PHI, i, j, k,
                                    2*alpha/(h*h) + 2*alpha/(h*h) + 2*alpha/(h*h) + kappaa);
                                fill_mat(PHI, i, j, k, PHI, i, j, k+1,
                                    -alpha/(h*h));
                                fill_mat(PHI, i, j, k, PHI, i, j+1, k,
                                    -alpha/(h*h));
                                fill_mat(PHI, i, j, k, PHI, i+1, j, k,
                                    -alpha/(h*h));
                                fill_rhs(PHI, i, j, k,
                                    -kappaa*3*pow(thetaold[i][j][k], 4));

                                break;
                            case 1:
                                get_ij_border(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                    i4, j4, k4, i5, j5, k5);

                                fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i4, j4, k4,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i5, j5, k5,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i, j, k,
                                    0.5*( a/h + beta + a/h + a/h +
                                        b*kappaa*h/2*4*pow(thetaold[i][j][k], 3) )
                                    + h/(2*tau));
                                fill_mat(THETA, i, j, k, PHI, i, j, k,
                                    -0.5*b*kappaa*h/2);
                                fill_rhs(THETA, i, j, k,
                                    0.5*( beta*thetab
                                        + b*kappaa*h/2*3*pow(thetaold[i][j][k], 4) )
                                    - 0.5*( -a/h*(theta[m-1][i1][j1][k1] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        - a/(2*h)*(theta[m-1][i2][j2][k2]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i3][j3][k3])
                                        - a/(2*h)*(theta[m-1][i4][j4][k4]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i5][j5][k5])
                                        + b*kappaa*h/2*(pow(theta[m-1][i][j][k], 4) - phi[m-1][i][j][k]) )
                                    + h/(2*tau)*theta[m-1][i][j][k]);

                                fill_mat(PHI, i, j, k, THETA, i, j, k,
                                    -kappaa*h/2*4*pow(thetaold[i][j][k], 3));
                                fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i4, j4, k4,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i5, j5, k5,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i, j, k,
                                    alpha/h + u[i][j][k] + alpha/h + alpha/h
                                    + kappaa*h/2);
                                fill_rhs(PHI, i, j, k,
                                    u[i][j][k]*phib - kappaa*h/2*3*pow(thetaold[i][j][k], 4));

                                break;
                            case 2:
                                get_ij_edge(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                    i4, j4, k4);

                                fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i4, j4, k4,
                                    -0.5*a/(2*h));
                                fill_mat(THETA, i, j, k, THETA, i, j, k,
                                    0.5*( a/h + beta + a/h + beta + a/h
                                        + b*kappaa*h/2*4*pow(thetaold[i][j][k], 3) )
                                    + h/(2*tau));
                                fill_mat(THETA, i, j, k, PHI, i, j, k,
                                    -0.5*b*kappaa*h/2);
                                fill_rhs(THETA, i, j, k,
                                    0.5*( beta*thetab + beta*thetab
                                        + b*kappaa*h/2*3*pow(thetaold[i][j][k], 4) )
                                    - 0.5*( -a/h*(theta[m-1][i1][j1][k1] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        - a/h*(theta[m-1][i2][j2][k2] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        - a/(2*h)*(theta[m-1][i3][j3][k3]
                                            - 2*theta[m-1][i][j][k] + theta[m-1][i4][j4][k4])
                                        + b*kappaa*h/2*(pow(theta[m-1][i][j][k], 4) - phi[m-1][i][j][k]) )
                                    + h/(2*tau)*theta[m-1][i][j][k]);

                                fill_mat(PHI, i, j, k, THETA, i, j, k,
                                    -kappaa*h/2*4*pow(thetaold[i][j][k], 3));
                                fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i4, j4, k4,
                                    -alpha/(2*h));
                                fill_mat(PHI, i, j, k, PHI, i, j, k,
                                    alpha/h + u[i][j][k] + alpha/h + u[i][j][k] + alpha/h
                                    + kappaa*h/2);
                                fill_rhs(PHI, i, j, k,
                                    u[i][j][k]*phib + u[i][j][k]*phib
                                    - kappaa*h/2*3*pow(thetaold[i][j][k], 4));

                                break;
                            case 3:
                                get_ij_corner(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3);

                                fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                    -0.5*a/h);
                                fill_mat(THETA, i, j, k, THETA, i, j, k,
                                    0.5*( a/h + beta + a/h + beta + a/h + beta
                                        + b*kappaa*h/2*4*pow(thetaold[i][j][k], 3) )
                                    + h/(2*tau));
                                fill_mat(THETA, i, j, k, PHI, i, j, k,
                                    -0.5*b*kappaa*h/2);
                                fill_rhs(THETA, i, j, k,
                                    0.5*( beta*thetab + beta*thetab + beta*thetab
                                        + b*kappaa*h/2*3*pow(thetaold[i][j][k], 4) )
                                    - 0.5*( -a/h*(theta[m-1][i1][j1][k1] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        - a/h*(theta[m-1][i2][j2][k2] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        - a/h*(theta[m-1][i3][j3][k3] - theta[m-1][i][j][k])
                                        + beta*(theta[m-1][i][j][k] - thetab)
                                        + b*kappaa*h/2*(pow(theta[m-1][i][j][k], 4) - phi[m-1][i][j][k]) )
                                    + h/(2*tau)*theta[m-1][i][j][k]);

                                fill_mat(PHI, i, j, k, THETA, i, j, k,
                                    -kappaa*h/2*4*pow(thetaold[i][j][k], 3));
                                fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                    -alpha/h);
                                fill_mat(PHI, i, j, k, PHI, i, j, k,
                                    alpha/h + u[i][j][k] + alpha/h + u[i][j][k]
                                    + alpha/h + u[i][j][k] + kappaa*h/2);
                                fill_rhs(PHI, i, j, k,
                                    u[i][j][k]*phib + u[i][j][k]*phib + u[i][j][k]*phib
                                    - kappaa*h/2*3*pow(thetaold[i][j][k], 4));

                                break;
                            }
                        }
                    }
                }

                {
                    //решение СЛАУ
                    ILU<Matrix> precond(A);
                    basic_iteration<double> iter(rhs, max_iter, iter_rtol);
                    bicgstab(A, x, rhs, precond(), iter);
                    if (iter.error_code()) {
                        cout << "\n not converge!\n";
                        exit(1);
                    }
                }

                for (int i = 0; i <= N; ++i) {
                    for (int j = 0; j <= N; ++j) {
                        for (int k = 0; k <= N; ++k) {
                            int s = get_index(THETA, i, j, k);
                            thetanew[i][j][k] = x[s];
                            s = get_index(PHI, i, j, k);
                            phinew[i][j][k] = x[s];
                        }
                    }
                }

                double maxdiff = 0.0;
                for (int i = 0; i <= N; ++i) {
                    for (int j = 0; j <= N; ++j) {
                        for (int k = 0; k <= N; ++k)
                            maxdiff = fmax(maxdiff, fabs(thetanew[i][j][k] - thetaold[i][j][k]));
                    }
                }
                if (maxdiff < eps_newton)
                    break;
                for (int i = 0; i <= N; ++i) {
                    for (int j = 0; j <= N; ++j) {
                        for (int k = 0; k <= N; ++k)
                            thetaold[i][j][k] = thetanew[i][j][k];
                    }
                }
            } // for iter

            cout << "\n";
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        theta[m][i][j][k] = thetanew[i][j][k];
                        phi[m][i][j][k] = phinew[i][j][k];
                    }
                }
            }

        } // for m


        //вычисление терминального условия
        if (cost_func == COST_FUNC_J1) {
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        p1old[i][j][k] = 0;
                        p2old[i][j][k] = 0;
                    }
                }
            }
        }
        else {
            cout << "(" << iteration_index << "/" << last_iteration << ") p2(T)\n";

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k)
                        p1old[i][j][k] = theta[tnum][i][j][k] - thetad;
                }
            }

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        int borders = count_borders(i, j, k);
                        int i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, i5, j5, k5;
                        switch (borders) {
                        case 0:
                            fill_mat_phi(i, j, k, i-1, j, k,
                                -alpha/(h*h));
                            fill_mat_phi(i, j, k, i, j-1, k,
                                -alpha/(h*h));
                            fill_mat_phi(i, j, k, i, j, k-1,
                                -alpha/(h*h));
                            fill_mat_phi(i, j, k, i, j, k,
                                2*alpha/(h*h) + 2*alpha/(h*h) + 2*alpha/(h*h) + kappaa);
                            fill_mat_phi(i, j, k, i, j, k+1,
                                -alpha/(h*h));
                            fill_mat_phi(i, j, k, i, j+1, k,
                                -alpha/(h*h));
                            fill_mat_phi(i, j, k, i+1, j, k,
                                -alpha/(h*h));
                            fill_rhs_phi(i, j, k,
                                kappaa*b*p1old[i][j][k]);
                            break;
                        case 1:
                            get_ij_border(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                i4, j4, k4, i5, j5, k5);
                            fill_mat_phi(i, j, k, i1, j1, k1,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i2, j2, k2,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i3, j3, k3,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i4, j4, k4,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i5, j5, k5,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + alpha/h
                                + kappaa*h/2);
                            fill_rhs_phi(i, j, k,
                                kappaa*h/2*b*p1old[i][j][k]);
                            break;
                        case 2:
                            get_ij_edge(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                i4, j4, k4);
                            fill_mat_phi(i, j, k, i1, j1, k1,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i2, j2, k2,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i3, j3, k3,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i4, j4, k4,
                                -alpha/(2*h));
                            fill_mat_phi(i, j, k, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + u[i][j][k] + alpha/h
                                + kappaa*h/2);
                            fill_rhs_phi(i, j, k,
                                kappaa*h/2*b*p1old[i][j][k]);
                            break;
                        case 3:
                            get_ij_corner(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3);
                            fill_mat_phi(i, j, k, i1, j1, k1,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i2, j2, k2,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i3, j3, k3,
                                -alpha/h);
                            fill_mat_phi(i, j, k, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + u[i][j][k]
                                + alpha/h + u[i][j][k] + kappaa*h/2);
                            fill_rhs_phi(i, j, k,
                                kappaa*h/2*b*p1old[i][j][k]);
                            break;
                        }
                    }
                }
            }

            {
                //решение СЛАУ
                for (int s = 0; s < neq_phi; ++s)
                    x_phi[s] = 0;
                ILU<Matrix> precond(A_phi);
                basic_iteration<double> iter(rhs_phi, max_iter, iter_rtol);
                bicgstab(A_phi, x_phi, rhs_phi, precond(), iter);
                if (iter.error_code()) {
                    cout << "\n not converge!\n";
                    exit(1);
                }
            }

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        int s = get_index(THETA, i, j, k);
                        p2old[i][j][k] = x_phi[s];
                    }
                }
            }

        } //else


        //инициализация градиента (первое слагаемое)
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    if (i != 0 && i != N && j != 0 && j != N && k != 0 && k != N)
                        continue;
                    double term = (phi[tnum][i][j][k] - phib)*p2old[i][j][k];
                    term /= 2;
                    gr[i][j][k] = term*tau;
                }
            }
        }

        //начальное приближение для итерационного метода решения СЛАУ
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    int s = get_index(THETA, i, j, k);
                    x[s] = p1old[i][j][k];
                    s = get_index(PHI, i, j, k);
                    x[s] = p2old[i][j][k];
                }
            }
        }

        for (int m = tnum-1; m >= 0; --m) {
            cout << "(" << iteration_index << "/" << last_iteration << ") m = " << m << "\n";

            //заполнение матрицы СЛАУ и вектора правых частей

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        int borders = count_borders(i, j, k);
                        int i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, i5, j5, k5;
                        switch (borders) {
                        case 0:
                            fill_mat(THETA, i, j, k, THETA, i-1, j, k,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, THETA, i, j-1, k,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, THETA, i, j, k-1,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, THETA, i, j, k,
                                1/tau + 0.5*( 2*a/(h*h) + 2*a/(h*h) + 2*a/(h*h)
                                    + 4*kappaa*pow(theta[m][i][j][k], 3)*b ));
                            fill_mat(THETA, i, j, k, THETA, i, j, k+1,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, THETA, i, j+1, k,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, THETA, i+1, j, k,
                                -0.5*a/(h*h));
                            fill_mat(THETA, i, j, k, PHI, i, j, k,
                                -0.5*4*kappaa*pow(theta[m][i][j][k], 3));
                            fill_rhs(THETA, i, j, k,
                                p1old[i][j][k]/tau
                                - 0.5*( -a/(h*h)*(p1old[i-1][j][k]
                                        - 2*p1old[i][j][k] + p1old[i+1][j][k])
                                    - a/(h*h)*(p1old[i][j-1][k]
                                        - 2*p1old[i][j][k] + p1old[i][j+1][k])
                                    - a/(h*h)*(p1old[i][j][k-1]
                                        - 2*p1old[i][j][k] + p1old[i][j][k+1])
                                    + 4*kappaa*pow(theta[m+1][i][j][k], 3)*
                                        (b*p1old[i][j][k] - p2old[i][j][k]) )
                                + 0.5*adj_src(m, i, j, k)
                                + 0.5*adj_src(m+1, i, j, k));

                            fill_mat(PHI, i, j, k, THETA, i, j, k,
                                -kappaa*b);
                            fill_mat(PHI, i, j, k, PHI, i-1, j, k,
                                -alpha/(h*h));
                            fill_mat(PHI, i, j, k, PHI, i, j-1, k,
                                -alpha/(h*h));
                            fill_mat(PHI, i, j, k, PHI, i, j, k-1,
                                -alpha/(h*h));
                            fill_mat(PHI, i, j, k, PHI, i, j, k,
                                2*alpha/(h*h) + 2*alpha/(h*h) + 2*alpha/(h*h) + kappaa);
                            fill_mat(PHI, i, j, k, PHI, i, j, k+1,
                                -alpha/(h*h));
                            fill_mat(PHI, i, j, k, PHI, i, j+1, k,
                                -alpha/(h*h));
                            fill_mat(PHI, i, j, k, PHI, i+1, j, k,
                                -alpha/(h*h));
                            fill_rhs(PHI, i, j, k,
                                0);

                            break;
                        case 1:
                            get_ij_border(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                i4, j4, k4, i5, j5, k5);

                            fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i4, j4, k4,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i5, j5, k5,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i, j, k,
                                0.5*( a/h + beta + a/h + a/h
                                    + 4*kappaa*h/2*pow(theta[m][i][j][k], 3)*b )
                                + h/(2*tau));
                            fill_mat(THETA, i, j, k, PHI, i, j, k,
                                -0.5*4*kappaa*h/2*pow(theta[m][i][j][k], 3));
                            fill_rhs(THETA, i, j, k,
                                h/(2*tau)*p1old[i][j][k]
                                - 0.5*( -a/h*(p1old[i1][j1][k1] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    - a/(2*h)*(p1old[i2][j2][k2]
                                        - 2*p1old[i][j][k] + p1old[i3][j3][k3])
                                    - a/(2*h)*(p1old[i4][j4][k4]
                                        - 2*p1old[i][j][k] + p1old[i5][j5][k5])
                                    + 4*kappaa*h/2*pow(theta[m+1][i][j][k], 3)*
                                        (b*p1old[i][j][k] - p2old[i][j][k]) )
                                + 0.5*h/2*adj_src(m, i, j, k)
                                + 0.5*h/2*adj_src(m+1, i, j, k));

                            fill_mat(PHI, i, j, k, THETA, i, j, k,
                                -kappaa*h/2*b);
                            fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i4, j4, k4,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i5, j5, k5,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + alpha/h
                                + kappaa*h/2);
                            fill_rhs(PHI, i, j, k,
                                0);

                            break;
                        case 2:
                            get_ij_edge(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3,
                                i4, j4, k4);

                            fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i4, j4, k4,
                                -0.5*a/(2*h));
                            fill_mat(THETA, i, j, k, THETA, i, j, k,
                                0.5*( a/h + beta + a/h + beta + a/h
                                    + 4*kappaa*h/2*pow(theta[m][i][j][k], 3)*b )
                                + h/(2*tau));
                            fill_mat(THETA, i, j, k, PHI, i, j, k,
                                -0.5*4*kappaa*h/2*pow(theta[m][i][j][k], 3));
                            fill_rhs(THETA, i, j, k,
                                h/(2*tau)*p1old[i][j][k]
                                - 0.5*( -a/h*(p1old[i1][j1][k1] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    - a/h*(p1old[i2][j2][k2] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    - a/(2*h)*(p1old[i3][j3][k3]
                                        - 2*p1old[i][j][k] + p1old[i4][j4][k4])
                                    + 4*kappaa*h/2*pow(theta[m+1][i][j][k], 3)*
                                        (b*p1old[i][j][k] - p2old[i][j][k]) )
                                + 0.5*h/2*adj_src(m, i, j, k)
                                + 0.5*h/2*adj_src(m+1, i, j, k));

                            fill_mat(PHI, i, j, k, THETA, i, j, k,
                                -kappaa*h/2*b);
                            fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i4, j4, k4,
                                -alpha/(2*h));
                            fill_mat(PHI, i, j, k, PHI, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + u[i][j][k] + alpha/h
                                + kappaa*h/2);
                            fill_rhs(PHI, i, j, k,
                                0);

                            break;
                        case 3:
                            get_ij_corner(i, j, k, i1, j1, k1, i2, j2, k2, i3, j3, k3);

                            fill_mat(THETA, i, j, k, THETA, i1, j1, k1,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i2, j2, k2,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i3, j3, k3,
                                -0.5*a/h);
                            fill_mat(THETA, i, j, k, THETA, i, j, k,
                                0.5*( a/h + beta + a/h + beta + a/h + beta
                                    + 4*kappaa*h/2*pow(theta[m][i][j][k], 3)*b )
                                + h/(2*tau));
                            fill_mat(THETA, i, j, k, PHI, i, j, k,
                                -0.5*4*kappaa*h/2*pow(theta[m][i][j][k], 3));
                            fill_rhs(THETA, i, j, k,
                                h/(2*tau)*p1old[i][j][k]
                                - 0.5*( -a/h*(p1old[i1][j1][k1] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    - a/h*(p1old[i2][j2][k2] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    - a/h*(p1old[i3][j3][k3] - p1old[i][j][k])
                                    + beta*p1old[i][j][k]
                                    + 4*kappaa*h/2*pow(theta[m+1][i][j][k], 3)*
                                        (b*p1old[i][j][k] - p2old[i][j][k]) )
                                + 0.5*h/2*adj_src(m, i, j, k)
                                + 0.5*h/2*adj_src(m+1, i, j, k));

                            fill_mat(PHI, i, j, k, THETA, i, j, k,
                                -kappaa*h/2*b);
                            fill_mat(PHI, i, j, k, PHI, i1, j1, k1,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i2, j2, k2,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i3, j3, k3,
                                -alpha/h);
                            fill_mat(PHI, i, j, k, PHI, i, j, k,
                                alpha/h + u[i][j][k] + alpha/h + u[i][j][k]
                                + alpha/h + u[i][j][k] + kappaa*h/2);
                            fill_rhs(PHI, i, j, k,
                                0);

                            break;
                        }
                    }
                }
            }

            {
                //решение СЛАУ
                ILU<Matrix> precond(A);
                basic_iteration<double> iter(rhs, max_iter, iter_rtol);
                bicgstab(A, x, rhs, precond(), iter);
                if (iter.error_code()) {
                    cout << "\n not converge!\n";
                    exit(1);
                }
            }

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        int s = get_index(THETA, i, j, k);
                        p1[i][j][k] = x[s];
                        s = get_index(PHI, i, j, k);
                        p2[i][j][k] = x[s];
                    }
                }
            }

            //прибавляем слагаемое к градиенту
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        if (i != 0 && i != N && j != 0 && j != N && k != 0 && k != N)
                            continue;
                        double term = (phi[m][i][j][k] - phib)*p2[i][j][k];
                        if (m == 0)
                            term /= 2;
                        gr[i][j][k] += term*tau;
                    }
                }
            }

            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    for (int k = 0; k <= N; ++k) {
                        p1old[i][j][k] = p1[i][j][k];
                        p2old[i][j][k] = p2[i][j][k];
                    }
                }
            }

        } // for m

        //вычисляем значение функционала
        double jj = 0;
        if (cost_func == COST_FUNC_J1) {
            for (int m = 0; m <= tnum; ++m) {
                double s = 0;
                for (int i = 0; i <= N; ++i) {
                    double term = 0;
                    for (int j = 0; j <= N; ++j) {
                        double term2 = 0;
                        for (int k = 0; k <= N; ++k) {
                            double term3 = pow(theta[m][i][j][k] - thetad, 2);
                            if (k == 0 || k == N)
                                term3 /= 2;
                            term2 += term3;
                        }
                        if (j == 0 || j == N)
                            term2 /= 2;
                        term += term2;
                    }
                    if (i == 0 || i == N)
                        term /= 2;
                    s += term;
                }
                if (m == 0 || m == tnum)
                    s /= 2;
                jj += s;
            }
            jj *= tau*h*h*h;
        }
        else { //cost_func == COST_FUNC_J2
            for (int i = 0; i <= N; ++i) {
                double term = 0;
                for (int j = 0; j <= N; ++j) {
                    double term2 = 0;
                    for (int k = 0; k <= N; ++k) {
                        double term3 = pow(theta[tnum][i][j][k] - thetad, 2);
                        if (k == 0 || k == N)
                            term3 /= 2;
                        term2 += term3;
                    }
                    if (j == 0 || j == N)
                        term2 /= 2;
                    term += term2;
                }
                if (i == 0 || i == N)
                    term /= 2;
                jj += term;
            }
            jj *= h*h*h;
        }
        ofstream fout_jj(jj_file_name(iteration_index).c_str());
        fout_jj.precision(17);
        fout_jj << jj;

        //вычисляем количество критических узлов
        int crit = 0;
        for (int i = 0; i <= ns; ++i) {
            for (int j = 0; j <= i; ++j) {
                if ((u[i][j][0] == umin && gr[i][j][0] > 0) ||
                        (u[i][j][0] == umax && gr[i][j][0] < 0))
                {
                    ++crit;
                }
            }
        }
        ofstream fout_crit(crit_file_name(iteration_index).c_str());
        fout_crit << crit;

        //выводим решение при z = L/2 в конечный момент времени
        ofstream fout_theta_phi(theta_phi_file_name(iteration_index).c_str());
        fout_theta_phi.precision(10);
        int kmid = N/2;
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j)
                fout_theta_phi << theta[tnum][i][j][kmid] << "   ";
            fout_theta_phi << "\n";
        }
        fout_theta_phi << "\n\n\n\n";
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j)
                fout_theta_phi << phi[tnum][i][j][kmid] << "   ";
            fout_theta_phi << "\n";
        }

        ofstream fout_log(log_file_name(iteration_index).c_str());
        fout_log << "iteration " << iteration_index << "\n";
        fout_log << "------------\n";

        //определяем, нужно ли уменьшить k
        if (kn > 1 && iteration_index > 1) {
            double min_jj;
            int min_iter;
            for (int iter = 1; iter < iteration_index; ++iter) {
                double jj1;
                ifstream fin_jj(jj_file_name(iter).c_str());
                fin_jj >> jj1;
                if (iter == 1 || jj1 < min_jj) {
                    min_jj = jj1;
                    min_iter = iter;
                }
            }
            if (jj >= min_jj) {
                //откатываем управление
                fout_log << "Returning to control at iteration " << min_iter - 1 << "\n";
                if (min_iter > 1)
                    read_control(min_iter - 1);
                else {  // min_iter == 1
                    for (int i = 0; i <= ns; ++i) {
                        for (int j = 0; j <= i; ++j)
                            u[i][j][0] = u_init_guess;
                    }
                    fill_control();
                }
                read_grad(min_iter); //считываем старый градиент
                //уменьшаем k
                kn /= 2;
                fout_log << "Reducing k to " << kn << "\n";
            }
        }

        //выводим k
        ofstream fout_k(k_file_name(iteration_index).c_str());
        fout_k << kn;

        //изменяем управление
        for (int i = 0; i <= ns; ++i) {
            for (int j = 0; j <= i; ++j)
                vis_u[i][j] = false;
        }
        for (int change_iter = 1; change_iter <= kn; ++change_iter) {
            double max_grad = 0;
            int ii = -1, jj = -1;
            for (int i = 0; i <= ns; ++i) {
                for (int j = 0; j <= i; ++j) {
                    if (!vis_u[i][j] &&
                        ((u[i][j][0] == umin && gr[i][j][0] > 0) ||
                            (u[i][j][0] == umax && gr[i][j][0] < 0)) &&
                        fabs(gr[i][j][0]) > max_grad)
                    {
                        max_grad = fabs(gr[i][j][0]);
                        ii = i;
                        jj = j;
                    }
                }
            }
            if (ii == -1) //не осталось критических узлов
                break;
            vis_u[ii][jj] = true;
            if (u[ii][jj][0] == umin)
                u[ii][jj][0] = umax;
            else
                u[ii][jj][0] = umin;
        }
        fill_control();

        //выводим управление
        ofstream fout(control_file_name(iteration_index).c_str());
        fout.precision(17);
        for (int i = 0; i <= ns; ++i) {
            for (int j = 0; j <= i; ++j) {
                if (vis_u[i][j])
                    fout << "   ";
                fout << u[i][j][0] << "   " << gr[i][j][0] << "\n";
            }
            fout << "\n";
        }

        //проверяем условие окончания итераций
        if (crit == 0) {
            fout_log << "No critical nodes\n";
            cout << "No critical nodes\n";
            break;
        }

        //проверяем зацикливание
        int the_same_iter = 0;
        for (int iter = iteration_index - 1; iter >= 1; --iter) {
            if (the_same_control(iter)) {
                the_same_iter = iter;
                break;
            }
        }
        if (the_same_iter != 0) {
            fout_log << "The same control at iterations " << the_same_iter <<
                " and " << iteration_index << "\n";
            cout << "The same control at iterations " << the_same_iter <<
                " and " << iteration_index << "\n";
            break;
        }

    } // for iteration_index

    delete_double_array4d(theta, tnum+1, N+1);
    delete_double_array4d(phi, tnum+1, N+1);

    delete_double_array3d(u, N+1);
    delete_double_array3d(phiinit, N+1);
    delete_double_array3d(thetaold, N+1);
    delete_double_array3d(thetanew, N+1);
    delete_double_array3d(phinew, N+1);
    delete_double_array3d(p1, N+1);
    delete_double_array3d(p2, N+1);
    delete_double_array3d(p1old, N+1);
    delete_double_array3d(p2old, N+1);
    delete_double_array3d(gr, N+1);

    delete_bool_array2d(vis_u, N+1);

    std::cout << "\nOK";
    getch();

    return 0;
}
