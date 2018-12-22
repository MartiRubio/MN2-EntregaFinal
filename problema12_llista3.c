// Author: Martí Rubio
// Codi @ https://github.com/MartiRubio/MN2-EntregaParcial

/*
Codi per al problema 26 de la llista 1 de Mètodes Numèrics II
*/

// Includes necessaris
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


double grad_result_g[2];
double grad_result[2];
double start[2];
double jacobian[2][2];
double jacobian_calculation[2][2];
double newton_actual[2];
double newton_anterior[2];
double tol;


double f(double x, double y)
{
    double result = pow(x, 4.0);

    result += 2.0*pow(x,3.0)*y;

    result += 3.0*pow(x,2.0)*pow(y,2.0);

    result += 2.0*x*pow(y,3.0);

    result += 2.0*pow(y,4.0);

    result -= 2.6*pow(x,3.0);

    result -= 3.2*pow(x,2.0)*y;

    result -= 2.6*x*pow(y,2.0);

    result -= 3.2*pow(y,3.0);

    result -= 7.22*pow(x,2.0);

    result -= 16.0*x*y;

    result -= 15.22*pow(y,2.0);

    result += 20.8*x;

    result += 25.6*y;

    result -= 5.94;

    return result;
}


double* grad_f(double x, double y)
{
    grad_result[0] = 4*pow(x, 3.0);

    grad_result[0] += 6.0*pow(x,2.0)*y;

    grad_result[0] += 6.0*x*pow(y,2.0);

    grad_result[0] += 2.0*pow(y,3.0);

    grad_result[0] -= 2.6*3*pow(x,2.0);

    grad_result[0] -= 3.2*2*x*y;

    grad_result[0] -= 2.6*pow(y,2.0);

    grad_result[0] -= 7.22*2*x;

    grad_result[0] -= 16.0*y;

    grad_result[0] += 20.8;

    grad_result[1] = 2.0*pow(x,3.0);

    grad_result[1] += 6.0*pow(x,2.0)*y;

    grad_result[1] += 6.0*x*pow(y,2.0);

    grad_result[1] += 8.0*pow(y,3.0);

    grad_result[1] -= 3.2*pow(x,2.0);

    grad_result[1] -= 2.6*x*2*y;

    grad_result[1] -= 3.2*3*pow(y,2.0);

    grad_result[1] -= 16.0*x;

    grad_result[1] -= 15.22*2*y;

    grad_result[1] += 25.6;

    return grad_result;
}

double g(double x, double y, double x0, double y0, double h)
{
    x = (x - x0);
    y = (y - y0);

    double result = pow(x,2.) + pow(y,2.) - pow(h, 2.);
    return result;
}


double* grad_g(double x, double y, double x0, double y0)
{
    grad_result_g[0] = 2.0*(x - x0);

    grad_result_g[1] = 2.0*(y - y0);

    return grad_result_g;
}

double find_solution()
{
    tol = pow(10.0, -12.0);
    double y = 0.0;
    double ant_y = 1.0;
    double new_y;
    while (fabs(y - ant_y) > tol){
        new_y = (y + ant_y) / 2.0;
        if (f(0, new_y)*f(0, y) > 0){
            y = new_y;
        }
        else{
            ant_y = y;
            y = new_y;
        }
    }
    return y;
}

double norma(double* vector)
{
    double result;
    result = pow(vector[0], 2.0) + pow(vector[1], 2.0);
    result = sqrt(result);

    return result;
}

double* get_first_guess(double* solution, double h)
{
    double* vector;
    vector = grad_f(solution[0], solution[1]);
    double norm = norma(vector);
    vector[0] = -h*vector[1]/norm;
    vector[1] = h*vector[0]/norm;
    start[0] = newton_anterior[0] + vector[0];
    start[1] = newton_anterior[1] + vector[1];
    //printf("origin:(%.12f,%.12f)\n", solution[0], solution[1]);
    //printf("start:(%.12f,%.12f)\n", vector[0], vector[1]);

    return start;
}

void inverse_jacobian(double x, double y, double x0, double y0)
{
    double* gradf;
    double* gradg;
    gradf = grad_f(x, y);
    gradg = grad_g(x, y, x0, y0);
    double det = 1./(gradf[0]*gradg[1] - gradf[1]*gradg[0]);
    jacobian[0][0] = gradg[1]*det;
    jacobian[0][1] = -gradf[1]*det;
    jacobian[1][0] = -gradg[0]*det;
    jacobian[1][1] = gradf[0]*det;
}

double* newton_method(double x, double y, double h)
{
    double* fg;
    double* gradf;
    double* gradg;
    double f_calc, g_calc;
    double** ccc;
    newton_anterior[0] = x;
    newton_anterior[1] = y;
    fg = get_first_guess(newton_anterior,h);
    newton_anterior[0] = fg[0];
    newton_anterior[1] = fg[1];

    inverse_jacobian(newton_anterior[0], newton_anterior[1], x, y);

    f_calc = f(newton_anterior[0], newton_anterior[1]);
    g_calc = g(newton_anterior[0], newton_anterior[1], x, y, h);
    int i = 0;
    newton_actual[0] = newton_anterior[0] - jacobian[0][0]*f_calc - jacobian[0][1]*g_calc;
    newton_actual[1] = newton_anterior[1] - jacobian[1][0]*f_calc - jacobian[1][1]*g_calc;
    while (fabs(f(newton_actual[0], newton_actual[1])) > pow(10.,-10.)){
        i++;
        // printf("%.12f\n", fabs(f(newton_actual[0], newton_actual[1])));
        newton_anterior[0] = newton_actual[0];
        newton_anterior[1] = newton_actual[1];
        inverse_jacobian(newton_anterior[0], newton_anterior[1], x, y);
        f_calc = f(newton_anterior[0], newton_anterior[1]);
        g_calc = g(newton_anterior[0], newton_anterior[1], x, y, h);
        newton_actual[0] = newton_anterior[0] - jacobian[0][0]*f_calc - jacobian[0][1]*g_calc;
        newton_actual[1] = newton_anterior[1] - jacobian[1][0]*f_calc - jacobian[1][1]*g_calc;
    }
    //printf("%d\n", i);
    return newton_actual;
}

int main()
{
    // Trobem un punt de la funció
    double x = 0;
    double y;
    double h = 0.01;
    double* solucio;
    FILE *fp;

    y = find_solution();
    printf("Un punt de la corba és: (0,%f)\n", y);

    solucio = newton_method(x, y, h);
    // printf("%.12f\n", solucio[0]);
    // printf("%.12f\n", solucio[1]);
    fp = fopen("./test.txt", "w+");
    for(int i = 0; i < 10; i++){
        fprintf(fp, "%.12f, %.12f \n", solucio[0], solucio[1]);
        solucio = newton_method(solucio[0], solucio[1], h);
    }
    fclose(fp);


    FILE *pipe_gp = popen("gnuplot -p", "w");
    fputs("plot  './test.txt' u 0:1 '\n", pipe_gp);
    pclose(pipe_gp);
    return 0;
}

