#include <omp.h>
#include <iostream>
using namespace std;


const int n = 2; //кол-во уравнений в системе

double F1(double y1, double y2) {
    return 2 * (y1 - y1 * y2);
}

double F2(double y1, double y2) {
    return -1 * (y2 - y1 * y2);
}

//определяем ф-ию вычисления правых частей
double func(double* y, double time, int i) {
	switch (i) {
    case 0: return F1(y[0], y[1]); break;
    case 1: return F2(y[0], y[1]); break;
	}
}

void Euler(double Y[], double& TIME, double T_MAX, double TAU) {
    
    cout << "Метод Эйлера:" << endl;

    cout << "Значения:" << endl;

    double yy[n] = { 0.0 }; // вспомогательный массив

    double t_nachalo = omp_get_wtime();

    do 
    {
        for (int i = 0; i < n; i++)
            yy[i] = Y[i] + TAU * func(Y, TIME, i);

        for (int i = 0; i < n; i++) {
            Y[i] = yy[i];
            if (abs(TIME - round(TIME)) < 1e-10) {
                cout << "y1 = " << Y[0] << "; " << "y2 = " << Y[1] << " " << "t = " << TIME << endl;
            }
        }

        TIME += TAU;

    } while (TIME <= T_MAX);

    double tk = omp_get_wtime();

    double tt = tk - t_nachalo;

    cout << "Время вычисления: " << tt << endl;

    cout << "" << endl;

}

void RK_2(double Y[], double& TIME, double T_MAX, double TAU) {

    cout << "RK_2:" << endl;

    cout << "Значения:" << endl;

    double y1_2[n] = { 0.0 }, //значения параметров в средней точке по времени
        ff[n] = { 0.0 }; // значения правых частей в средней точке по времени
    
    double t_nachalo = omp_get_wtime();

    do
    {
        // 4.1 организуем цикл по числу параметров и вычисляем значения параметров в средней точке шага по времени
        for (int i = 0; i < n; i++)
            y1_2[i] = Y[i] + 0.5 * TAU * func(Y, TIME, i);

        // 4.2 организуем цикл по числу параметров и вычисляем значения правых частей в средней точке шага по времени
        for (int i = 0; i < n; i++)
            ff[i] = func(y1_2, TIME + 0.5 * TAU, i);

        //4.3 организуем цикл по числу параметров и вычисляем новые значения параметров
        for (int i = 0; i < n; i++) {
            Y[i] += TAU * ff[i];
            if (abs(TIME - round(TIME)) < 1e-10) {
                cout << "y1 = " << Y[0] << "; " << "y2 = " << Y[1] << " " << "t = " << TIME << endl;
            }
        }
 
        // возвращаемся в начальную точку, делаем полный шаг по времени, но в правую часть половинчатую
        //делаем полный шаг по времени и проверяем условия итерационного процесса

        TIME += TAU;

    } while (TIME <= T_MAX);

    double tk = omp_get_wtime();

    double tt = tk - t_nachalo;

    cout << "Время вычисления: " << tt << endl;

    cout << "" << endl;
}

void Predictor_Correktor(double Y[], double& TIME, double T_MAX, double TAU) {

    cout << "Предиктор корректор:" << endl;

    cout << "Значения:" << endl;

    double yy_pr[n] = { 0.0 }, //прогнозируемое значение параметров (грубые)
        ff_pr[n] = { 0.0 }, // правые части в k-ый момент времени
        f_prog[n] = { 0.0 }; // прогнозируемые правые части (грубые)
    
    double t_nachalo = omp_get_wtime();

    do
    {
        // этап 1: прогноз
        // 4.1 организуем цикл по числу параметров и вычисляем значения правых частей в k-ый момент времени
        for (int i = 0; i < n; i++)
            ff_pr[i] = func(Y, TIME, i);

        // 4.2 организуем цикл по числу параметров и вычисляем прогнозируемые значения параметров
        for (int i = 0; i < n; i++)
            yy_pr[i] = Y[i] + TAU * ff_pr[i];

        //4.3 организуем цикл по числу параметров и вычисляем прогнозируемые значения правых частей (грубые значения)
        for (int i = 0; i < n; i++)
            f_prog[i] = func(yy_pr, TIME + TAU, i);

        // этап 2: коррекция
        //4.4 организуем цикл по числу параметров и вычисляем новые значения параметров
        for (int i = 0; i < n; i++) {
            Y[i] += 0.5 * TAU * (ff_pr[i] + f_prog[i]);
            if (abs(TIME - round(TIME)) < 1e-10) {
                cout << "y1 = " << Y[0] << "; " << "y2 = " << Y[1] << " " << "t = " << TIME << endl;
            }
        }

        //делаем шаг по времени и проверяем условия итерационного процесса
        TIME += TAU;
    } while (TIME <= T_MAX);

    double tk = omp_get_wtime();

    double tt = tk - t_nachalo;


    cout << "Время вычисления: " << tt << endl;

    cout << "" << endl;
}

void RK_4(double Y[], double& TIME, double T_MAX, double TAU) {

    cout << "RK_4:" << endl;

    cout << "Значения:" << endl;

    const int m = 4; // размерность массива коэффициентов приближения для каждого параметра
    double yy_rk_4[n] = { 0.0 }, // значения параметров в смещенных точках, когда старые + 1/2 предыдущего
        r[m][n] = { 0.0 }; // для коэффициентов (при таком раскладе весь массив заполняется нулями)
    
    double t_nachalo = omp_get_wtime();

    do
    {

        // 4.1 организуем цикл по числу параметров и вычисляем коэффициенты первого приближения
        for (int i = 0; i < n; i++)
            r[0][i] = TAU * func(Y, TIME, i);
        // будет вычисляться 00, 01, 02, 03

        // 4.2 организуем цикл по времени и вычисляем значения параметров в смещённых точках
        for (int i = 0; i < n; i++)
            yy_rk_4[i] = Y[i] + 0.5 * r[0][i]; // подготовили значения параметров в смещённых точках

        // 4.3 организуем цикл по числу параметров и вычисляем коэффициенты второго приближения для всех параметров по расчётным формулам
        for (int i = 0; i < n; i++)
            r[1][i] = TAU * func(yy_rk_4, TIME + 0.5 * TAU, i);

        // 4.4 организуем цикл по числу параметров и вычисляем значения параметров в смещённых точках
        for (int i = 0; i < n; i++)
            yy_rk_4[i] = Y[i] + 0.5 * r[1][i];
        // подготовка к вычислению коэффициента третьего приближения

        // 4.5 организуем цикл по числу параметров и вычисляем коэффициенты третьего приближения для всех параметров по расчётным формулам
        for (int i = 0; i < n; i++)
            r[2][i] = TAU * func(yy_rk_4, TIME + 0.5 * TAU, i);

        // 4.6 организуем цикл по числу параметров и вычисляем значения параметров в смещённых точках
        for (int i = 0; i < n; i++)
            yy_rk_4[i] = Y[i] + r[2][i];
        // подготовка к вычислению коэффициента четвёртого приближения

        // 4.7 организуем цикл по числу параметров и вычисляем коэффициенты четвертого приближения для всех параметров по расчётным формулам
        for (int i = 0; i < n; i++)
            r[3][i] = TAU * func(yy_rk_4, TIME + TAU, i);

        //4.8 организуем цикл по числу параметров и вычисляем новые значения параметров
        for (int i = 0; i < n; i++) {
            Y[i] += (r[0][i] + 2.0 * r[1][i] + 2.0 * r[2][i] + r[3][i]) / 6.0;
            if (abs(TIME - round(TIME)) < 1e-10) {
                cout << "y1 = " << Y[0] << "; " << "y2 = " << Y[1] << " " << "t = " << TIME << endl;
            }
        }

        TIME += TAU;

    } while (TIME <= T_MAX);

    double tk = omp_get_wtime();

    double tt = tk - t_nachalo;

    cout << "Время вычисления: " << tt << endl;

    cout << "" << endl;

}

void Implicit_Euler(double Y[], double& TIME, double T_MAX, double TAU) {

    cout << "Неявный метод Эйлера:" << endl;

    cout << "Значения:" << endl;

    double h = 0.000001;
    double b[n] = { 0.0 },  //вспомогательный массив для хранения правых частей в k момент времени
        a[n][n] = { 0.0 }, //квадратичный массив для хранения основной матрицы
        p[n] = { 0.0 }; //массив приращений

    double t_nachalo = omp_get_wtime();



    do
    {
        //вычисляем правые части в текущей точке
        for (int i = 0; i < n; i++)
            b[i] = -func(Y, TIME, i);

        for (int i = 0; i < n; i++) { //иду по строчкам
            for (int j = 0; j < n; j++) { //иду по столбикам

                double Y_plus_h[n] = { 0.0 }; //вспомогательный массив

                for (int k = 0; k < n; k++)
                    Y_plus_h[k] = Y[k];
                Y_plus_h[j] += h; // добавляем h к j-й переменной

                a[i][j] = (func(Y_plus_h, TIME, i) - func(Y, TIME, i)) / h;  //заполнений главной матрицы
            }
        }

        for (int i = 0; i < n; i++)
            a[i][i] -= 1.0 / TAU; //корректируем главную диагональ

        double main_determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0]; //рассчет главного определителя

        if (!main_determinant) { //проверка на ненулевой определитель
            cout << "Определитель равен 0 :(" << endl;
            exit(1);
        }

        double determinant[n] = { 0.0 };

        determinant[0] = b[0] * a[1][1] - b[1] * a[0][1]; //рассчет побочных определителей
        determinant[1] = b[1] * a[0][0] - b[0] * a[1][0];

        for (int i = 0; i < n; i++)
            p[i] = determinant[i] / main_determinant;

        for (int i = 0; i < n; i++) {
            Y[i] += p[i];
            if (abs(TIME - round(TIME)) < 1e-10) {
                cout << "y1 = " << Y[0] << "; " << "y2 = " << Y[1] << " " << "t = " << TIME << endl;
            }
        }

        TIME += TAU;

    } while (TIME <= T_MAX);

    double tk = omp_get_wtime();

    double tt = tk - t_nachalo;

    cout << "Время вычисления: " << tt << endl;
}


int main() {

    setlocale(LC_ALL, "");

        {// метод Эйлера
            double y[n] = { 1.0, 3.0 },
                t_0 = 0.0, t_max = 10.0, tau = 0.001;
            double time = t_0;

            Euler(y, time, t_max, tau);
        }
    
        {// метод Рунге Кутто 2
            double y[n] = { 1.0, 3.0 },
                t_0 = 0.0, t_max = 10.0, tau = 0.001;
            double time = t_0;

            RK_2(y, time, t_max, tau);
        }

        { // метод предиктор-корректор
            double y[n] = { 1.0, 3.0 },
                t_0 = 0.0, t_max = 10.0, tau = 0.001;
            double time = t_0;

            Predictor_Correktor(y, time, t_max, tau);
        }

        {// метод Рунге Кутто 4
            double y[n] = { 1.0, 3.0 },
                t_0 = 0.0, t_max = 10.0, tau = 0.001;
            double time = t_0;

            RK_4(y, time, t_max, tau);
        }

        {// неявный метод Эйлера
            double y[n] = { 1.0, 3.0 },
                t_0 = 0.0, t_max = 10.0, tau = 0.001;
            double time = t_0;

            Implicit_Euler(y, time, t_max, tau);
        }
}


