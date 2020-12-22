#include <iostream>
using namespace std;

void Out_of_system_of_equations(double** Matrix_A, double* Vector_b, int rank)
{
    for (int i = 0; i < rank; i++)
    {
        for (int j = 0; j < rank; j++)
        {
            cout << Matrix_A[i][j] << "*x" << j;
            if (j < rank - 1)
                cout << " + ";
        }
        cout << " = " << Vector_b[i] << endl;
    }
    return;
}
double* Solution_by_the_Gaussian_method(double** Matrix_A, double* Vector_b, int rank)
{
    double* Vector_X, max;
    int k, index;
    const double eps = 0.000001;  // точность
    Vector_X = new double[rank];
    k = 0;
    while (k < rank)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(Matrix_A[k][k]);
        index = k;
        for (int i = k + 1; i < rank; i++)
        {
            if (abs(Matrix_A[i][k]) > max)
            {
                max = abs(Matrix_A[i][k]);
                index = i;
            }
        }        
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        // Перестановка строк
        for (int j = 0; j < rank; j++)
        {
            double temp = Matrix_A[k][j];
            Matrix_A[k][j] = Matrix_A[index][j];
            Matrix_A[index][j] = temp;
        }
        double temp = Vector_b[k];
        Vector_b[k] = Vector_b[index];
        Vector_b[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < rank; i++)
        {
            double temp = Matrix_A[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] / temp;
            Vector_b[i] = Vector_b[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] - Matrix_A[k][j];
            Vector_b[i] = Vector_b[i] - Vector_b[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = rank - 1; k >= 0; k--)
    {
        Vector_X[k] = Vector_b[k];
        for (int i = 0; i < k; i++)
            Vector_b[i] = Vector_b[i] - Matrix_A[i][k] * Vector_X[k];
    }
    return Vector_X;
    delete[]Vector_X;
}
double Error(double* Vector_X, double* Vector_solution_b, int rank)
{
    double max_ZNAM_V_x = abs(Vector_X[0]);// поиск максимального значения знаменателя формулы погрешности
    for (int i = 0; i < rank; i++)
    {
        if (abs(Vector_X[i]) > max_ZNAM_V_x)
        {
            max_ZNAM_V_x = abs(Vector_X[i]);
        }
    }
    double max_CHISL_V_x = abs(Vector_solution_b[0] - Vector_X[0]);// поиск максимального значения ЧИСЛИТЕЛЯ формулы погрешности
    for (int i = 0; i < rank; i++)
    {
        if (abs(Vector_solution_b[i] - Vector_X[i]) > max_CHISL_V_x)
        {
            max_CHISL_V_x = abs(Vector_solution_b[i] - Vector_X[i]);
        }
    }
    double ERROR;
    ERROR = max_CHISL_V_x / max_ZNAM_V_x;
    return ERROR;
}
int main()
{
    setlocale(LC_ALL, "rus");

    double** Matrix_A, * Vector_b, * Vector_X;
    int rank;
    
    cout << "Введите количество уравнений: ";
    cin >> rank;
    Matrix_A = new double* [rank];//deleted
    Vector_b = new double[rank];//deleted
    for (int i = 0; i < rank; i++)
    {
        Matrix_A[i] = new double[rank];//deleted
        for (int j = 0; j < rank; j++)
        {
            cout << "Matrix_A[" << i << "][" << j << "]= ";
            cin >> Matrix_A[i][j];
        }
    }    
    for (int i = 0; i < rank; i++)
    {
        cout << "Vector_b[" << i << "]= ";
        cin >> Vector_b[i];
    }
    cout << endl << endl << endl;
    cout << "Введённая матрица А:" << endl;
    for (int i = 0; i < rank; i++)
    {
        for (int j = 0; j < rank; j++)
        {
            cout << Matrix_A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Введённый вектор b:" << endl;
    for (int i = 0; i < rank; i++)
    {
        cout << Vector_b[i] << endl;        
    }
    double** Matrix_A_original, * Vector_b_original;

    //Сохранение исходных данных матрицы А и вектора b
    Matrix_A_original = new double* [rank];//deleted
    Vector_b_original = new double[rank];//deleted
    for (int i = 0; i < rank; i++)
    {
        Matrix_A_original[i] = new double[rank];//deleted
        for (int j = 0; j < rank; j++)
        {
            Matrix_A_original[i][j] = Matrix_A[i][j];
        }
    }
    for (int i = 0; i < rank; i++)
    {
        Vector_b_original[i] = Vector_b[i];
    }
    Out_of_system_of_equations(Matrix_A, Vector_b, rank);
    ////////////////////////////////////////////////////////////////////////
    Vector_X = Solution_by_the_Gaussian_method(Matrix_A, Vector_b, rank);
    for (int i = 0; i < rank; i++)
        cout << "x[" << i << "]=" << Vector_X[i] << endl;   
    //////////////////////////////////////////////////////////////////
    double* residual_vector;
    residual_vector = new double[rank];//deleted

    double max_RV = 0.000000;
    for (int i = 0; i < rank; i++)
    {
        residual_vector[i] = 0;
        for (int j = 0; j < rank; j++)
        {
            residual_vector[i] += Matrix_A_original[i][j] * Vector_X[j];
        }
        residual_vector[i] -= Vector_b_original[i];

        if (residual_vector[i] > max_RV)
        {
            max_RV = residual_vector[i];
        }
    }

    cout << "Вектор невязки :" << endl;
    cout << max_RV << endl;  

    /////////////////////////////////////////////////////////////////////
    double* Vector_solution_b = new double[rank];//deleted
    for (int i = 0; i < rank; i++)
    {
        Vector_solution_b[i] = 0;
        for (int j = 0; j < rank; j++)
        {
            Vector_solution_b[i] += Matrix_A_original[i][j] * Vector_X[j];
        }
    }
    cout << "новое решение: " << endl;
    for (int i = 0; i < rank; i++)
    {
        cout << Vector_solution_b[i] << endl;
    }
    ///////////////////////////////////////////////////////////////////////
    double* Vector_error_X = new double[rank];//deleted
    Vector_error_X = Solution_by_the_Gaussian_method(Matrix_A_original, Vector_solution_b, rank);
    ///////////////////////////////////////////////////////////////////////
    double error;
    error = Error(Vector_X, Vector_solution_b, rank);
    cout << "Погрешность :" << endl;
    cout << error << endl;
    //////////////////////////////////////////////////////////////////////
    for (int i = 0; i < rank; i++)
    {
        delete[]Matrix_A[i];        
    }
    delete[]Matrix_A;

    delete[]Vector_b;

    delete[]residual_vector;

    for (int i = 0; i < rank; i++)
    {
        delete[]Matrix_A_original[i];
    }
    delete[]Matrix_A_original;

    delete[]Vector_b_original;

    delete[]Vector_solution_b;

    delete[]Vector_error_X;

    return 0;
}