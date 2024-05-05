#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Необязательный флаг для дебага (можно удалить)
#define DEBUG 0

// Фиксированная размерность матрицы (обрабатывается препроцессором, просто alias)
#define N 10

// Тип вектора
typedef double vec[N];

// Тип матрицы
typedef double mat[N][N];

void print_mat(mat m);
// Скалярное произведение
double dot(vec v, vec u) {
    double res = 0;
    for(int i = 0; i < N; i++) res += v[i]*u[i];

    return res;
}

// Прибавление вектора
void add(vec v, vec u) {
    for(int i = 0; i < N; i++) {
        v[i] += u[i];
    }
}

// Вычитание вектора
void sub(vec v, vec u) {
    for(int i = 0; i < N; i++) {
        v[i] -= u[i];
    }
}

// Умножение на скаляр
void mul(vec v, double a) {
    for(int i = 0; i < N; i++) {
        v[i] *= a;
    }
}

// Умножение матрицы на вектор (res = m*v)
void apply(mat m, vec v, vec res) {   
    for(int i = 0; i < N; i++) {
        res[i] = 0;
        for(int j = 0; j < N; j++) {
            res[i] += m[i][j] * v[j]; 
        }
    }
}

// Умножение матриц (res = m*n)
void multiply(mat m, mat n, mat res) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            res[i][j] = 0;
            for(int k = 0; k < N; k++) {
                res[i][j] += m[i][k] * n[k][j];
            }
        }
    }
}

// Копирование матрицы (copy = m)
void copy_mat(mat m, mat copy) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            copy[i][j] = m[i][j];
        }
    }
}

// Получение j-ого столбца матрицы (res = m[*][j])
void get_column(mat m, int j, vec res) {
    for(int i = 0; i < N; i++) {
        res[i] = m[i][j];
    }
}

// Получение i-ой строки матрицы (res = m[i][*])
void get_row(mat m, int i, vec res) {
    for(int j = 0; j < N; j++) {
        res[j] = m[i][j];
    }
}

// Задание i-ого столбца матрицы (m[*][i] = col)
void set_column(mat m, int i, vec col) {
    for(int j = 0; j < N; j++) {
        m[j][i] = col[j];
    }
}

void normalize(mat m) {
    vec v;
    for(int i = 0; i < N; i++) {
        get_column(m, i, v);
        double norm = sqrt(dot(v,v));
        mul(v, 1/norm);
        set_column(m, i, v);
    }
}

// Ортогонализация Грамма-Шмидта
void ortogonalize(mat m) {
    vec v;
    vec u;
    for(int i = 0; i < N; i++) {
        get_column(m, i, v);
        for(int j = 0; j < i; j++) {
            get_column(m, j, u);

            double c = dot(v, u)/dot(u, u);
            mul(u, c);
            sub(v, u);
        }
        set_column(m, i, v);
    }
    normalize(m);
}

// QR-разложение (m = q*r)
void qr_decomposition(mat m, mat q, mat r) {
    copy_mat(m, q);
    ortogonalize(q);

    vec v;
    vec u;
    for(int j = 0; j < N; j++) {
        get_column(m, j, v);
        for(int i = 0; i <= j; i++) {
            get_column(q, i, u);
            r[i][j] = dot(v, u);
        }
        for(int i = j+1; i < N; i++) r[i][j] = 0;
    }
}

// Транспонирование матрицы
void transpose(mat m) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            double d = m[i][j];
            m[i][j] = m[j][i];
            m[j][i] = d;
        }
    }
}

// Обращение верхне-треугольной матрицы (inv = r^{-1}) 
void inverse_r(mat r, mat inv) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            inv[i][j] = i == j ? 1./r[i][i] : 0;
        }
    }

    vec v;
    vec u;
    for(int j = 0; j < N; j++) {
        for(int i = j-1; i >= 0; i--) {
            get_row(r, i, v);
            get_column(inv, j, u);
            inv[i][j] = -dot(v, u)/r[i][i];
        }
    }
}

// Обращение матрицы (inv = m^{-1})
void inverse(mat m, mat inv) {
    mat q;
    mat r;
    qr_decomposition(m, q, r);
    mat inv_r;
    inverse_r(r, inv_r);
    transpose(q);
    multiply(inv_r, q, inv);
}

// Решение СЛАУ методом QR-разложения
void solve(mat m, vec v, vec sol) {
    mat inv;
    inverse(m, inv);
    apply(inv, v, sol);
}

// Чтение матрицы из файла
void mat_from_file(char* filename, mat m) {
	double d;
	FILE* f;
	if ((f = fopen(filename, "r")) == NULL) return;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fscanf(f, "%lf,", &d);
			m[i][j] = d;
		}
	}
	fclose(f);
}

// Красивая печать матрицы
void print_mat(mat m) {
    for(int j = 0; j < N; j++) {
        printf("-----------|");
    }
    printf("\n");
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            printf("%10.4lf |", m[i][j]);
        }
        printf("\n");
        for(int j = 0; j < N; j++) {
            printf("-----------|");
        }
        printf("\n");
    }
}

// Красивая печать вектора
void print_vec(vec v) {
    for(int j = 0; j < N; j++) {
        printf("%10.4lf |", v[j]);
    }
    printf("\n");
}

int main(){
    // Матрица системы
    mat m;
    mat_from_file("lab_3.txt", m);
    printf("Coefficients: \n");
    print_mat(m);

    #if 1
    // Решение по умолчанию (x[i] = 2*i + 1)
    vec x;
    for(int i = 0; i < N; i++) {
        x[i] = 2*i+1;
    }

    // Свободный столбец уравнения (b = m*x)
    vec b;
    apply(m, x, b);
    printf("Constant terms: \n");
    print_vec(b);

    // Решение с помощью QR-разложения
    vec sol;
    solve(m, b, sol);
    printf("Solution: \n");
    print_vec(sol);
    #endif 

    // Тесты
    #if 0
    ortogonalize(m);
    print_mat(m);
    for(int i = 0; i < N; i++){
        for(int j = i; j < N; j++){
            vec v;
            vec u;
            get_column(m, i, v);
            get_column(m, j, u);
            printf("%d:%d - %lf\n", i, j, dot(u,v));
        }
    }
    #endif

    #if 0
    mat q;
    mat r;
    qr_decomposition(m, q, r);
    print_mat(m);
    print_mat(q);
    print_mat(r);

    mat inv;
    inverse_r(r, inv);
    print_mat(inv);
    #endif

    #if 0
    transpose(m);
    print_mat(m);
    #endif

    #if 0
    mat inv;
    inverse(m, inv);
    print_mat(inv);
    #endif
}