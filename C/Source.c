#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
/*Задача:
Реализовать структуру matrix, которая включает в себя следующие поля:
Количество строк
Количество столбцов
Двумерный массив
Максимальный размер структуры – 10x10.

При этом, реализовать следующие методы:

bool init_mtx( matrix* mtx, int n_rows, int n_cols);

bool sum_mtx( matrix* mtx1, matrix* mtx2, matrix* res);

bool product_mtx( matrix* mtx1, matrix* mtx2, matrix* res);

void print_mtx( matrix* mtx1, matrix* mtx2, matrix* res);

bool det_mtx( matrix* mtx, double* det );

bool inv_mtx( matrix* mtx, matrix* res );

 Каждая функция возвращает true или false для отработки ошибок.
 */

typedef struct matrix {
	int row;
	int col;
	double arr[10][10];
} matrix;


double generateGaussianNoise(double sigma) {
	double u1 = ((double)rand() / RAND_MAX);
	double u2 = ((double)rand() / RAND_MAX);

	double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14 * u2);

	return z0 * sigma;
}


bool clear_mtx(matrix* mtx) {
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			mtx->arr[i][j] = 0;
		}
	}
}


bool init_mtx(matrix* mtx, int n_rows, int n_cols)
{
	if ((n_rows <= 0) || (n_cols <= 0) || (n_rows > 10) || (n_cols > 10)) {
		return false;
	}
	mtx->col = n_cols;
	mtx->row = n_rows;
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			mtx->arr[i][j] = 0;
		}
	}
	return true;
}


void print_mtx(matrix* mtx) {
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			printf("%f", mtx->arr[i][j]);
		}
		printf("\n");
	}
}


bool sum_mtx(matrix* mtx1, matrix* mtx2, matrix* res)
{
	if ((mtx1->col != mtx2->col) || (mtx1->row != mtx2->row)) {
		return false;
	}
	res->col = mtx1->col;
	res->row = mtx1->row;
	clear_mtx(res);
	for (int i = 0; i < mtx1->row; i++) {
		for (int j = 0; j < mtx1->col; j++) {
			res->arr[i][j] = mtx1->arr[i][j] + mtx2->arr[i][j];
		}
	}
	return true;
}


bool diff_mtx(matrix* mtx1, matrix* mtx2, matrix* res)
{
	if ((mtx1->col != mtx2->col) || (mtx1->row != mtx2->row)) {
		return false;
	}
	res->col = mtx1->col;
	res->row = mtx1->row;
	clear_mtx(res);
	for (int i = 0; i < mtx1->row; i++) {
		for (int j = 0; j < mtx1->col; j++) {
			res->arr[i][j] = mtx1->arr[i][j] - mtx2->arr[i][j];
		}
	}
	return true;
}


bool product_mtx(matrix* mtx1, matrix* mtx2, matrix* res)
{

	if (mtx1->col != mtx2->row) {
		return false;
	}
	res->col = mtx2->col;
	res->row = mtx1->row;
	clear_mtx(res);
	for (int i = 0; i < res->row; i++) {
		for (int j = 0; j < res->col; j++) {
			double sum = 0;
			for (int k = 0; k < mtx1->col; k++) {
				sum += mtx1->arr[i][k] * mtx2->arr[k][j];
			}
			res->arr[i][j] = sum;

		}
	}
	return true;
}


void transp_mtx(matrix* mtx, matrix* res)
{
	res->col = mtx->row;
	res->row = mtx->col;
	clear_mtx(res);
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			res->arr[j][i] = mtx->arr[i][j];
		}
	}
}

bool det_mtx(matrix* mtx, double* det)
{
	matrix U;
	U.col = mtx->col;
	U.row = mtx->row;
	matrix L;
	L.col = mtx->col;
	L.row = mtx->row;
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			U.arr[i][j] = 0;
			L.arr[i][j] = 0;
		}
		L.arr[i][i] = 1;
	}
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			if (i <= j) {
				double sum_LU = 0;
				for (int k = 0; k < i; k++)
				{
					sum_LU += L.arr[i][k] * U.arr[k][j];
				}
				U.arr[i][j] = mtx->arr[i][j] - sum_LU;
			}
			else {
				double sum_LU = 0;
				for (int k = 0; k < i; k++)
				{
					sum_LU += L.arr[i][k] * U.arr[k][j];
				}
				L.arr[i][j] = (mtx->arr[i][j] - sum_LU) / U.arr[j][j];
			}
		}
	}
	double pr_L = 1;
	double pr_U = 1;
	for (int i = 0; i < mtx->row; i++) {
		pr_L *= L.arr[i][i];
		pr_U *= U.arr[i][i];
	}
	*det = pr_L * pr_U;
	return true;
}


bool inv_mtx(matrix* mtx, matrix* res)
{
	matrix mtx_t;
	transp_mtx(mtx, &mtx_t);
	double det = 0;
	det_mtx(mtx, &det);
	res->col = res->row = mtx->col;
	clear_mtx(res);
	for (int i = 0; i < mtx->row; i++) {
		for (int j = 0; j < mtx->col; j++) {
			res->arr[i][j] = mtx_t.arr[i][j] / det;
		}
	}
	return true;
}


bool create_I(int n, matrix* I)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				I->arr[i][j] = 1;
			}
			else
			{
				I->arr[i][j] = 0;
			}
		}
	}
	I->col = I->row = n;
	return true;
}


typedef struct Kalman_filter
{
	matrix X; //x, y, Vx, Vy
	matrix P; //матрица недоверия вектору состояния Х
} Kalman_filter;

void default_init_kf(Kalman_filter* kf) {
	kf->X.arr[0][0] = 0; // x{0}
	kf->X.arr[1][0] = 0; // y{0}
	kf->X.arr[2][0] = 0; // Vx{0}
	kf->X.arr[3][0] = 0; // Vy{0}
	kf->X.col = 1;
	kf->X.row = 4;

	kf->P.col = 4;
	kf->P.row = 4;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			kf->P.arr[i][j] = 0;
		}
	}
	kf->P.arr[0][0] = 1.e9;
	kf->P.arr[1][1] = 1.e9;
	kf->P.arr[2][2] = 1.e9;
	kf->P.arr[3][3] = 1.e9;
}

void init_kf(Kalman_filter* kf, double std_x, double std_y, double std_vx, double std_vy, double x0, double y0, double vx0, double vy0)
{
	kf->X.arr[0][0] = x0; // x{0}
	kf->X.arr[1][0] = y0; // y{0}
	kf->X.arr[2][0] = vx0; // Vx{0}
	kf->X.arr[3][0] = vy0; // Vy{0}
	kf->X.col = 1;
	kf->X.row = 4;

	kf->P.col = 4;
	kf->P.row = 4;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			kf->P.arr[i][j] = 0;
		}
	}
	kf->P.arr[0][0] = std_x * std_x;
	kf->P.arr[1][1] = std_y * std_y;
	kf->P.arr[2][2] = std_vx * std_vx;
	kf->P.arr[3][3] = std_vy * std_vy;
}


void predict(Kalman_filter* kf, double dt, double std_a) 
{
	matrix F; //матрица перехода
	F.col = 4;
	F.row = 4;
	for(int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				F.arr[i][j] = 1;
			}
			else {
				F.arr[i][j] = 0;
			}
		}
	}
	F.arr[0][2] = dt;
	F.arr[1][3] = dt;
	matrix X_predict;
	product_mtx(&F, &kf->X, &X_predict);
	kf->X = X_predict;
	matrix Q;
	matrix F_transp;
	transp_mtx(&F, &F_transp);
	product_mtx(&F, &F_transp, &Q);
	for (int i = 0; i < Q.row; i++) {
		for (int j = 0; j < Q.col; j++) {
			Q.arr[i][j] *= (std_a * std_a);
		}
	}
	matrix P_predict;
	matrix FP;
	product_mtx(&F, &kf->P, &FP);
	matrix FPFt;
	product_mtx(&FP, &F_transp, &FPFt);
	sum_mtx(&FPFt, &Q, &P_predict);
	kf->P = P_predict;
	return;
}

void update(Kalman_filter* kf, matrix* Z, matrix* R)
{
	matrix H;
	H.row = 2;
	H.col = 4;
	for (int i = 0; i < H.row; i++) {
		for (int j = 0; j < H.col; j++) {
			if (i == j) {
				H.arr[i][j] = 1;
			}
			else
			{
				H.arr[i][j] = 0;
			}
		}
	}
	matrix Hx;
	product_mtx(&H, &kf->X, &Hx);
	matrix Y;
	diff_mtx(Z, &Hx, &Y);

	matrix HP;
	product_mtx(&H, &kf->P, &HP);
	matrix Ht;
	transp_mtx(&H, &Ht);
	matrix HPHt;
	product_mtx(&HP, &Ht, &HPHt);
	matrix S;
	sum_mtx(&HPHt, R, &S);

	matrix S_inv;
	inv_mtx(&S, &S_inv);
	matrix PH_t;
	product_mtx(&kf->P, &Ht, &PH_t);
	matrix K;
	product_mtx(&PH_t, &S_inv, &K);

	matrix KY;
	product_mtx(&K, &Y, &KY);
	matrix X_res;
	sum_mtx(&kf->X, &KY, &X_res);
	kf->X = X_res;

	matrix KH;
	product_mtx(&K, &H, &KH);
	matrix I;
	create_I(KH.row, &I);
	matrix I_KH;
	diff_mtx(&I, &KH, &I_KH);
	matrix P_res;
	product_mtx(&I_KH, &kf->P, &P_res);
	kf->P = P_res;

}


int main()
{
	double x0 = 53.2 ;
	double y0 = 18;
	double vx = 230;
	double vy = 30;
	double dt = 12;
	matrix Z_array[100];
	matrix R_array[100];
	for (int i = 0; i < 10; i++) {
		double t = dt * i;
		double noise_x = generateGaussianNoise(300);
		double noise_y = generateGaussianNoise(300);

		Z_array[i].row = 2;
		Z_array[i].col = 1;
		Z_array[i].arr[0][0] = x0 + (vx * t) + noise_x;
		Z_array[i].arr[1][0] = y0 + (vy * t) + noise_y;
		R_array[i].row = R_array[i].col = 2;
		clear_mtx(&R_array[i]);
		R_array[i].arr[0][0] = 1;
		R_array[i].arr[1][1] = 1;
	}
	Kalman_filter kf;
	init_kf(&kf, 1000, 1000, 1000, 1000, 0, 0, 0, 0);
	for (int i = 0; i < 10; i++) {
		predict(&kf, dt, 1);
		update(&kf, &Z_array[i], &R_array[i]);
		printf("%d\n", i);
		print_mtx(&kf.X);
	}
	print_mtx(&kf.X);
	matrix test;
	init_mtx(&test, 5, 6);
	print_mtx(&test);
	return 0;
}