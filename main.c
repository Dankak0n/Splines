#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>



#define MAX_POINTS 100

typedef struct {
    double x, y;
} point;

typedef struct {
    double A, B, C, D;
    point begin, end;
} segment_of_spline;

typedef struct {
    int size;
    segment_of_spline * list_y, * list_x;
} spline;

double get_dist_pp(point a, point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

const double INIT_APP = 0.5; //Начальные значения иксов
const double EPS = 0.00001;
const int CNT = 100; //Количество итераций

typedef double ** num_type;

void allocate_memory(void *** matrix, int row, int col, int size_type, int size_ptr) {
    (*matrix) = malloc(size_ptr * (row));
    for (int i = 0; i < row; i++) {
        (*matrix)[i] = malloc(col * size_type);
    }
}

num_type mat_mat_mul(num_type matrix_a, int row_a, int col_a,
                     num_type matrix_b, int row_b, int col_b) {
    num_type ans;
    allocate_memory(&ans, row_a, col_b, sizeof(double), sizeof(double *));
    for (int i = 0; i < row_a; ++i) {
        for (int j = 0; j < col_b; ++j) {
            ans[i][j] = 0;
            for (int k = 0; k < col_a; ++k) {
                ans[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
    return ans;
}

num_type mat_num_mul(num_type matrix, int row, int col, double number) {
    num_type ans;
    allocate_memory(&ans, row, col, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[i][j] = matrix[i][j] * number;
        }
    }
    return ans;
}

num_type mat_mat_add(num_type matrix_a, num_type matrix_b, int row, int col) {
    num_type ans;
    allocate_memory(&ans, row, col, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[i][j] = matrix_a[i][j] + matrix_b[i][j];
        }
    }
    return ans;
}

num_type get_transposed_matrix(num_type matrix, int row, int col) {
    num_type ans;
    allocate_memory(&ans, col, row, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[j][i] = matrix[i][j];
        }
    }
    return ans;
}

num_type get_cofactor(num_type matrix, int order, int row, int col) {
    num_type ans;
    allocate_memory(&ans, order - 1, order - 1, sizeof(double), sizeof(double *));
    for (int i = 0; i < order - 1; ++i) {
        for (int j = 0; j < order - 1; ++j) {
            int x, y;
            x = i + (i >= row);
            y = j + (j >= col);
            ans[i][j] = matrix[x][y];
        }
    }
    return ans;
}

double get_determinant(num_type matrix, int order) {
    if (order == 1) {
        return matrix[0][0];
    }
    double ans = 0;
    for (int i = 0; i < order; ++i) {
        int sign = ((i & 1) ? -1 : 1);
        ans += (double)sign * matrix[0][i] * get_determinant(get_cofactor(matrix, order, 0, i), order - 1);
    }
    return ans;
}

num_type get_union_matrix(num_type matrix, int order) {
    num_type ans;
    allocate_memory(&ans, order, order, sizeof(double), sizeof(double *));
    if (order == 1) {
        ans[0][0] = matrix[0][0];
        return ans;
    }
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            int sign = ((i + j & 1) ? -1 : 1);
            ans[i][j] = (double)sign * get_determinant(get_cofactor(matrix, order, i, j), order - 1);
        }
    }
    ans = get_transposed_matrix(ans, order, order);
    return ans;
}

num_type get_inversed_matrix(num_type matrix, int order) {
    num_type ans;
    allocate_memory(&ans, order, order, sizeof(double), sizeof(double *));
    ans = mat_num_mul(get_union_matrix(matrix, order), order, order, 1.0 / get_determinant(matrix, order));
    return ans;
}

num_type get_func_value() {

}

num_type get_mat_yakobi(double t1, double t2, segment_of_spline * sx1, segment_of_spline * sy1,
                                              segment_of_spline * sx2, segment_of_spline * sy2) {
    num_type ans;
    allocate_memory(&ans, 2, 2, sizeof(double), sizeof(double *));
    ans[0][1] = ans[1][0] = 0.0;
    ans[0][0] = pow(sx1->A * pow(t1, 2) + sx1->B * t1 + sx1->C, 2) +
                pow(sy1->A * pow(t1, 2) + sy1->B * t1 + sy1->C, 2);
    ans[1][1] = pow(sx2->A * pow(t2, 2) + sx2->B * t2 + sx2->C, 2) +
                pow(sy2->A * pow(t2, 2) + sy2->B * t2 + sy2->C, 2);
    ans[0][0] *= 2.0;
    ans[1][1] *= 2.0;
    return ans;
}


void output_matrix(num_type matrix, int row, int col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (j == 0) {
                if (row == 1)
                    printf("( ");
                else if (i == 0)
                    printf("/ ");
                else if (i == row - 1)
                    printf("\\ ");
                else
                    printf("| ");
            }
            printf("%.6lf", matrix[i][j]);
            if (i == row - 1 && j == col - 1)
                printf("  ");
            else
                printf(", ");
            if (j == col - 1) {
                if (row == 1)
                    printf(")");
                else if (i == 0)
                    printf("\\");
                else if (i == row - 1)
                    printf("/");
                else
                    printf("|");
            }
        }
        puts("");
    }
}





void init(int * count_of_segments, point p_arr[], FILE * f) {
    fprintf(stderr, "init...\n");
    if (f != NULL) {
        fscanf(f, "%d", count_of_segments);
    } else {
        scanf("%d", count_of_segments);
    }

    for (int i = 0; i <= *count_of_segments; ++i) {
        if (f != NULL) {
            fscanf(f, "%lf %lf", &(p_arr[i].x), &(p_arr[i].y));
        } else {
            scanf("%lf %lf", &(p_arr[i].x), &(p_arr[i].y));
        }
    }
}

void set_values(int count_of_segments,
                point p_arr[],
                double h[],
                double parameter_t[],
                double diagonal_matrix[MAX_POINTS][MAX_POINTS]) {
    parameter_t[0] = h[0] = 0.0;
    for (int i = 0; i < count_of_segments; ++i) {
        h[i + 1] = get_dist_pp(p_arr[i], p_arr[i + 1]);
    }
    for (int i = 1; i <= count_of_segments; ++i) {
        parameter_t[i] = parameter_t[i - 1] + h[i];
    }
    for (int i = 1; i < count_of_segments; ++i) {
        diagonal_matrix[i - 1][i] =
        diagonal_matrix[i][i - 1] = h[i] / 6.0;
        diagonal_matrix[i][i] = (h[i] + h[i + 1]) / 3.0;
    }
}

segment_of_spline * tridiagonal_matrix_algorithm(int count_of_segments,
                                                 double h[],
                                                 double diagonal_matrix[MAX_POINTS][MAX_POINTS],
                                                 double parameter_t[],
                                                 double current_points[]) {
    double gamma[MAX_POINTS] = {0},
           B[MAX_POINTS] = {0},
           U[MAX_POINTS] = {0},
           V[MAX_POINTS] = {0};
    for (int i = 1; i < count_of_segments; ++i) {
        B[i] = (current_points[i + 1] - current_points[i]) / h[i + 1] -
               (current_points[i] - current_points[i - 1]) / h[i];
    }

    V[1] = diagonal_matrix[1][2] / (-diagonal_matrix[1][1]);
    U[1] = (-B[1]) / (-diagonal_matrix[1][1]);

    for (int i = 2; i < count_of_segments; ++i) {
        V[i] =  diagonal_matrix[i][i + 1] /
               (-diagonal_matrix[i][i] - diagonal_matrix[i][i - 1] * V[i - 1]);
        U[i] = (diagonal_matrix[i][i - 1] * U[i - 1] - B[i]) /
               (-diagonal_matrix[i][i] - diagonal_matrix[i][i - 1] * V[i - 1]);
    }

    gamma[count_of_segments - 1] = U[count_of_segments - 1];

    for (int i = count_of_segments - 2; i > 0; --i) {
        gamma[i] = V[i] * gamma[i + 1] + U[i];
    }

    segment_of_spline * ans;
    ans = (segment_of_spline *)malloc(sizeof(segment_of_spline) * count_of_segments);
    for (int i = 0; i < count_of_segments; ++i) {
        ans[i].A = (gamma[i + 1] - gamma[i]) / (6.0 * h[i + 1]);
        ans[i].B = (gamma[i] * parameter_t[i + 1] - gamma[i + 1] * parameter_t[i]) / (2.0 * h[i + 1]);
        ans[i].C = (6.0 * current_points[i + 1] -
                    6.0 * current_points[i] +
                    3.0 * gamma[i + 1] * pow(parameter_t[i], 2) -
                    3.0 * gamma[i] * pow(parameter_t[i + 1], 2) +
                    gamma[i] * pow(h[i + 1], 2) -
                    gamma[i + 1] * pow(h[i + 1], 2))
        / (6.0 * h[i + 1]);
        ans[i].D = (6.0 * current_points[i] * parameter_t[i + 1] -
                    6.0 * current_points[i + 1] * parameter_t[i] +
                    gamma[i] * pow(parameter_t[i + 1], 3) -
                    gamma[i + 1] * pow(parameter_t[i], 3) +
                    pow(h[i + 1], 2) * parameter_t[i] * gamma[i + 1] -
                    pow(h[i + 1], 2) * parameter_t[i + 1] * gamma[i])
        / (6.0 * h[i + 1]);
    }
    return ans;
}

spline get_spline(int count_of_segments, point p_arr[]) {
    fprintf(stderr, "getting spline...\n");
    double h[MAX_POINTS],
           parameter_t[MAX_POINTS],
           diagonal_matrix[MAX_POINTS][MAX_POINTS];

    set_values(count_of_segments, p_arr, h, parameter_t, diagonal_matrix);

    double x_arr[MAX_POINTS], y_arr[MAX_POINTS];
    for (int i = 0; i <= count_of_segments; ++i) {
        x_arr[i] = p_arr[i].x;
        y_arr[i] = p_arr[i].y;
    }
    spline ans;
    ans.size = count_of_segments;
    ans.list_x = tridiagonal_matrix_algorithm(count_of_segments, h, diagonal_matrix, parameter_t, x_arr);
    ans.list_y = tridiagonal_matrix_algorithm(count_of_segments, h, diagonal_matrix, parameter_t, y_arr);

    for (int i = 0; i < count_of_segments; ++i) {
        ans.list_x[i].begin = p_arr[i];
        ans.list_x[i].end = p_arr[i + 1];

        ans.list_y[i].begin = p_arr[i];
        ans.list_y[i].end = p_arr[i + 1];
    }

    return ans;
}



void iteration(double t[]) {
    num_type m_yakobi = get_mat_yakobi(n, vector_x, vector_func);
    ans = mat_mat_add(vector_x, mat_num_mul(mat_mat_mul(get_inversed_matrix(m_yakobi, n), n, n, get_func_value(vector_x, vector_func, n), n, 1), n, 1, -1), n, 1);
}

double get_dist_seg_seg(segment_of_spline) {
    double t[2] = {0};

    for (int i = 0; i < CNT; ++i) {
        iteration(t);
    }
}

double get_dist_ss(spline * sp1, spline * sp2) {

    double t1[MAX_POINTS] = {0},
           t2[MAX_POINTS] = {0};
    for (int i = 1; i <= sp1->size; ++i) {
        t1[i] = t1[i - 1] +
        get_dist_pp((point){sp1->list_x[i - 1].begin.x, sp1->list_x[i - 1].begin.y},
                    (point){sp1->list_x[i - 1].end.x, sp1->list_x[i - 1].end.y});
    }

    for (int i = 1; i <= sp2->size; ++i) {
        t2[i] = t2[i - 1] +
        get_dist_pp((point){sp2->list_x[i - 1].begin.x, sp2->list_x[i - 1].begin.y},
                    (point){sp2->list_x[i - 1].end.x, sp2->list_x[i - 1].end.y});
    }


    double ans = 100000;
    for (int i = 0; i < sp1->size; ++i) {
        for (int j = 0; j < sp2->size; ++j) {
            double cur_dist = get_dist_seg_seg(&(sp1->list_x[i]), &(sp1->list_y[i]),
                                               &(sp2->list_x[j]), &(sp2->list_y[j]));
            ans = (cur_dist < ans ? cur_dist : ans);
        }
    }

    return ans;
}

void output_spline(spline * sp, FILE * f) {
    fprintf(stderr, "transfering info to python script...\n");

    double parameter_t[MAX_POINTS] = {0};
    for (int i = 1; i <= sp->size; ++i) {
        parameter_t[i] = parameter_t[i - 1] +
        get_dist_pp((point){sp->list_x[i - 1].begin.x, sp->list_x[i - 1].begin.y},
                    (point){sp->list_x[i - 1].end.x, sp->list_x[i - 1].end.y});
    }

    if (f != NULL)
        fprintf(f, "%d\n", sp->size);
    else
        printf("%d\n", sp->size);
    for (int i = 0; i <= sp->size; ++i) {
        if (f != NULL)
            fprintf(f, "%.15lf ", parameter_t[i]);
        else
            printf("%.15lf ", parameter_t[i]);
    }
    if (f != NULL)
        fprintf(f, "\n");
    else
        printf("\n");
    for (int i = 0; i < sp->size; ++i) {
        char mask[] = "(%.15lf)*t^3+(%.15lf)*t^2+(%.15lf)*t+(%.15lf)\n";
        if (f != NULL) {
            fprintf(f, mask, sp->list_x[i].A, sp->list_x[i].B, sp->list_x[i].C, sp->list_x[i].D);
            fprintf(f, mask, sp->list_y[i].A, sp->list_y[i].B, sp->list_y[i].C, sp->list_y[i].D);
        } else {
            printf(mask, sp->list_x[i].A, sp->list_x[i].B, sp->list_x[i].C, sp->list_x[i].D);
            printf(mask, sp->list_y[i].A, sp->list_y[i].B, sp->list_y[i].C, sp->list_y[i].D);
        }
    }
}

int main() {
    FILE * fr = fopen("input.txt", "r");
    FILE * fw = fopen("D://Dankakon//pythonProject/to_draw.txt", "w");

    int count_of_segments1, count_of_segments2;
    point p_arr1[MAX_POINTS], p_arr2[MAX_POINTS];

    init(&count_of_segments1, p_arr1, fr);
    spline sp1 = get_spline(count_of_segments1, p_arr1);
    output_spline(&sp1, fw);

    init(&count_of_segments2, p_arr2, fr);
    spline sp2 = get_spline(count_of_segments2, p_arr2);
    output_spline(&sp2, fw);

    double dist = get_dist_ss(&sp1, &sp2);
    fprintf(stderr, ">>>dist = %.10lf\n", dist);

    fprintf(stderr, ">>>is_intersect = %s\n", (dist < EPS ? "true" : "false"));

    fclose(fw);
    fclose(fr);

    fprintf(stderr, "calling python scrypt...");
    system("python D://Dankakon//pythonProject/main.py");
    return 0;
}
