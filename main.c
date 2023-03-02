#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

#define INFINITY 10000
#define MAX_POINTS 100
#define EPS 0.0001
#define CNT 100

typedef struct {
    double x, y;
} point;

typedef struct {
    point first, second;
} two_points;

typedef struct {
    double A, B, C, D;
    point begin, end;
} segment_of_spline;

typedef struct {
    int size;
    segment_of_spline * list_y, * list_x;
} spline;


double get_dist_p_p(point a, point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
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
        h[i + 1] = get_dist_p_p(p_arr[i], p_arr[i + 1]);
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


double get_der(segment_of_spline * seg, double number, int pr) {
    if (pr == 0) {
        return (double)((seg->A) * pow(number, 3) + (seg->B) * pow(number, 2) + (seg->C) * number + (seg->D));
    }
    if (pr == 1) {
        return (double)(3.0 * (seg->A) * pow(number, 2) + 2.0 * (seg->B) * number + (seg->C));
    }
    if (pr == 2) {
        return (double)(6.0 * (seg->A) * number + 2.0 * (seg->B));
    }
}

two_points get_dist_seg_seg(segment_of_spline * segx1, segment_of_spline * segy1, double t1_first, double t1_second,
                            segment_of_spline * segx2, segment_of_spline * segy2, double t2_first, double t2_second) {
    double t[2] = {(t1_first + t1_second) / 2, (t2_first + t2_second) / 2};
    for (int i = 0; i < CNT; ++i) {
        double t1 = t[0], t2 = t[1];

        t[0] = t1 - (((get_der(segx2, t2, 0) - get_der(segx1, t1, 0)) *
                       get_der(segx2, t2, 2) + pow(get_der(segx2, t2, 1), 2) +
                      (get_der(segy2, t2, 0) - get_der(segy1, t1, 0)) *
                       get_der(segy2, t2, 2) + pow(get_der(segy2, t2, 1), 2)) *
                      ((get_der(segx1, t1, 0) - get_der(segx2, t2, 0)) *
                       get_der(segx1, t1, 1) +
                      (get_der(segy1, t1, 0) - get_der(segy2, t2, 0)) *
                       get_der(segy1, t1, 1)) +
                      ((get_der(segx1, t1, 1) * get_der(segx2, t2, 1) +
                        get_der(segy1, t1, 1) * get_der(segy2, t2, 1)) *
                      ((get_der(segx2, t2, 0) - get_der(segx1, t1, 0)) *
                        get_der(segx2, t2, 1) +
                       (get_der(segy2, t2, 0) - get_der(segy1, t1, 0)) *
                        get_der(segy2, t2, 1)))) /
                        (
                        ((get_der(segx1, t1, 0) - get_der(segx2, t2, 0)) *
                          get_der(segx1, t1, 2) + pow(get_der(segx1, t1, 1), 2) +

                         (get_der(segy1, t1, 0) - get_der(segy2, t2, 0)) *
                          get_der(segy1, t1, 2) + pow(get_der(segy1, t1, 1), 2)) *

                        ((get_der(segx2, t2, 0) - get_der(segx1, t1, 0)) *
                          get_der(segx2, t2, 2) + pow(get_der(segx2, t2, 1), 2) +

                         (get_der(segy2, t2, 0) - get_der(segy1, t1, 0)) *
                          get_der(segy2, t2, 2) + pow(get_der(segy2, t2, 1), 2)) -
                          pow(get_der(segx1, t1, 1) * get_der(segx2, t2, 1) +
                              get_der(segy1, t1, 1) * get_der(segy2, t2, 1), 2)
                         );

        t[1] = t2 -   (((get_der(segx1, t1, 0) - get_der(segx2, t2, 0)) *
                       get_der(segx1, t1, 2) + pow(get_der(segx1, t1, 1), 2) +
                      (get_der(segy1, t1, 0) - get_der(segy2, t2, 0)) *
                       get_der(segy1, t1, 2) + pow(get_der(segy1, t1, 1), 2)) *
                      ((get_der(segx2, t2, 0) - get_der(segx1, t1, 0)) *
                       get_der(segx2, t2, 1) +
                      (get_der(segy2, t2, 0) - get_der(segy1, t1, 0)) *
                       get_der(segy2, t2, 1)) +
                      ((get_der(segx2, t2, 1) * get_der(segx1, t1, 1) +
                        get_der(segy2, t2, 1) * get_der(segy1, t1, 1)) *
                      ((get_der(segx1, t1, 0) - get_der(segx2, t2, 0)) *
                        get_der(segx1, t1, 1) +
                       (get_der(segy1, t1, 0) - get_der(segy2, t2, 0)) *
                        get_der(segy1, t1, 1)))) /
                        (
                        ((get_der(segx2, t2, 0) - get_der(segx1, t1, 0)) *
                          get_der(segx2, t2, 2) + pow(get_der(segx2, t2, 1), 2) +

                         (get_der(segy2, t2, 0) - get_der(segy1, t1, 0)) *
                          get_der(segy2, t2, 2) + pow(get_der(segy2, t2, 1), 2)) *

                        ((get_der(segx1, t1, 0) - get_der(segx2, t2, 0)) *
                          get_der(segx1, t1, 2) + pow(get_der(segx1, t1, 1), 2) +

                         (get_der(segy1, t1, 0) - get_der(segy2, t2, 0)) *
                          get_der(segy1, t1, 2) + pow(get_der(segy1, t1, 1), 2)) -
                          pow(get_der(segx2, t2, 1) * get_der(segx1, t1, 1) +
                              get_der(segy2, t2, 1) * get_der(segy1, t1, 1), 2)
                         );

        if (t[0] < t1_first) {
            t[0] = t1_first;
        }
        if (t[0] > t1_second) {
            t[0] = t1_second;
        }
        if (t[1] < t2_first) {
            t[1] = t2_first;
        }
        if (t[1] > t2_second) {
            t[1] = t2_second;
        }
    }

    two_points ans;
    ans.first.x = (segx1->A) * pow(t[0], 3) + (segx1->B) * pow(t[0], 2) + (segx1->C) * t[0] + (segx1->D);
    ans.first.y = (segy1->A) * pow(t[0], 3) + (segy1->B) * pow(t[0], 2) + (segy1->C) * t[0] + (segy1->D);

    ans.second.x = (segx2->A) * pow(t[1], 3) + (segx2->B) * pow(t[1], 2) + (segx2->C) * t[1] + (segx2->D);
    ans.second.y = (segy2->A) * pow(t[1], 3) + (segy2->B) * pow(t[1], 2) + (segy2->C) * t[1] + (segy2->D);

    return ans;
}

two_points get_dist_s_s(spline * sp1, spline * sp2) {
    double t1[MAX_POINTS] = {0},
           t2[MAX_POINTS] = {0};
    for (int i = 1; i <= sp1->size; ++i) {
        t1[i] = t1[i - 1] +
        get_dist_p_p((point){sp1->list_x[i - 1].begin.x, sp1->list_x[i - 1].begin.y},
                    (point){sp1->list_x[i - 1].end.x, sp1->list_x[i - 1].end.y});
    }

    for (int i = 1; i <= sp2->size; ++i) {
        t2[i] = t2[i - 1] +
        get_dist_p_p((point){sp2->list_x[i - 1].begin.x, sp2->list_x[i - 1].begin.y},
                    (point){sp2->list_x[i - 1].end.x, sp2->list_x[i - 1].end.y});
    }

    double min_dist = INFINITY;
    two_points ans;
    for (int i = 0; i < sp1->size; ++i) {
        for (int j = 0; j < sp2->size; ++j) {
            two_points now = get_dist_seg_seg(&(sp1->list_x[i]), &(sp1->list_y[i]), t1[i], t1[i + 1],
                                              &(sp2->list_x[j]), &(sp2->list_y[j]), t2[j], t2[j + 1]);
            double cur_dist = get_dist_p_p(now.first, now.second);
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                ans = now;
            }
        }
    }

    return ans;
}

void output_spline(spline * sp, FILE * f) {
    fprintf(stderr, "transfering info to python script...\n");

    double parameter_t[MAX_POINTS] = {0};
    for (int i = 1; i <= sp->size; ++i) {
        parameter_t[i] = parameter_t[i - 1] +
        get_dist_p_p((point){sp->list_x[i - 1].begin.x, sp->list_x[i - 1].begin.y},
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

    two_points pp = get_dist_s_s(&sp1, &sp2);
    double dist = get_dist_p_p(pp.first, pp.second);
    fprintf(stderr, ">>>dist = %.10lf\n>>>points: (%lf ,  %lf)  (%lf ,  %lf)\n", dist, pp.first.x, pp.first.y,  pp.second.x, pp.second.y);

    if (dist < EPS) {
        fprintf(stderr, ">>>is_intersect = true:(%lf,  %lf)\n", pp.first.x, pp.first.y);
    } else {
        fprintf(stderr, ">>>is_intersect = false\n");
    }

    fclose(fw);
    fclose(fr);

    fprintf(stderr, "calling python scrypt...");
    system("python D://Dankakon//pythonProject/main.py");
    return 0;
}
