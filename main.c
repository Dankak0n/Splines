#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define MAX_POINTS 100
#define EPS 0.0001
#define CNT 100

typedef struct {
    double x, y;
} point;

typedef struct {
    double begin, end;
} pair_of_double;

typedef struct {
    point begin, end;
} segment;

typedef struct {
    double A, B, C, D;
} cooficients;

typedef struct {
    cooficients x, y;
    segment borders;
    pair_of_double t_param;
} piece_of_spline;

typedef struct {
    int size;
    piece_of_spline * list;
} spline;

const char input_file_path[] = "D://Dankakon//Programming//Cbiba//input.txt";
const char to_draw_file_path[] = "D://Dankakon//pythonProject//to_draw.txt";
const char python_script_path[] = "D://Dankakon//pythonProject//main.py";

double get_segment_length(segment arg) {
    fprintf(stderr, "   +[ getting segment length ]\n");
    double diff_x = arg.begin.x - arg.end.x,
           diff_y = arg.begin.y - arg.end.y;
    return sqrt(pow(diff_x, 2) + pow(diff_y, 2));
}

void init(int * count_of_segments, point coordinates[], FILE * f) {
    fprintf(stderr, "+[ initialization ]\n");

    if (f != NULL) {
        fscanf(f, "%d", count_of_segments);
    } else {
        scanf("%d", count_of_segments);
    }

    for (int i = 0; i <= *count_of_segments; ++i) {
        if (f != NULL) {
            fscanf(f, "%lf %lf", &(coordinates[i].x), &(coordinates[i].y));
        } else {
            scanf("%lf %lf", &(coordinates[i].x), &(coordinates[i].y));
        }
    }
}

cooficients * get_cooficients(int count_of_segments, spline * sp, double coordinates[]) {
    fprintf(stderr, "+[ getting cooficients ]\n");

    double h[MAX_POINTS] = {0}, diagonal_matrix[MAX_POINTS][MAX_POINTS] = {0};
    for (int i = 0; i < count_of_segments; ++i) {
        h[i + 1] = get_segment_length(sp->list[i].borders);
    }

    for (int i = 1; i < count_of_segments; ++i) {
        diagonal_matrix[i - 1][i] =
        diagonal_matrix[i][i - 1] = h[i] / 6.0;
        diagonal_matrix[i][i] = (h[i] + h[i + 1]) / 3.0;
    }

    double gamma[MAX_POINTS] = {0}, B[MAX_POINTS] = {0}, U[MAX_POINTS] = {0}, V[MAX_POINTS] = {0};
    for (int i = 1; i < count_of_segments; ++i) {
        B[i] = (coordinates[i + 1] - coordinates[i]) / h[i + 1] -
               (coordinates[i] - coordinates[i - 1]) / h[i];
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

    double parameter_t[MAX_POINTS] = {0};
    for (int i = 0; i < count_of_segments; ++i) {
        parameter_t[i + 1] = sp->list[i].t_param.end;
    }

    cooficients * ans;
    ans = (cooficients *)malloc(sizeof(cooficients) * count_of_segments);
    for (int i = 0; i < count_of_segments; ++i) {
        ans[i].A = (gamma[i + 1] - gamma[i]) / (6.0 * h[i + 1]);
        ans[i].B = (gamma[i] * parameter_t[i + 1] - gamma[i + 1] * parameter_t[i]) / (2.0 * h[i + 1]);
        ans[i].C = (6.0 * coordinates[i + 1] -
                    6.0 * coordinates[i] +
                    3.0 * gamma[i + 1] * pow(parameter_t[i], 2) -
                    3.0 * gamma[i] * pow(parameter_t[i + 1], 2) +
                    gamma[i] * pow(h[i + 1], 2) -
                    gamma[i + 1] * pow(h[i + 1], 2))
                    / (6.0 * h[i + 1]);
        ans[i].D = (6.0 * coordinates[i] * parameter_t[i + 1] -
                    6.0 * coordinates[i + 1] * parameter_t[i] +
                    gamma[i] * pow(parameter_t[i + 1], 3) -
                    gamma[i + 1] * pow(parameter_t[i], 3) +
                    pow(h[i + 1], 2) * parameter_t[i] * gamma[i + 1] -
                    pow(h[i + 1], 2) * parameter_t[i + 1] * gamma[i])
                    / (6.0 * h[i + 1]);
    }
    return ans;
}

spline get_spline(int count_of_segments, point coordinates[]) {
    fprintf(stderr, "+[ getting spline ]\n");

    spline ans;
    ans.size = count_of_segments;
    ans.list = (piece_of_spline *)malloc(sizeof(piece_of_spline) * count_of_segments);
    for (int i = 0; i < count_of_segments; ++i) {
        ans.list[i].borders.begin = coordinates[i];
        ans.list[i].borders.end = coordinates[i + 1];
        ans.list[i].t_param.begin = (i == 0 ? 0.0 : ans.list[i - 1].t_param.end);
        ans.list[i].t_param.end = ans.list[i].t_param.begin + get_segment_length(ans.list[i].borders);
    }

    double x_arr[MAX_POINTS], y_arr[MAX_POINTS];
    for (int i = 0; i <= count_of_segments; ++i) {
        x_arr[i] = coordinates[i].x;
        y_arr[i] = coordinates[i].y;
    }

    cooficients * x_cf = get_cooficients(count_of_segments, &ans, x_arr);
    cooficients * y_cf = get_cooficients(count_of_segments, &ans, y_arr);

    for (int i = 0; i < count_of_segments; ++i) {
        ans.list[i].x = x_cf[i];
        ans.list[i].y = y_cf[i];
    }
    return ans;
}

double get_der(cooficients cf, double number, int pr) {
//    fprintf(stderr, "   +[ getting derivative ]\n");

    if (pr == 2) {
        return (double)(6.0 * (cf.A) * number + 2.0 * (cf.B));
    }
    if (pr == 1) {
        return (double)(3.0 * (cf.A) * pow(number, 2) + 2.0 * (cf.B) * number + (cf.C));
    }
    return (double)((cf.A) * pow(number, 3) + (cf.B) * pow(number, 2) + (cf.C) * number + (cf.D));
}

double do_some_magic(piece_of_spline * seg1, piece_of_spline * seg2, double t1, double t2) {
    double magic_var1 = get_der(seg2->x, t2, 0) - get_der(seg1->x, t1, 0),
           magic_var2 = get_der(seg2->y, t2, 0) - get_der(seg1->y, t1, 0),
           magic_var3 = get_der(seg1->x, t1, 1) * get_der(seg2->x, t2, 1) +
                        get_der(seg1->y, t1, 1) * get_der(seg2->y, t2, 1);
    return t1 -  ((magic_var1 * get_der(seg2->x, t2, 2) +
                   magic_var2 * get_der(seg2->y, t2, 2) +
                   pow(get_der(seg2->x, t2, 1), 2) +
                   pow(get_der(seg2->y, t2, 1), 2)) *
                 ((-magic_var1) * get_der(seg1->x, t1, 1) +
                  (-magic_var2) * get_der(seg1->y, t1, 1)) +
                  (magic_var3 *
                  (magic_var1 * get_der(seg2->x, t2, 1) +
                   magic_var2 * get_der(seg2->y, t2, 1)))) /
                (((-magic_var1) * get_der(seg1->x, t1, 2) +
                  (-magic_var2) * get_der(seg1->y, t1, 2) +
                   pow(get_der(seg1->x, t1, 1), 2) +
                   pow(get_der(seg1->y, t1, 1), 2)) *
                  (magic_var1 * get_der(seg2->x, t2, 2) +
                   magic_var2 * get_der(seg2->y, t2, 2) +
                   pow(get_der(seg2->x, t2, 1), 2) +
                   pow(get_der(seg2->y, t2, 1), 2)) -
                   pow(magic_var3, 2));
}

segment get_pieces_dist(piece_of_spline * seg1, piece_of_spline * seg2) {
    fprintf(stderr, "+[ getting pieces distance ]\n");

    double t[2] = {(seg1->t_param.begin + seg1->t_param.end) / 2, (seg2->t_param.begin + seg2->t_param.end) / 2};
    for (int i = 0; i < CNT; ++i) {
        double t1 = t[0], t2 = t[1];

        t[0] = do_some_magic(seg1, seg2, t1, t2);
        t[1] = do_some_magic(seg2, seg1, t2, t1);

        if (t[0] < seg1->t_param.begin) {
            t[0] = seg1->t_param.begin;
        }
        if (t[0] > seg1->t_param.end) {
            t[0] = seg1->t_param.end;
        }
        if (t[1] < seg2->t_param.begin) {
            t[1] = seg2->t_param.begin;
        }
        if (t[1] > seg1->t_param.end) {
            t[1] = seg1->t_param.end;
        }
    }

    segment ans;
    ans.begin.x = (seg1->x.A) * pow(t[0], 3) + (seg1->x.B) * pow(t[0], 2) + (seg1->x.C) * t[0] + (seg1->x.D);
    ans.begin.y = (seg1->y.A) * pow(t[0], 3) + (seg1->y.B) * pow(t[0], 2) + (seg1->y.C) * t[0] + (seg1->y.D);

    ans.end.x = (seg2->x.A) * pow(t[1], 3) + (seg2->x.B) * pow(t[1], 2) + (seg2->x.C) * t[1] + (seg2->x.D);
    ans.end.y = (seg2->y.A) * pow(t[1], 3) + (seg2->y.B) * pow(t[1], 2) + (seg2->y.C) * t[1] + (seg2->y.D);

    return ans;
}

segment get_closest_points(spline * sp1, spline * sp2) {
    fprintf(stderr, "+[ getting minimum distance ]\n");

    double min_dist = INFINITY;
    segment ans;
    for (int i = 0; i < sp1->size; ++i) {
        for (int j = 0; j < sp2->size; ++j) {
            segment now = get_pieces_dist(&(sp1->list[i]), &(sp2->list[j]));
            double cur_dist = get_segment_length(now);
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                ans = now;
            }
        }
    }
    return ans;
}

void output_spline_info(spline * sp, FILE * f) {
    if (f != NULL) {
        fprintf(stderr, "+[ transfering info to python script ]\n");
    } else {
        fprintf(stderr, "+[ outputting ]\n");
    }

    if (f != NULL) {
        fprintf(f, "%d\n%.15lf ", sp->size, 0.0);
        for (int i = 0; i < sp->size; ++i) {
            fprintf(f, "%.15lf ", sp->list[i].t_param.end);
        }
        fprintf(f, "\n");
        for (int i = 0; i < sp->size; ++i) {
            char mask[] = "(%.15lf)*t^3+(%.15lf)*t^2+(%.15lf)*t+(%.15lf)\n";
            cooficients cf;
            cf = sp->list[i].x;
            fprintf(f, mask, cf.A, cf.B, cf.C, cf.D);
            cf = sp->list[i].y;
            fprintf(f, mask, cf.A, cf.B, cf.C, cf.D);
        }
    } else {
        for (int i = 0; i < sp->size; ++i) {
            printf("%lf < t%d < %lf\n", sp->list[i].t_param.begin, i + 1, sp->list[i].t_param.end);
        }
        printf("\n");
        for (int i = 0; i < sp->size; ++i) {
            printf("(");
            char mask[] = "(%lf)*t^3+(%lf)*t^2+(%lf)*t+(%lf)";
            cooficients cf;

            cf = sp->list[i].x;
            printf(mask, cf.A, cf.B, cf.C, cf.D);
            printf(", ");

            cf = sp->list[i].y;
            printf(mask, cf.A, cf.B, cf.C, cf.D);
            printf(")\n");
        }
    }
}

void draw_spline(FILE * fr, FILE * fw) {
    fclose(fr);
    if (fw != NULL) {
        fclose(fw);
        fprintf(stderr, "+[ calling python script ]\n");
        char sys_commad[] = "python ";
        strcat(sys_commad, python_script_path);
        system(sys_commad);
    } else {
        fprintf(stderr, "-[ not able to draw ]\n");
    }
}

void output_splines_distance(spline * sp1, spline * sp2) {
    fprintf(stderr, "+[ outputting splines distance ]\n");

    segment closest_points = get_closest_points(sp1, sp2);
    double dist = get_segment_length(closest_points);
    printf("> distance = %.15lf\n", dist);
    printf("> is intersect = %s\n", (dist < EPS ? "true" : "false"));
    if (dist < EPS) {
        printf("> point of intersection = (%.15lf , %.15lf)\n", closest_points.begin.x, closest_points.begin.y);
    } else {
        printf("> closest points = (%.15lf , %.15lf), (%.15lf , %.15lf)\n", closest_points.begin.x, closest_points.begin.y,
                                                                            closest_points.end.x, closest_points.end.y);
    }
}

int main() {
    FILE * fr = fopen(input_file_path, "r");
    FILE * fw = fopen(to_draw_file_path, "w");

    int count_of_segments1, count_of_segments2;
    point p_arr1[MAX_POINTS], p_arr2[MAX_POINTS];

    init(&count_of_segments1, p_arr1, fr);
    init(&count_of_segments2, p_arr2, fr);

    spline sp1 = get_spline(count_of_segments1, p_arr1);
    spline sp2 = get_spline(count_of_segments2, p_arr2);

    output_spline_info(&sp1, fw);
    output_spline_info(&sp2, fw);

    output_splines_distance(&sp1, &sp2);

    draw_spline(fr, fw);
    fprintf(stderr, "+[ end ]\n");
    return 0;
}
