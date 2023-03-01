import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
f = open('D://Dankakon//pythonProject/to_draw.txt', 'r')
spline_cnt = f.readline()
while spline_cnt != "":
    t_arr = list(map(float, (f.readline()).split()))
    for i in range(int(spline_cnt)):
        t = sp.Symbol('t')
        function_x = sp.sympify(str(f.readline()))
        function_y = sp.sympify(str(f.readline()))
        interval = np.arange(t_arr[i], t_arr[i + 1], 0.1)
        x_values = [function_x.subs(t, value) for value in interval]
        y_values = [function_y.subs(t, value) for value in interval]
        plt.plot(x_values, y_values)
    spline_cnt = f.readline()
plt.show()
f.close()
