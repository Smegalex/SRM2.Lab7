from sympy.abc import x, y
import math
from sympy import log, simplify, sqrt
from sympy.plotting import plot_implicit
import copy

sym_list = [x, y]


def spec_subst(func, subst: dict) -> float:
    if isinstance(func, int) or isinstance(func, float):
        return func

    result = simplify(func.subs(subst)).evalf()

    return result


def error_calc(q: float, xN: list, x: list) -> float:
    return (q/(1-q))*max(abs(xN[0]-x[0]), abs(xN[1]-x[1]))


def cycle(f1, f2, df1_div_dx1, df1_div_dx2, df2_div_dx1, df2_div_dx2, sym_list: list, start_approx: list, e: float, maxIter: int = 1000) -> list:
    # Область G
    G = {abs(sym_list[0]-start_approx[0]): start_approx[0],
         abs(sym_list[1]-start_approx[1]): start_approx[1]}

    dfs = [[df1_div_dx1, df1_div_dx2], [df2_div_dx1, df2_div_dx2]]

    Gx = 2*start_approx[0]
    Gy = 2*start_approx[1]

    # Перевірка збіжності
    max_row_sum = float("-inf")
    for row in dfs:
        substifiable = dict(zip(sym_list, [Gx, Gy]))

        d1 = abs(spec_subst(row[0], substifiable))
        d2 = abs(spec_subst(row[1], substifiable))
        row_sum = d1 + d2

        if row_sum > max_row_sum:
            max_row_sum = row_sum
    q = None
    if max_row_sum < 1:
        print("Розв'язання рівняння методом простої ітерації збігається якщо всі подальші наближення будуть лежати в області G.")
        q = max_row_sum
    else:
        raise ValueError(
            "Рівняння з наданим початковим наближенням не є збіжним за використання методу простих ітерацій.")

    k = 0
    error = float("inf")
    x = copy.deepcopy(start_approx)
    while (True):
        # print(f"Ітерація: {k}\n")
        # print(f"x1 = {x[0]}")
        # print(f"x2 = {x[1]}\n\n")

        if (error <= e):
            # print(f"Ітерація завершена на кроці {k}.")
            # print(f"x1 = {x[0]}")
            # print(f"x2 = {x[1]}\n\n")
            return x, k
        substifiable = dict(zip(sym_list, x))
        xN = [spec_subst(f1, substifiable), spec_subst(f2, substifiable)]

        # print(f"f1 = {xN[0]}")
        # print(f"f2 = {xN[1]}\n\n")

        substifiable = dict(zip(sym_list, xN))
        if list(G.keys())[0].subs(substifiable) > list(G.values())[0] or list(G.keys())[1].subs(substifiable) > list(G.values())[1]:
            raise ValueError(f"Наближення {k} не лежить в межах області G.")

        error = error_calc(q, xN, x)

        # print(f"e = {error}")
        # print("".ljust(40, "_"))

        k += 1
        x = copy.deepcopy(xN)
        if (k >= maxIter):
            raise TimeoutError(
                "Кількість ітерацій перевищує максимальну дозволену.")


if __name__ == "__main__":
    F1 = x**2-2*log(y, 10)-1
    F2 = x**2-x*y+1

    f1 = sqrt(2*log(y, 10)+1)
    f2 = x+1/x

    plot1 = plot_implicit(F1, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="red")
    plot2 = plot_implicit(F2, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="green")
    plot1.append(plot2[0])
    # plot1.show()

    df1_div_dx1 = -((2*log(y, 10)/(x**2))+(1/(x**2)))
    df1_div_dx2 = (2/x)*(1/(log(10)*y))
    df2_div_dx1 = 1-(1/(x**2))
    df2_div_dx2 = 0
    e = 0.0001

    start_approx = [1.25, 2.05]

    cycle(f1, f2, df1_div_dx1, df1_div_dx2, df2_div_dx1,
          df2_div_dx2, sym_list, start_approx, e)
