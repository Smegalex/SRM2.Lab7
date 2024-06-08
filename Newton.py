from sympy.abc import x, y
import math
from sympy import Matrix, log, simplify
from sympy.plotting import plot_implicit
import copy


def evalf_callable(el): return el.evalf()


sym_list = [x, y]


def subst_matrix(matrix: list, val_list: list) -> list:
    return Matrix(copy.deepcopy(matrix)
                  ).xreplace(dict(zip(sym_list, val_list))).applyfunc(simplify).applyfunc(evalf_callable).tolist()


def det_2x2(matrix: list) -> float:

    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def cycle(f1, f2, df1_div_dx1, df1_div_dx2, df2_div_dx1, df2_div_dx2, sym_list: list, start_approx: list, e: float, maxIter: int = 1000):

    xs = start_approx
    k = 0
    xsN = []

    A1 = [[f1, df1_div_dx2], [f2, df2_div_dx2]]
    A2 = [[df1_div_dx1, f1], [df2_div_dx1, f2]]
    J = [[df1_div_dx1, df1_div_dx2], [df2_div_dx1, df2_div_dx2]]

    detA1 = det_2x2(subst_matrix(A1, xs))
    detA2 = det_2x2(subst_matrix(A2, xs))
    detJ = det_2x2(subst_matrix(J, xs))

    # oldSubs = dict(zip(sym_list, xs))
    # print("".ljust(40, "_"))
    # print(f"iteration: 0\nx1 = {xs[0]}\nx2 = {xs[1]}\n\n")
    # print(
    #     f"f1 = {f1.subs(oldSubs).evalf(n=5)}\nf2 = {f2.subs(oldSubs).evalf(n=5)}\n\n")
    # print(
    #     f"df1/dx1 = {df1_div_dx1.subs(oldSubs).evalf(n=5).evalf(n=5)}\ndf2/dx1 = {df2_div_dx1.subs(oldSubs).evalf(n=5)}\n\n")
    # print(
    #     f"df1/dx2 = {df1_div_dx2.subs(oldSubs).evalf(n=5)}\ndf2/dx2 = {df2_div_dx2.subs(oldSubs).evalf(n=5)}\n\n")
    # print(f"detA1 = {detA1}\ndetA2 = {detA2}\n")
    # print(f"detJ = {detJ}\n")
    while True:
        xsN.append(xs[0]-(detA1/detJ))
        xsN.append(xs[1]-(detA2/detJ))

        maxError = float("-inf")
        for i in range(len(xsN)):
            error = abs(xsN[i]-xs[i])
            if error >= maxError:
                maxError = error

        if maxError <= e:
            # print("    End:")
            # print(xsN)
            return xsN, k

        k += 1
        detA1 = det_2x2(subst_matrix(A1, xsN))
        detA2 = det_2x2(subst_matrix(A2, xsN))
        detJ = det_2x2(subst_matrix(J, xsN))

        # newSubs = dict(zip(sym_list, xsN))
        # print("".ljust(40, "_"))
        # print(f"iteration: {k}\nx1 = {xsN[0]}\nx2 = {xsN[1]}\n\n")
        # print(
        #     f"f1 = {f1.subs(newSubs).evalf(n=5)}\nf2 = {f2.subs(newSubs).evalf(n=5)}\n\n")
        # print(
        #     f"df1/dx1 = {df1_div_dx1.subs(newSubs).evalf(n=5)}\ndf2/dx1 = {df2_div_dx1.subs(newSubs).evalf(n=5)}\n\n")
        # print(
        #     f"df1/dx2 = {df1_div_dx2.subs(newSubs).evalf(n=5)}\ndf2/dx2 = {df2_div_dx2.subs(newSubs).evalf(n=5)}\n\n")
        # print(f"detA1 = {detA1}\ndetA2 = {detA2}\n\n")
        # print(f"detJ = {detJ}\n\n\n")

        xs = xsN
        xsN = []

        if (k >= maxIter):
            raise TimeoutError(
                "Кількість ітерацій перевищує максимальну дозволену.")


if __name__ == "__main__":
    f1 = x**2-2*log(y, 10)-1
    f2 = x**2-x*y+1

    plot1 = plot_implicit(f1, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="red")
    plot2 = plot_implicit(f2, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="green")
    plot1.append(plot2[0])
    # plot1.show()

    df1_div_dx1 = 2*x
    df1_div_dx2 = -2/(log(10)*y)
    df2_div_dx1 = 2*x
    df2_div_dx2 = -x
    e = 0.0001

    start_approx = [1.25, 2.05]

    cycle(f1, f2, df1_div_dx1, df1_div_dx2, df2_div_dx1,
          df2_div_dx2, sym_list, start_approx, e)
