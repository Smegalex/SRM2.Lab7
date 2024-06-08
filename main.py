from sympy.abc import x, y
from sympy import log, sqrt, pprint
from sympy.plotting import plot_implicit

from Newton import cycle as nwcl
from Iterations import cycle as itcl


sym_list = [x, y]


if __name__ == "__main__":
    F1 = x**2-2*log(y, 10)-1
    F2 = x**2-x*y+1

    print(f"Маємо систему рівнянь з двома змінними:")
    pprint(F1)
    print("= 0\n")
    pprint(F2)
    print("= 0\n")


    plot1 = plot_implicit(F1, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="red")
    plot2 = plot_implicit(F2, (x, 1, 1.5), (y, 1.8, 2.2),
                          show=False, line_color="green")
    plot1.append(plot2[0])
    # print("Графічний розв'язок рівнянь виглядає як:")
    # plot1.show()
    
    e = 0.0001

    df1_div_dx1 = 2*x
    df1_div_dx2 = -2/(log(10)*y)
    df2_div_dx1 = 2*x
    df2_div_dx2 = -x

    
    start_approx = [1.25, 2.05]
    print(f"Початкове наближення обрали як:")
    print(f"x1 = {start_approx[0]}\nx2 = {start_approx[1]}")
    print("".ljust(40, "_"),"\n")


    newtSol, newtK = nwcl(F1, F2, df1_div_dx1, df1_div_dx2, df2_div_dx1,
          df2_div_dx2, sym_list, start_approx, e)
    
    print(f"Методом Н'ютона рівняння розв'язали за {newtK} ітерацій та отримали:")
    print(f"x1 = {newtSol[0]}\nx2 = {newtSol[1]}")
    print("".ljust(40, "_"),"\n")
    

    f1 = sqrt(2*log(y, 10)+1)
    f2 = x+1/x

    print("Систему рівнянь привели до вигляду:")
    print("x =")
    pprint(f1)
    print("")
    print("y =")
    pprint(f2)

    df1_div_dx1 = -((2*log(y, 10)/(x**2))+(1/(x**2)))
    df1_div_dx2 = (2/x)*(1/(log(10)*y))
    df2_div_dx1 = 1-(1/(x**2))
    df2_div_dx2 = 0


    iterSol, iterK = itcl(f1, f2, df1_div_dx1, df1_div_dx2, df2_div_dx1,
          df2_div_dx2, sym_list, start_approx, e)
    
    print(f"Методом простих ітерацій рівняння розв'язали за {iterK} ітерацій та отримали:")
    print(f"x1 = {iterSol[0]}\nx2 = {iterSol[1]}\n\n")
