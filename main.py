import math

import sympy as sp
x = sp.symbols('x')


Newton_Raphson_iterations_limit = 100
secant_iterations_limit = 100

def Bisection_Method(polinom,start_point,end_point,epsilon):
    """
    Finds root of function by bisection method and returns number of iterations and the root
    :param polinom: polinom to be solved
    :param start_point: start point
    :param end_point: end point
    :param epsilon: epsilon when to stop
    :return: number of iterations and the result
    """
    f = sp.lambdify(x, polinom)
    if not f(start_point)*f(end_point)<0:
        return None,None
    c = (start_point+end_point)/2
    i = 0
    max_iteration_num = -(math.log(epsilon/(end_point-start_point))/math.log(2))
    while abs(end_point-start_point) >= epsilon:
        print(str(i+1) + "   a: " + str(start_point) + " b: " + str(end_point) + " c: " + str(c) + " f(a): " + str(f(start_point)) + " f(b): " + str(f(end_point)) + " f(c): " + str(f(c)))
        if f(c)*f(start_point)<0:
            end_point = c
            c = (start_point+c)/2
        elif f(c)*f(end_point)<0:
            start_point = c;
            c = (end_point + c) / 2
        else:
            if f(start_point) == 0:
                return i,start_point
            elif f(end_point) == 0:
                return i,end_point
            else:
                return i,c
        i += 1
        if(i > max_iteration_num+1):
            print("This function is not suitable for biselection method")
            return None,None
    print(str(i + 1) + "   a: " + str(start_point) + " b: " + str(end_point) + " c: " + str(c) + " f(a): " + str(f(start_point)) + " f(b): " + str(f(end_point)) + " f(c): " + str(f(c)))
    print("Number of iterations: " + str(i+1) + "  The solution: " + str(c))
    return i+1,c

def Newton_Raphson(polinomToIterate,start_point,end_point,epsilon):
    """
    Finds root of function by Newton Raphson method and returns number of iterations and the root
    :param polinom: polinom to be solved
    :param start_point: start point
    :param end_point: end point
    :param epsilon: epsilon when to stop
    :return: number of iterations and the result
    """
    polinomDiff = sp.diff(polinomToIterate, x)
    iterX = (start_point+end_point)/2
    f = sp.lambdify(x, polinomToIterate)
    fDiff = sp.lambdify(x, polinomDiff)
    i = 0
    while abs(f(iterX)) >= epsilon:
        print(str(i+1) + "   Xr = " + str(iterX) + " f(Xr) = " + str(f(iterX)) + " f'(Xr) = " + str(fDiff(iterX)))
        if fDiff(iterX) == 0:
            print("This function is not suitable for Newton Raphson method")
            return None,None
        iterX = iterX - (f(iterX)/fDiff(iterX))
        i += 1
        if(i>Newton_Raphson_iterations_limit):
            print("This function is not suitable for Newton Raphson method")
            return None, None
    print(str(i + 1) + "   Xr = " + str(iterX) + " f(Xr) = " + str(f(iterX)) + " f'(Xr) = " + str(fDiff(iterX)))
    print("Number of iterations: " + str(i + 1) + "  The solution: " + str(iterX))
    return i+1,iterX
def secant_method(polinom,start_point,end_point,epsilon):
    """
    Finds root of function by secant method and returns number of iterations and the root
    :param polinom: polinom to be solved
    :param start_point: start point
    :param end_point: end point
    :param epsilon: epsilon when to stop
    :return: number of iterations and the result
    """
    iterX = start_point
    iterXNext = end_point
    f = sp.lambdify(x, polinom)
    i = 0
    while abs(f(iterX)) >= epsilon:
        print(str(i + 1) + "   Xr = " + str(iterX) + " Xr+1 = " + str(f(iterXNext)) + " f(Xr) = " + str(f(iterX)))
        if f(iterX) - f(iterXNext) == 0:
            print("This function is not suitable for secant method")
            return None, None
        temp = iterXNext
        iterXNext = (iterX*f(iterXNext)-iterXNext*f(iterX))/(f(iterXNext)-f(iterX))
        iterX = temp
        i += 1
        if (i > secant_iterations_limit):
            print("This function is not suitable for secant method")
            return None, None
    print(str(i + 1) + "   Xr = " + str(iterX) + " Xr+1 = " + str(f(iterXNext)) + " f(Xr) = " + str(f(iterX)))
    print("Number of iterations: " + str(i + 1) + "  The solution: " + str(iterX))
    return i + 1, iterX

choice = -1
while choice != '0':
    choice = input("Enter 1 for bisection method, 2 for Neuton Rafson, 3 for secant method or 0 to exit: ")
    if choice == '1':
        iteration_arr = []
        solutions_arr = []
        f = sp.lambdify(x, polinom)
        prev_f = f(start_point)
        can_be_solved = True
        doingDiff = False
        polinomDiff = sp.diff(polinom,x)
        for k in range(2):
            i = start_point
            prev_i = i
            i += size_of_small_section
            while i <= end_point:
                curr_f = f(i)
                if prev_f*curr_f<0:
                    if not doingDiff:
                        num_of_iterations,solution = Bisection_Method(polinom,prev_i,i,epsilon)
                    else:
                        num_of_iterations, solution = Bisection_Method(polinomDiff, prev_i, i, epsilon)
                    if num_of_iterations == None:
                        can_be_solved = False
                        break
                    if not doingDiff:
                        iteration_arr.append(num_of_iterations)
                        solutions_arr.append(solution)
                    else:
                        if abs(solution) < epsilon:
                            iteration_arr.append(num_of_iterations)
                            solutions_arr.append(solution)
                        else:
                            print("Solution is not 0, so it's not valid solution")
                prev_i = i
                i += size_of_small_section
                prev_f = f(prev_i)
            if not can_be_solved:
                break
            f = sp.diff(polinom,x)
            f = sp.lambdify(x, f)
            prev_f = f(start_point)
            if not doingDiff:
                print("---Working on differential of function---")
                doingDiff = True
        if can_be_solved:
            i = 0
            while i < len(iteration_arr):
                print("Solution " + str(i + 1) + " - Number of iterations: " + str(
                    iteration_arr[i]) + ", Solution: " + str(solutions_arr[i]))
                i += 1
    elif choice == '2':
        iteration_arr = []
        solutions_arr = []
        f = sp.lambdify(x, polinom)
        prev_f = f(start_point)
        can_be_solved = True
        doingDiff = False
        polinomDiff = sp.diff(polinom, x)
        for k in range(2):
            i = start_point
            prev_i = i
            i += size_of_small_section
            while i <= end_point:
                curr_f = f(i)
                if prev_f * curr_f < 0:
                    if not doingDiff:
                        num_of_iterations, solution = Newton_Raphson(polinom, prev_i, i, epsilon)
                    else:
                        num_of_iterations, solution = Newton_Raphson(polinomDiff, prev_i, i, epsilon)
                    if num_of_iterations == None:
                        can_be_solved = False
                        break
                    if not doingDiff:
                        iteration_arr.append(num_of_iterations)
                        solutions_arr.append(solution)
                    else:
                        if abs(solution) < epsilon:
                            iteration_arr.append(num_of_iterations)
                            solutions_arr.append(solution)
                        else:
                            print("Solution is not 0, so it's not valid solution")
                prev_i = i
                i += size_of_small_section
                prev_f = f(prev_i)
            if not can_be_solved:
                break
            f = sp.diff(polinom, x)
            f = sp.lambdify(x, f)
            prev_f = f(start_point)
            if not doingDiff:
                print("---Working on differential of function---")
                doingDiff = True
        if can_be_solved:
            i = 0
            while i < len(iteration_arr):
                print("Solution " + str(i + 1) + " - Number of iterations: " + str(
                    iteration_arr[i]) + ", Solution: " + str(solutions_arr[i]))
                i += 1
    elif choice == '3':
        iteration_arr = []
        solutions_arr = []
        f = sp.lambdify(x, polinom)
        prev_f = f(start_point)
        can_be_solved = True
        doingDiff = False
        polinomDiff = sp.diff(polinom, x)
        for k in range(2):
            i = start_point
            prev_i = i
            i += size_of_small_section
            while i <= end_point:
                curr_f = f(i)
                if prev_f * curr_f < 0:
                    if not doingDiff:
                        num_of_iterations, solution = secant_method(polinom, prev_i, i, epsilon)
                    else:
                        num_of_iterations, solution = secant_method(polinomDiff, prev_i, i, epsilon)
                    if num_of_iterations == None:
                        can_be_solved = False
                        break
                    if not doingDiff:
                        iteration_arr.append(num_of_iterations)
                        solutions_arr.append(solution)
                    else:
                        if abs(solution) < epsilon:
                            iteration_arr.append(num_of_iterations)
                            solutions_arr.append(solution)
                        else:
                            print("Solution is not 0, so it's not valid solution")
                prev_i = i
                i += size_of_small_section
                prev_f = f(prev_i)
            if not can_be_solved:
                break
            f = sp.diff(polinom, x)
            f = sp.lambdify(x, f)
            prev_f = f(start_point)
            if not doingDiff:
                print("---Working on differential of function---")
                doingDiff = True
        if can_be_solved:
            i = 0
            while i < len(iteration_arr):
                print("Solution " + str(i + 1) + " - Number of iterations: " + str(
                    iteration_arr[i]) + ", Solution: " + str(solutions_arr[i]))
                i += 1
    elif choice != '0':
        print("Invalid input")
