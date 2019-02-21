''' Preemptive formulation for the interval-based resource
constrained scheduling problem.
Written by David Torres Sanchez, 2019: d.torressanchez@lancaster.ac.uk'''

import time
import gurobipy as grb


def get_constants(it):
    def _floor(x, y):
        try:
            return floor(x / y)
        except ZeroDivisionError:
            return 0

    from lcalg import process
    from math import floor
    start = time.time()
    n, T, K, p, R, r, E, V, A, LB, ES, LS = process(it)
    duration = time.time() - start
    V = range(1, n)
    p_minus = [min(dur, floor(dur / 2)) for dur in p]
    N = int(sum([_floor(p[i], p_minus[i]) for i in range(len(p))]))
    # Updated number of events for preemptions
    E = range(N - 2)
    return n, N, T, K, p, p_minus, R, r, E, V, A, LB, ES, LS, duration


def generate_constraints(it):

    def _create_variables():
        z = model.addVars([(i, e) for i in V for e in range(-1, N - 2)],
                          name="z", vtype="B")
        a = model.addVars([(i, e) for i in V for e in E], name="a",
                          vtype="C", lb=0.0)
        t = model.addVars([e for e in E], name="t", vtype="C",
                          lb=0.0, ub=T)
        C_max = model.addVar(name="C_max", vtype="C",
                                                 lb=int(LB), ub=T, obj=1)
        model.update()
        return z, a, t, C_max

    def _add_oee_constraints():
        model.addConstrs((grb.quicksum(z[i, e] for e in E) >= 1 for i in V),
                         name="(42)")
        model.addConstr(t[E[-1]] == C_max, name="(43)")

        model.addConstr(t[0] == 0, name="(44)")

        model.addConstrs(
            (t[e + 1] - t[e] >= 0 for e in E if e != N - 3), name="(45)")

        model.addConstrs((p[i] >= t[f] - t[e] >=
                          p_minus[i] * ((z[i, e] - z[i, e - 1]) -
                                        (z[i, f] - z[i, f - 1]) - 1)
                          for i in V
                          for e in E
                          for f in range(e + 1, len(E)) if e > 0 and f > 0),
                         name="(46)")
        model.addConstrs((grb.quicksum(z[i, e_1] for e_1 in range(e)) -
                          (e * (1 - (z[i, e] - z[i, e - 1]))) <= 0
                          for i in V for e in E if e != 0), name="(47)")
        model.addConstrs((grb.quicksum(z[i, e_1] for e_1 in range(e, N - 2)) -
                          ((N - 2 - e) * (1 + (z[i, e] - z[i, e - 1]))) <= 0
                          for i in V for e in E if e != 0), name="(48)")
        model.addConstrs((z[i, e] +
                          grb.quicksum(z[j, e_1] for e_1 in range(e + 1)) -
                          (e * (1 - z[i, e])) <= 1
                          for (i, j) in A for e in E), name="(49)")
        model.addConstrs((grb.quicksum(r[i][k] * z[i, e] for i in V) <= R[k]
                          for k in range(K)
                          for e in E), name="(50)")
        model.addConstrs((z[i, -1] == 0 for i in V), name="(51)")
        model.addConstrs((ES['%s' % i] * z[i, e] <= t[e]
                          for i in V for e in E), name="(52.1)")
        model.addConstrs((LS['%s' % i] * (z[i, e] - z[i, e - 1]) +
                          LS['%s' % (n - 1)] * (1 - (z[i, e] - z[i, e - 1])) >=
                          t[e] for i in V for e in E), name="(52.2)")
        model.addConstr(ES['%s' % n], "<=", rhs=C_max, name="(53.1)")
        model.addConstr(LS['%s' % n], ">=", rhs=C_max, name="(53.2)")
        return

    def _add_preemption_constraints():
        # Preemption constraints
        BIGM = 1e5
        model.addConstrs((a[i, 0] == 0 for i in V), name="(1.2.5a)")  # init
        model.addConstrs((a[i, e - 1] <= a[i, e] for i in V
                          for e in E if e != 0),
                         name="(1.2.5b)")
        model.addConstrs((a[i, e] <= a[i, e - 1] +
                          BIGM * (z[i, e] + z[i, e - 1])
                          for i in V for e in E if e != 0),
                         name="(1.2.5c)")
        model.addConstrs((a[i, e] >= a[i, e - 1] + (t[e] - t[e - 1]) -
                          BIGM * (1 - z[i, e])
                          for i in V for e in E if e != 0),
                         name="(1.2.5d)")
        model.addConstrs((a[i, e] >= a[i, e - 1] + (t[e] - t[e - 1]) -
                          BIGM * (1 - z[i, e - 1])
                          for i in V for e in E if e != 0),
                         name="(1.2.5d)")
        model.addConstrs((a[i, E[-1]] == p[i]
                          for i in V), name="(1.2.5f)")  # last
        model.addConstrs((a[i, e] <= p[i] for i in V for e in E),
                         name="(1.2.5g)")  # upper bound
        model.update()
        return

    model = grb.Model("Linear Program %s" % it)
    model.setParam('TimeLimit', 5 * 60)
    # load constants
    n, N, T, K, p, p_minus, R, r, E, V, A, LB, ES, LS, duration_prec = \
        get_constants(it)
    # Create variables
    z, a, t, C_max = _create_variables()
    # Constraints
    _add_oee_constraints()
    _add_preemption_constraints()
    return model, E[-1], duration_prec


def optimise(model, it, duration_prec):
    start = time.time()
    model.optimize()
    runtime = time.time() - start
    if model.getAttr("Status") == 3:  # Infeasible
        model.computeIIS()
        model.write("./output/model_%s.ilp" % it)
        return
    else:
        model.write("./output/solution_%s.sol" % it)
        file_out = open("./output/j301_results.txt", "a")
        file_out.write('%s \t %s \t %s \t %s \t %s \n ' %
                       (it, model.objVal, runtime,
                        model.Runtime, duration_prec))
        file_out.close()
        return
        # t_last = model.getVarByName("t[%s]" % last)
        # print t_last.x
        # raw_input()


def main():
    for it in range(1, 11):
        model, last, duration_prec = generate_constraints(it)
        optimise(model, it, duration_prec)


if __name__ == main():
    main()
