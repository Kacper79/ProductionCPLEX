#include <ProductionCPLEX.h>

ILOSTLBEGIN

struct Instance {
    int K = 0; // number of products
    int T = 0; // number of periods

    std::vector<double> hc;   // holding cost serviceables
    std::vector<double> hcr;  // holding cost recoverables
    std::vector<double> pc;   // production cost
    std::vector<double> pcr;  // remanufacturing cost
    std::vector<double> sc;   // setup cost production
    std::vector<double> scr;  // setup cost remanufacturing

    std::vector<double> tp;   // time production per unit
    std::vector<double> tpr;  // time remanufacturing per unit
    std::vector<double> ts;   // time setup production
    std::vector<double> tsr;  // time setup remanufacturing

    std::vector<double> ct;   // capacity production
    std::vector<double> ctr;  // capacity remanufacturing

    std::vector<std::vector<double>> d; // demand per product per period
    std::vector<std::vector<double>> r; // returns per product per period

    std::vector<double> y0;   // initial serviceable inventory
    std::vector<double> yr0;  // initial recoverable inventory
    bool useValidIneq = true; // enable valid inequalities
};

struct Schedule {
    int product = -1;

    // original schedule cost used in the master objective
    double cost = 0.0;

    // capacity consumptions in each period for master constraints
    std::vector<double> prodCap;
    std::vector<double> remCap;

    // detailed plan (kept so this can later be inspected / reused)
    std::vector<double> Q, Qr, Y, Yr, gamma, gammar;

    bool isDummy = false;
};

struct PricingResult {
    bool improving = false;
    double reducedCost = 0.0;
    Schedule column;
};

struct MasterResult {
    bool solved = false;
    double objective = 0.0;
    std::vector<double> sigma;  // duals of choose-one constraints
    std::vector<double> piProd; // duals of production capacities
    std::vector<double> piRem;  // duals of reman capacities
    std::vector<std::vector<double>> deltaValues;
};

double max2(double a, double b) { return (a > b) ? a : b; }
double min2(double a, double b) { return (a < b) ? a : b; }

Instance buildReferenceInstance() {

    Instance ins;

    ins.K = 4;
    ins.T = 5;

    ins.hc = { 1, 1, 1, 1 };
    ins.hcr = { 0.5, 0.5, 0.5, 0.5 };

    ins.pc = { 0, 0, 0, 0 };
    ins.pcr = { 0, 0, 0, 0 };

    ins.sc = { 500, 500, 500, 500 };
    ins.scr = { 500, 500, 500, 500 };

    ins.tp = { 1, 1, 1, 1 };
    ins.tpr = { 1, 1, 1, 1 };
    ins.ts = { 20, 20, 20, 20 };
    ins.tsr = { 20, 20, 20, 20 };

    ins.ct = { 300, 300, 300, 300, 300 };
    ins.ctr = { 300, 300, 300, 300, 300 };

    ins.y0 = { 0, 0, 0, 0 };
    ins.yr0 = { 0, 0, 0, 0 };

    ins.useValidIneq = true;

    ins.d = {
        { 40,  70,  60,  70,  60},
        {100,  60, 180, 170, 110},
        { 80,  90,  90, 140, 140},
        { 30,  80,  50,  80,  90}
    };

    ins.r = {
        {20, 30, 30, 20, 30},
        {50, 40, 60, 40, 40},
        {40, 60, 50, 30, 40},
        {30, 30, 30, 50, 20}
    };

    return ins;
}

// Maximum possible production quantity given available capacity
std::vector<std::vector<double>> computeBigM(const Instance& ins) {
    std::vector<std::vector<double>> M(ins.K, std::vector<double>(ins.T, 0.0));

    for (int k = 0; k < ins.K; ++k) {
        for (int t = 0; t < ins.T; ++t) {
            const double numerator = max2(0.0, max2(ins.ct[t] - ins.ts[k], ins.ctr[t] - ins.tsr[k]));
            const double denominator = min2(ins.tp[k], ins.tpr[k]);

            if (denominator <= 0.0) {
                throw std::runtime_error("Parameters tp, tpr must be > 0");
            }
            M[k][t] = numerator / denominator;
        }
    }
    return M;
}

// Dummy schedule with high penalty cost (ensures feasibility)
Schedule makeDummySchedule(const Instance& ins, int k, double penalty = 1.0e5) {
    Schedule s;
    s.product = k;
    s.cost = penalty;
    s.prodCap.assign(ins.T, 0.0);
    s.remCap.assign(ins.T, 0.0);
    s.Q.assign(ins.T, 0.0);
    s.Qr.assign(ins.T, 0.0);
    s.Y.assign(ins.T, 0.0);
    s.Yr.assign(ins.T, 0.0);
    s.gamma.assign(ins.T, 0.0);
    s.gammar.assign(ins.T, 0.0);
    s.isDummy = true;
    return s;
}

// Solves SLULSP-RM for one product to find improving column
PricingResult solvePricingSubproblem( const Instance& ins, const std::vector<std::vector<double>>& M, int k, double sigma, const std::vector<double>& piProd, const std::vector<double>& piRem, double rcTolerance = -1e-5) {
    IloEnv env;
    PricingResult result;

    try {
        IloModel model(env);

        std::vector<IloNumVar> Q(ins.T), Qr(ins.T), Y(ins.T), Yr(ins.T);
        std::vector<IloBoolVar> gamma(ins.T), gammar(ins.T);

        for (int t = 0; t < ins.T; ++t) {
            Q[t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
            Qr[t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
            Y[t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
            Yr[t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
            gamma[t] = IloBoolVar(env);
            gammar[t] = IloBoolVar(env);
        }

        // Reduced-cost objective from the paper:
        // original cost - sigma_k - sum_t pi_t * prodCap_t - sum_t pi^r_t * remCap_t
        IloExpr originalCost(env);
        IloExpr reducedCostExpr(env);

        for (int t = 0; t < ins.T; ++t) {
            originalCost += ins.hc[k] * Y[t] + ins.hcr[k] * Yr[t]
                + ins.pc[k] * Q[t] + ins.pcr[k] * Qr[t]
                + ins.sc[k] * gamma[t] + ins.scr[k] * gammar[t];

            reducedCostExpr += ins.hc[k] * Y[t] + ins.hcr[k] * Yr[t]
                + ins.pc[k] * Q[t] + ins.pcr[k] * Qr[t]
                + ins.sc[k] * gamma[t] + ins.scr[k] * gammar[t]
                - piProd[t] * (ins.tp[k] * Q[t] + ins.ts[k] * gamma[t])
                - piRem[t] * (ins.tpr[k] * Qr[t] + ins.tsr[k] * gammar[t]);
        }
        reducedCostExpr -= sigma;

        model.add(IloMinimize(env, reducedCostExpr));

        // inventory balances
        model.add(ins.y0[k] + Q[0] + Qr[0] == ins.d[k][0] + Y[0]);
        for (int t = 1; t < ins.T; ++t) {
            model.add(Y[t - 1] + Q[t] + Qr[t] == ins.d[k][t] + Y[t]);
        }

        model.add(ins.yr0[k] + ins.r[k][0] == Qr[0] + Yr[0]);
        for (int t = 1; t < ins.T; ++t) {
            model.add(Yr[t - 1] + ins.r[k][t] == Qr[t] + Yr[t]);
        }

        // linking constraints
        for (int t = 0; t < ins.T; ++t) {
            model.add(Q[t] <= M[k][t] * gamma[t]);
            model.add(Qr[t] <= M[k][t] * gammar[t]);
        }

        // optional VI
        if (ins.useValidIneq) {
            // (18) for one product k
            for (int t = 0; t <= ins.T - 2; ++t) {
                for (int p = 1; p <= ins.T - 1 - t; ++p) {
                    IloExpr rhs(env);
                    for (int s = t + 1; s <= t + p; ++s) {
                        rhs += ins.d[k][s];
                        rhs -= M[k][s] * gamma[s];
                        rhs -= M[k][s] * gammar[s];
                    }
                    model.add(Y[t] >= rhs);
                    rhs.end();
                }
            }

            // (19) for one product k
            for (int t = 1; t < ins.T; ++t) {
                for (int p = 1; p <= t; ++p) {
                    IloExpr rhs(env);
                    for (int s = t - p; s <= t; ++s) {
                        rhs += ins.r[k][s];
                        rhs -= M[k][s] * gammar[s];
                    }
                    model.add(Yr[t] >= rhs);
                    rhs.end();
                }
            }
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);

        if (!cplex.solve()) {
            throw std::runtime_error("Pricing subproblem was not solved successfully.");
        }

        result.reducedCost = cplex.getObjValue();
        result.improving = (result.reducedCost < rcTolerance); //a schedule is improving if reduced cost is negative

        if (result.improving) {
            Schedule col;
            col.product = k;
            col.cost = cplex.getValue(originalCost); //the original (positive) cost is important while the reduced (negative) cost is being minimized (both are proportional to each other)
            col.prodCap.assign(ins.T, 0.0);
            col.remCap.assign(ins.T, 0.0);
            col.Q.assign(ins.T, 0.0);
            col.Qr.assign(ins.T, 0.0);
            col.Y.assign(ins.T, 0.0);
            col.Yr.assign(ins.T, 0.0);
            col.gamma.assign(ins.T, 0.0);
            col.gammar.assign(ins.T, 0.0);
            col.isDummy = false;

            for (int t = 0; t < ins.T; ++t) {
                col.Q[t] = cplex.getValue(Q[t]);
                col.Qr[t] = cplex.getValue(Qr[t]);
                col.Y[t] = cplex.getValue(Y[t]);
                col.Yr[t] = cplex.getValue(Yr[t]);
                col.gamma[t] = cplex.getValue(gamma[t]);
                col.gammar[t] = cplex.getValue(gammar[t]);
                col.prodCap[t] = ins.tp[k] * col.Q[t] + ins.ts[k] * col.gamma[t];
                col.remCap[t] = ins.tpr[k] * col.Qr[t] + ins.tsr[k] * col.gammar[t];
            }
            result.column = std::move(col);
        }

        originalCost.end();
        reducedCostExpr.end();
        env.end();
        return result;
    }
    catch (...) {
        env.end();
        throw;
    }
}

// Solves LP relaxation of master problem (Set Partitioning Problem)
MasterResult solveRelaxedMasterProblemSPP(const Instance& datasetInstance,const std::vector<std::vector<Schedule>>& schedulesPool) {
    IloEnv env;
    MasterResult res;

    try {
        IloModel model(env);

        std::vector<std::vector<IloNumVar>> delta(datasetInstance.K);
        IloExpr objective(env);

        for (int k = 0; k < datasetInstance.K; ++k) {
            delta[k].reserve(static_cast<std::size_t>(schedulesPool[k].size()));
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                std::ostringstream name;
                name << "delta_k" << (k + 1) << "_n" << (n + 1);

                delta[k].push_back(IloNumVar(env, 0.0, 1.0, ILOFLOAT, name.str().c_str()));
                objective += schedulesPool[k][n].cost * delta[k].back(); // min Z = sum_k sum_n_from_S nc_kn * delta_kn
            }
        }
        model.add(IloMinimize(env, objective));

        IloRangeArray chooseOne(env);
        IloRangeArray prodCapConstr(env);
        IloRangeArray remCapConstr(env);

        // (23) exactly one schedule per product
        for (int k = 0; k < datasetInstance.K; ++k) {
            IloExpr lhs(env);
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                lhs += delta[k][n];
            }

            std::ostringstream cname;
            cname << "chooseOne_k" << (k + 1);

            chooseOne.add(IloRange(env, 1.0, lhs, 1.0, cname.str().c_str()));
            lhs.end();
        }
        model.add(chooseOne);

        // (24) capacity constraints for production
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < datasetInstance.K; ++k) {
                for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                    lhs += schedulesPool[k][n].prodCap[t] * delta[k][n];
                }
            }

            std::ostringstream cname;
            cname << "prodCap_t" << (t + 1);

            prodCapConstr.add(IloRange(env, -IloInfinity, lhs, datasetInstance.ct[t], cname.str().c_str()));
            lhs.end();
        }
        model.add(prodCapConstr);

        // (25) capacity constraints for remanufacturing
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < datasetInstance.K; ++k) {
                for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                    lhs += schedulesPool[k][n].remCap[t] * delta[k][n];
                }
            }
            std::ostringstream cname;
            cname << "remCap_t" << (t + 1);

            remCapConstr.add(IloRange(env, -IloInfinity, lhs, datasetInstance.ctr[t], cname.str().c_str()));
            lhs.end();
        }
        model.add(remCapConstr);

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);

        if (!cplex.solve()) {
            throw std::runtime_error("Restricted master LP was not solved successfully.");
        }

        // Extract dual variables for pricing
        res.solved = true;
        res.objective = cplex.getObjValue();
        res.sigma.resize(datasetInstance.K, 0.0);
        res.piProd.resize(datasetInstance.T, 0.0);
        res.piRem.resize(datasetInstance.T, 0.0);
        res.deltaValues.resize(datasetInstance.K);

        for (int k = 0; k < datasetInstance.K; ++k) {
            res.sigma[k] = cplex.getDual(chooseOne[k]);
            res.deltaValues[k].resize(schedulesPool[k].size(), 0.0);
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                res.deltaValues[k][n] = cplex.getValue(delta[k][n]);
            }
        }
        for (int t = 0; t < datasetInstance.T; ++t) {
            res.piProd[t] = cplex.getDual(prodCapConstr[t]);
            res.piRem[t] = cplex.getDual(remCapConstr[t]);
        }

        objective.end();
        env.end();
        return res;
    }
    catch (...) {
        env.end();
        throw;
    }
}

// Solves integer master problem over generated columns (upper bound)
std::vector<std::vector<double>> solveIntegerMasterProblemSPP( const Instance& datasetInstance, const std::vector<std::vector<Schedule>>& schedulesPool, double& objectiveValue) {
    IloEnv env;

    try {
        IloModel model(env);
        std::vector<std::vector<IloIntVar>> delta(datasetInstance.K);
        IloExpr objective(env);

        for (int k = 0; k < datasetInstance.K; ++k) {
            delta[k].reserve(static_cast<std::size_t>(schedulesPool[k].size()));
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                std::ostringstream name;
                name << "delta_k" << (k + 1) << "_n" << (n + 1);

                delta[k].push_back(IloNumVar(env, 0.0, 1.0, ILOINT, name.str().c_str()));
                objective += schedulesPool[k][n].cost * delta[k].back(); // min Z = sum_k sum_n_from_S nc_kn * delta_kn
            }
        }
        model.add(IloMinimize(env, objective));

        // Exactly one schedule per product
        for (int k = 0; k < datasetInstance.K; ++k) {
            IloExpr lhs(env);
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) lhs += delta[k][n];
            model.add(lhs == 1.0);
            lhs.end();
        }

        // Capacity constraints
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhsProd(env), lhsRem(env);
            for (int k = 0; k < datasetInstance.K; ++k) {
                for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                    lhsProd += schedulesPool[k][n].prodCap[t] * delta[k][n];
                    lhsRem += schedulesPool[k][n].remCap[t] * delta[k][n];
                }
            }
            model.add(lhsProd <= datasetInstance.ct[t]);
            model.add(lhsRem <= datasetInstance.ctr[t]);
            lhsProd.end();
            lhsRem.end();
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);

        if (!cplex.solve()) {
            throw std::runtime_error("Integer master over generated columns could not be solved.");
        }

        objectiveValue = cplex.getObjValue();
        std::vector<std::vector<double>> vals(datasetInstance.K);

        for (int k = 0; k < datasetInstance.K; ++k) {
            vals[k].resize(schedulesPool[k].size(), 0.0);

            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                vals[k][n] = cplex.getValue(delta[k][n]);
            }
        }

        objective.end();
        env.end();
        return vals;
    }
    catch (...) {
        env.end();
        throw;
    }
}

// Printing schedule details
void printScheduleBrief(const Schedule& s) {
    std::cout << "  cost=" << s.cost << (s.isDummy ? "  [dummy]" : "") << "\n";
    std::cout << "  Q   : ";
    for (double v : s.Q) std::cout << std::setw(12) << v;
    std::cout << "\n  Qr  : ";
    for (double v : s.Qr) std::cout << std::setw(12) << v;

    std::cout << "\n  Y  : ";
    for (double v : s.Y) std::cout << std::setw(12) << v;
    std::cout << "\n  Yr  : ";
    for (double v : s.Yr) std::cout << std::setw(12) << v;

    std::cout << "\n  gam : ";
    for (double v : s.gamma) std::cout << std::setw(12) << v;
    std::cout << "\n  gamr: ";
    for (double v : s.gammar) std::cout << std::setw(12) << v;
    std::cout << "\n";
}

#ifndef BRANCHING

struct BranchConstraint{
    int k = -1;      // product index
    int n = -1;      // original column index in pool[k]
    int value = -1;  // 0 or 1
};

struct RestrictedBranchAndBoundResult {
    bool foundFeasible = false;
    double bestObjective = IloInfinity;
    std::vector<std::vector<double>> bestDeltaValues;
    int exploredNodes = 0;
};

struct ActiveColumnsFromCG {
    std::vector<std::vector<int>> activeCols;  // original column indices from pool
};

ActiveColumnsFromCG extractPositiveDeltaColumns(const Instance& ins, const std::vector<std::vector<Schedule>>& pool, const MasterResult& relaxedMaster, double tol = 1e-9)
{
    ActiveColumnsFromCG activeSet;
    activeSet.activeCols.resize(ins.K);

    for (int k = 0; k < ins.K; ++k) {
        for (int n = 0; n < pool[k].size(); ++n) {
            activeSet.activeCols[k].push_back(n);
        }
    }

    return activeSet;
}

bool isIntegralRestrictedMasterSolution(const MasterResult& master, const ActiveColumnsFromCG& activeSet,double tol = 1e-9)
{
    for (int k = 0; k < activeSet.activeCols.size(); ++k) {
        for (int originalCol : activeSet.activeCols[k]) {
            double v = master.deltaValues[k][originalCol];
            if (v > tol && v < 1.0 - tol) {
                return false;
            }
        }
    }
    return true;
}

bool findFractionalDeltaOnPositiveSupport(const MasterResult& master, const ActiveColumnsFromCG& activeSet, int& branchK, int& branchOriginalCol, double& branchValue, double tol = 1e-9)
{
    for (int k = 0; k < activeSet.activeCols.size(); ++k) {
        for (int originalCol : activeSet.activeCols[k]) {
            double v = master.deltaValues[k][originalCol];
            if (v > tol && v < 1.0 - tol) {
                branchK = k;
                branchOriginalCol = originalCol;
                branchValue = v;
                return true;
            }
        }
    }
    return false;
}

bool areBranchingConstraintsConsistent(const Instance& datasetInstance, const std::vector<BranchConstraint>& branchingConstraints)
{
    std::vector<std::vector<int>> fixed(datasetInstance.K);

    for (int k = 0; k < datasetInstance.K; ++k) {
        fixed[k].assign(10000, -1);
    }

    std::vector<int> chosenOne(datasetInstance.K, -1);

    for (const auto& bc : branchingConstraints) {
        if (bc.k < 0 || bc.n < 0 || (bc.value != 0 && bc.value != 1)) {
            return false;
        }

        if (bc.value == 1) {
            if (chosenOne[bc.k] != -1 && chosenOne[bc.k] != bc.n) {
                return false;
            }
            chosenOne[bc.k] = bc.n;
        }
    }

    for (std::size_t i = 0; i < branchingConstraints.size(); ++i) {
        for (std::size_t j = i + 1; j < branchingConstraints.size(); ++j) {
            if (branchingConstraints[i].k == branchingConstraints[j].k && branchingConstraints[i].n == branchingConstraints[j].n && branchingConstraints[i].value != branchingConstraints[j].value) {
                return false;
            }
        }
    }

    return true;
}

MasterResult solveSPPWithBranching(const Instance& datasetInstance, const std::vector<std::vector<Schedule>>& schedulesPool, const ActiveColumnsFromCG& activeSet, const std::vector<BranchConstraint>& branchingConstraints)
{
    IloEnv env;
    MasterResult res;

    try {
        if (!areBranchingConstraintsConsistent(datasetInstance, branchingConstraints)) {
            res.solved = false;
            env.end();
            return res;
        }

        IloModel model(env);

        std::vector<std::vector<IloNumVar>> deltaLocal(datasetInstance.K);
        IloExpr objective(env);

        for (int k = 0; k < datasetInstance.K; ++k) {
            const auto& cols = activeSet.activeCols[k];
            deltaLocal[k].reserve(cols.size());

            for (std::size_t j = 0; j < cols.size(); ++j) {
                int originalCol = cols[j];

                std::ostringstream name;
                name << "delta_k" << (k + 1) << "_col" << (originalCol + 1);

                deltaLocal[k].push_back(IloNumVar(env, 0.0, 1.0, ILOFLOAT, name.str().c_str()));
                objective += schedulesPool[k][originalCol].cost * deltaLocal[k].back();
            }
        }

        model.add(IloMinimize(env, objective));

        // choose one schedule per product
        IloRangeArray chooseOne(env);
        for (int k = 0; k < datasetInstance.K; ++k) {
            IloExpr lhs(env);
            for (std::size_t j = 0; j < deltaLocal[k].size(); ++j) {
                lhs += deltaLocal[k][j];
            }
            std::ostringstream cname;
            cname << "chooseOne_k" << (k + 1);

            chooseOne.add(IloRange(env, 1.0, lhs, 1.0, cname.str().c_str()));
            lhs.end();
        }
        model.add(chooseOne);

        // production capacity
        IloRangeArray prodCapConstr(env);
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < datasetInstance.K; ++k) {
                const auto& cols = activeSet.activeCols[k];
                for (std::size_t j = 0; j < cols.size(); ++j) {
                    int originalCol = cols[j];
                    lhs += schedulesPool[k][originalCol].prodCap[t] * deltaLocal[k][j];
                }
            }
            std::ostringstream cname;
            cname << "prodCap_t" << (t + 1);

            prodCapConstr.add(IloRange(env, -IloInfinity, lhs, datasetInstance.ct[t], cname.str().c_str()));
            lhs.end();
        }
        model.add(prodCapConstr);

        // rem capacity
        IloRangeArray remCapConstr(env);
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < datasetInstance.K; ++k) {
                const auto& cols = activeSet.activeCols[k];
                for (std::size_t j = 0; j < cols.size(); ++j) {
                    int originalCol = cols[j];
                    lhs += schedulesPool[k][originalCol].remCap[t] * deltaLocal[k][j];
                }
            }
            std::ostringstream cname;
            cname << "remCap_t" << (t + 1);

            remCapConstr.add(IloRange(env, -IloInfinity, lhs, datasetInstance.ctr[t], cname.str().c_str()));
            lhs.end();
        }
        model.add(remCapConstr);

        for (const auto& bc : branchingConstraints) {
            int k = bc.k;
            int originalCol = bc.n;

            int localIndex = -1;
            for (int j = 0; j < static_cast<int>(activeSet.activeCols[k].size()); ++j) {
                if (activeSet.activeCols[k][j] == originalCol) {
                    localIndex = j;
                    break;
                }
            }

            if (localIndex < 0) {
                // branch on index from outside support - invalid
                res.solved = false;
                objective.end();
                env.end();
                return res;
            }

            if (bc.value == 0) {
                model.add(deltaLocal[k][localIndex] == 0.0);
            }
            else {
                model.add(deltaLocal[k][localIndex] == 1.0);

                // if the column is 1, then all others must be 0
                for (int j = 0; j < static_cast<int>(activeSet.activeCols[k].size()); ++j) {
                    if (j == localIndex) continue;
                    model.add(deltaLocal[k][j] == 0.0);
                }
            }
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

        if (!cplex.solve()) {
            res.solved = false;
            objective.end();
            env.end();
            return res;
        }

        res.solved = true;
        res.objective = cplex.getObjValue();
        res.sigma.resize(datasetInstance.K, 0.0);
        res.piProd.resize(datasetInstance.T, 0.0);
        res.piRem.resize(datasetInstance.T, 0.0);

        res.deltaValues.resize(datasetInstance.K);
        for (int k = 0; k < datasetInstance.K; ++k) {
            res.deltaValues[k].assign(schedulesPool[k].size(), 0.0);
        }

        for (int k = 0; k < datasetInstance.K; ++k) {
            res.sigma[k] = cplex.getDual(chooseOne[k]);

            const auto& cols = activeSet.activeCols[k];
            for (std::size_t j = 0; j < cols.size(); ++j) {
                int originalCol = cols[j];
                res.deltaValues[k][originalCol] = cplex.getValue(deltaLocal[k][j]);
            }
        }

        for (int t = 0; t < datasetInstance.T; ++t) {
            res.piProd[t] = cplex.getDual(prodCapConstr[t]);
            res.piRem[t] = cplex.getDual(remCapConstr[t]);
        }

        objective.end();
        env.end();
        return res;
    }
    catch (...) {
        env.end();
        throw;
    }
}

void restrictedBranchAndBound(const Instance& datasetInstance, const std::vector<std::vector<Schedule>>& schedulesPool, const ActiveColumnsFromCG& activeSet,const std::vector<BranchConstraint>& currentConstraints, RestrictedBranchAndBoundResult& globalBest, double tol = 1e-9)
{
    MasterResult nodeResult = solveSPPWithBranching(datasetInstance, schedulesPool, activeSet, currentConstraints);

    globalBest.exploredNodes++;

    // infeasible node
    if (!nodeResult.solved) {
        return;
    }

    // bound pruning
    if (globalBest.foundFeasible && nodeResult.objective >= globalBest.bestObjective - tol) {
        return;
    }

    // integral node -> incumbent
    if (isIntegralRestrictedMasterSolution(nodeResult, activeSet, tol)) {
        globalBest.foundFeasible = true;
        globalBest.bestObjective = nodeResult.objective;
        globalBest.bestDeltaValues = nodeResult.deltaValues;
        return;
    }

    int branchK = -1;
    int branchOriginalCol = -1;
    double branchValue = 0.0;

    bool found = findFractionalDeltaOnPositiveSupport(nodeResult, activeSet, branchK,  branchOriginalCol, branchValue, tol);

    if (!found) {
        return;
    }

    // LEFT branch: delta(k,n) = 0
    {
        std::vector<BranchConstraint> leftConstraints = currentConstraints;
        leftConstraints.push_back({ branchK, branchOriginalCol, 0 });

        restrictedBranchAndBound(datasetInstance, schedulesPool, activeSet, leftConstraints, globalBest, tol);
    }

    // RIGHT branch: delta(k,n) = 1
    {
        std::vector<BranchConstraint> rightConstraints = currentConstraints;
        rightConstraints.push_back({ branchK, branchOriginalCol, 1 });

        restrictedBranchAndBound(datasetInstance, schedulesPool, activeSet, rightConstraints, globalBest, tol);
    }
}

RestrictedBranchAndBoundResult solveBranchAndBound(const Instance& datasetInstance, const std::vector<std::vector<Schedule>>& schedulesPool, const MasterResult& relaxedMasterFromCG, double tol = 1e-9){
    ActiveColumnsFromCG activeSet = extractPositiveDeltaColumns(datasetInstance, schedulesPool, relaxedMasterFromCG, tol);

    RestrictedBranchAndBoundResult result;
    std::vector<BranchConstraint> rootConstraints;

    for (int k = 0; k < datasetInstance.K; ++k) {
        if (activeSet.activeCols[k].empty()) {
            return result;
        }
    }

    restrictedBranchAndBound(datasetInstance, schedulesPool, activeSet, rootConstraints, result, tol);

    return result;
}
#endif

//Fixing
struct FixedSetupsFromCG {
    // -1 = no fixing
    //  0 = fix do 0
    //  1 = fix do 1
    std::vector<std::vector<int>> gammaFix;
    std::vector<std::vector<int>> gammarFix;
};

struct FixedFullResult {
    bool solved = false;
    bool optimal = false;

    double objective = 0.0;
    double bestBound = 0.0;
    double mipGap = 0.0;

    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> Qr;
    std::vector<std::vector<double>> Y;
    std::vector<std::vector<double>> Yr;
    std::vector<std::vector<double>> gamma;
    std::vector<std::vector<double>> gammar;
};

FixedSetupsFromCG extractSetupFixingsFromIntegralProducts(const Instance& datasetInstance, const MasterResult& master, const std::vector<std::vector<Schedule>>& schedulesPool, double tol = 1e-6){
    FixedSetupsFromCG fix;
    fix.gammaFix.assign(datasetInstance.K, std::vector<int>(datasetInstance.T, -1));
    fix.gammarFix.assign(datasetInstance.K, std::vector<int>(datasetInstance.T, -1));

    for (int k = 0; k < datasetInstance.K; ++k) {
        int chosenCol = -1;

        for (int n = 0; n < master.deltaValues[k].size(); ++n) {
            if (std::abs(master.deltaValues[k][n] - 1.0) <= tol) {
                chosenCol = n;
                break;
            }
        }

        // fix only when product is integral
        if (chosenCol < 0) {
            continue;
        }

        const Schedule& s = schedulesPool[k][chosenCol];

        for (int t = 0; t < datasetInstance.T; ++t) {
            fix.gammaFix[k][t] = (s.gamma[t] > 0.5 ? 1 : 0);
            fix.gammarFix[k][t] = (s.gammar[t] > 0.5 ? 1 : 0);
        }
    }

    return fix;
}

void printSetupFixingsFromCG(const Instance& datasetInstance, const FixedSetupsFromCG& fix){
    std::cout << "\n========== Setup fixings extracted from integral CG products ==========\n";

    for (int k = 0; k < datasetInstance.K; ++k) {
        bool anyFix = false;
        for (int t = 0; t < datasetInstance.T; ++t) {
            if (fix.gammaFix[k][t] != -1 || fix.gammarFix[k][t] != -1) {
                anyFix = true;
                break;
            }
        }

        if (!anyFix) {
            std::cout << "Product " << (k + 1) << ": no fixing\n";
            continue;
        }

        std::cout << "Product " << (k + 1) << ":\n";
        std::cout << "  gamma : ";
        for (int t = 0; t < datasetInstance.T; ++t) {
            std::cout << std::setw(6) << fix.gammaFix[k][t];
        }
        std::cout << "\n";

        std::cout << "  gammar: ";
        for (int t = 0; t < datasetInstance.T; ++t) {
            std::cout << std::setw(6) << fix.gammarFix[k][t];
        }
        std::cout << "\n";
    }

    std::cout << "======================================================================\n";
}

FixedFullResult solveOriginalModelWithFixings(const Instance& datasetInstance, const std::vector<std::vector<double>>& M, const FixedSetupsFromCG& fixings, bool useValidIneq = true, double timeLimitSeconds = 60.0){
    IloEnv env;
    FixedFullResult result;

    try {
        IloModel model(env);

        std::vector<std::vector<IloNumVar>> Q(datasetInstance.K), Qr(datasetInstance.K);
        std::vector<std::vector<IloNumVar>> Y(datasetInstance.K), Yr(datasetInstance.K);
        std::vector<std::vector<IloBoolVar>> gamma(datasetInstance.K), gammar(datasetInstance.K);

        for (int k = 0; k < datasetInstance.K; ++k) {
            Q[k].reserve(datasetInstance.T);
            Qr[k].reserve(datasetInstance.T);
            Y[k].reserve(datasetInstance.T);
            Yr[k].reserve(datasetInstance.T);
            gamma[k].reserve(datasetInstance.T);
            gammar[k].reserve(datasetInstance.T);

            for (int t = 0; t < datasetInstance.T; ++t) {
                Q[k].push_back(IloNumVar(env, 0.0, IloInfinity, ILOINT));
                Qr[k].push_back(IloNumVar(env, 0.0, IloInfinity, ILOINT));
                Y[k].push_back(IloNumVar(env, 0.0, IloInfinity, ILOINT));
                Yr[k].push_back(IloNumVar(env, 0.0, IloInfinity, ILOINT));
                gamma[k].push_back(IloBoolVar(env));
                gammar[k].push_back(IloBoolVar(env));
            }
        }

        IloExpr objective(env);

        for (int k = 0; k < datasetInstance.K; ++k) {
            for (int t = 0; t < datasetInstance.T; ++t) {
                objective += datasetInstance.hc[k] * Y[k][t];
                objective += datasetInstance.hcr[k] * Yr[k][t];
                objective += datasetInstance.pc[k] * Q[k][t];
                objective += datasetInstance.pcr[k] * Qr[k][t];
                objective += datasetInstance.sc[k] * gamma[k][t];
                objective += datasetInstance.scr[k] * gammar[k][t];
            }
        }

        model.add(IloMinimize(env, objective));

        // serviceables balance
        for (int k = 0; k < datasetInstance.K; ++k) {
            model.add(datasetInstance.y0[k] + Q[k][0] + Qr[k][0] == datasetInstance.d[k][0] + Y[k][0]);

            for (int t = 1; t < datasetInstance.T; ++t) {
                model.add(Y[k][t - 1] + Q[k][t] + Qr[k][t] == datasetInstance.d[k][t] + Y[k][t]);
            }
        }

        // recoverables balance
        for (int k = 0; k < datasetInstance.K; ++k) {
            model.add(datasetInstance.yr0[k] + datasetInstance.r[k][0] == Qr[k][0] + Yr[k][0]);

            for (int t = 1; t < datasetInstance.T; ++t) {
                model.add(Yr[k][t - 1] + datasetInstance.r[k][t] == Qr[k][t] + Yr[k][t]);
            }
        }

        // production capacity
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);

            for (int k = 0; k < datasetInstance.K; ++k) {
                lhs += datasetInstance.tp[k] * Q[k][t] + datasetInstance.ts[k] * gamma[k][t];
            }

            model.add(lhs <= datasetInstance.ct[t]);
            lhs.end();
        }

        // remanufacturing capacity
        for (int t = 0; t < datasetInstance.T; ++t) {
            IloExpr lhs(env);

            for (int k = 0; k < datasetInstance.K; ++k) {
                lhs += datasetInstance.tpr[k] * Qr[k][t] + datasetInstance.tsr[k] * gammar[k][t];
            }

            model.add(lhs <= datasetInstance.ctr[t]);
            lhs.end();
        }

        // linking
        for (int k = 0; k < datasetInstance.K; ++k) {
            for (int t = 0; t < datasetInstance.T; ++t) {
                model.add(Q[k][t] <= M[k][t] * gamma[k][t]);
                model.add(Qr[k][t] <= M[k][t] * gammar[k][t]);
            }
        }

        // optional valid inequalities
        if (useValidIneq) {
            for (int k = 0; k < datasetInstance.K; ++k) {
                // VI (18)
                for (int t = 0; t <= datasetInstance.T - 2; ++t) {
                    for (int p = 1; p <= datasetInstance.T - 1 - t; ++p) {
                        IloExpr rhs(env);

                        for (int s = t + 1; s <= t + p; ++s) {
                            rhs += datasetInstance.d[k][s];
                            rhs -= M[k][s] * gamma[k][s];
                            rhs -= M[k][s] * gammar[k][s];
                        }

                        model.add(Y[k][t] >= rhs);
                        rhs.end();
                    }
                }

                // VI (19)
                for (int t = 1; t < datasetInstance.T; ++t) {
                    for (int p = 1; p <= t; ++p) {
                        IloExpr rhs(env);

                        for (int s = t - p; s <= t; ++s) {
                            rhs += datasetInstance.r[k][s];
                            rhs -= M[k][s] * gammar[k][s];
                        }

                        model.add(Yr[k][t] >= rhs);
                        rhs.end();
                    }
                }
            }
        }

        // ===== FIXINGS FROM CG =====
        for (int k = 0; k < datasetInstance.K; ++k) {
            for (int t = 0; t < datasetInstance.T; ++t) {
                if (fixings.gammaFix[k][t] != -1) {
                    model.add(gamma[k][t] == fixings.gammaFix[k][t]);
                }
                if (fixings.gammarFix[k][t] != -1) {
                    model.add(gammar[k][t] == fixings.gammarFix[k][t]);
                }
            }
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::TimeLimit, timeLimitSeconds);

        if (!cplex.solve()) {
            objective.end();
            env.end();
            return result;
        }

        result.solved = true;
        result.objective = cplex.getObjValue();
        result.bestBound = cplex.getBestObjValue();
        result.mipGap = cplex.getMIPRelativeGap();
        result.optimal = (cplex.getStatus() == IloAlgorithm::Optimal);

        result.Q.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));
        result.Qr.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));
        result.Y.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));
        result.Yr.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));
        result.gamma.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));
        result.gammar.assign(datasetInstance.K, std::vector<double>(datasetInstance.T, 0.0));

        for (int k = 0; k < datasetInstance.K; ++k) {
            for (int t = 0; t < datasetInstance.T; ++t) {
                result.Q[k][t] = cplex.getValue(Q[k][t]);
                result.Qr[k][t] = cplex.getValue(Qr[k][t]);
                result.Y[k][t] = cplex.getValue(Y[k][t]);
                result.Yr[k][t] = cplex.getValue(Yr[k][t]);
                result.gamma[k][t] = cplex.getValue(gamma[k][t]);
                result.gammar[k][t] = cplex.getValue(gammar[k][t]);
            }
        }

        objective.end();
        env.end();
        return result;
    }
    catch (...) {
        env.end();
        throw;
    }
}

void printFixedFullModelResultSummary(const FixedFullResult& res)
{
    std::cout << "\n========== Full Model with setup fixings from CG ==========\n";

    if (!res.solved) {
        std::cout << "No feasible solution found.\n";
        std::cout << "==========================================================\n";
        return;
    }

    std::cout << "Objective  = " << res.objective << "\n";
    std::cout << "BestBound  = " << res.bestBound << "\n";
    std::cout << "MIP Gap    = " << 100.0 * res.mipGap << "%\n";
    std::cout << "Status     = " << (res.optimal ? "Optimal" : "Time-limited / non-proven") << "\n";
    std::cout << "==========================================================\n";
}

int main() {
    try {
        const Instance ins = buildReferenceInstance();
        const auto M = computeBigM(ins);

        // column pool: one list of schedules per product
        std::vector<std::vector<Schedule>> pool(ins.K);
        for (int k = 0; k < ins.K; ++k) {
            pool[k].push_back(makeDummySchedule(ins, k));
        }

        const int maxIterations = 100;
        const double rcTolerance = -1e-6;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "============================================================\n";
        std::cout << "Column generation for CLSP-RM-SS (SPP + SLULSP-RM pricing)\n";
        std::cout << "============================================================\n\n";

        MasterResult master;
        int iter = 0;
        for (; iter < maxIterations; ++iter) {
            master = solveRelaxedMasterProblemSPP(ins, pool);

            std::cout << "Iteration " << (iter + 1) << "\n";
            std::cout << "  RMP LP objective = " << master.objective << "\n";
            std::cout << "  sigma: ";
            for (double v : master.sigma) std::cout << std::setw(16) << v;
            std::cout << "\n  piProd: ";
            for (double v : master.piProd) std::cout << std::setw(16) << v;
            std::cout << "\n  piRem : ";
            for (double v : master.piRem) std::cout << std::setw(16) << v;
            std::cout << "\n";

            // Generate new columns for each product
            bool addedAnyColumn = false;
            for (int k = 0; k < ins.K; ++k) {
                PricingResult pr = solvePricingSubproblem(ins, M, k, master.sigma[k], master.piProd, master.piRem, rcTolerance);
                std::cout << "  Product " << (k + 1) << " pricing rc = " << pr.reducedCost;

                if (pr.improving) {
                    pool[k].push_back(pr.column);
                    addedAnyColumn = true;
                    std::cout << "  -> added column #" << pool[k].size();
                }
                std::cout << "\n";
            }

            std::cout << "  Columns per product: ";
            for (int k = 0; k < ins.K; ++k) std::cout << std::setw(8) << pool[k].size();
            std::cout << "\n\n";

            if (!addedAnyColumn) {
                std::cout << "No improving columns found. CG terminates.\n\n";
                break;
            }
        }

        if (iter == maxIterations) {
            std::cout << "Reached maxIterations = " << maxIterations << ".\n\n";
        }

        // Lower bound from relaxed master
        std::cout << "Lower bound (from relaxing problem) " << master.objective << "\n\n";

        // Upper bound from integer master
        double integerMasterObj = std::numeric_limits<double>::quiet_NaN();
        auto chosen = solveIntegerMasterProblemSPP(ins, pool, integerMasterObj);
        std::cout << "Upper Bound (from solving master problem over generated columns): " << integerMasterObj << "\n\n"; //?

        // Printing delta values

        std::cout << "========== Relaxed SPP deltas (found schedules) ==========\n";
        for (int i = 0; i < master.deltaValues.size(); i++) {
            std::cout << "\nProduct " << i + 1 << " has " << master.deltaValues[i].size() << " schedules\n";

            for (int j = 0; j < master.deltaValues[i].size(); j++) {
                if (master.deltaValues[i][j] <= 0.001) continue;
                std::cout << "\tcol: " << j+1 <<  "\tdelta: " << master.deltaValues[i][j] << "\t cost: " << pool[i][j].cost << "\n";
            }
        }

        std::cout << "========== Integral deltas (found schedules) ==========\n";
        for (int i = 0; i < chosen.size(); i++) {
            std::cout << "\nProduct " << i + 1 << " has " << chosen[i].size() << " schedules\n";

            for (int j = 0; j < chosen[i].size(); j++) {
                //if (chosen[i][j] <= 0.001) continue;
                std::cout << "\tcol: " << j + 1 << "\tdelta: " << chosen[i][j] << "\t cost: " << pool[i][j].cost << "\n";
            }
        }

#ifdef PRINT_SCHEDULES
        std::cout << "\n========== Printing schedules ==========\n";
        for (int i = 0; i < pool.size(); i++) {
            for (int j = 0; j < pool[i].size(); j++) {
                std::cout << "Product " << i+1 << " - schedule " << j+1 << ":\n";
                printScheduleBrief(pool[i][j]);
            }
        }
#endif
#ifdef BRANCHING
        std::cout << "\n========== B&B only on >0 deltas from CG ==========\n";
        RestrictedBranchAndBoundResult restrictedBB = solveBranchAndBound(ins, pool, master);

        if (!restrictedBB.foundFeasible) {
            std::cout << "  No feasible integral SPP solution exists on this positive-delta support.\n";
        }
        else {
            std::cout << "  Best objective = " << restrictedBB.bestObjective << "\n";
            std::cout << "  Explored nodes = " << restrictedBB.exploredNodes << "\n";

            for (int k = 0; k < ins.K; ++k) {
                for (std::size_t n = 0; n < restrictedBB.bestDeltaValues[k].size(); ++n) {
                    if (restrictedBB.bestDeltaValues[k][n] > 0.5) {
                        std::cout << "  Product " << (k + 1) << " -> column #" << (n + 1) << "  cost = " << pool[k][n].cost << "\n";
                    }
                }
            }
        }
#endif

        FixedSetupsFromCG fixings = extractSetupFixingsFromIntegralProducts(ins, master, pool);
        printSetupFixingsFromCG(ins, fixings);

        FixedFullResult fixedMilp = solveOriginalModelWithFixings(ins, M, fixings, true, 60.0);
        printFixedFullModelResultSummary(fixedMilp);

        return 0;
    }
    catch (const IloException& e) {
        std::cerr << "Concert/CPLEX exception: " << e << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Standard exception: " << e.what() << "\n";
    }
    catch (...) {
        std::cerr << "Unknown exception.\n";
    }
    return 1;
}
