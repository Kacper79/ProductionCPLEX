#include <ProductionCPLEX.h>

ILOSTLBEGIN

struct Instance {
    int K = 0;
    int T = 0;

    std::vector<double> hc, hcr;
    std::vector<double> pc, pcr;
    std::vector<double> sc, scr;

    std::vector<double> tp, tpr;
    std::vector<double> ts, tsr;

    std::vector<double> ct, ctr;

    std::vector<std::vector<double>> d;
    std::vector<std::vector<double>> r;

    std::vector<double> y0, yr0;
    bool useValidIneq = true;
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

std::vector<std::vector<double>> computeBigM(const Instance& ins) {
    std::vector<std::vector<double>> M(ins.K, std::vector<double>(ins.T, 0.0));
    for (int k = 0; k < ins.K; ++k) {
        for (int t = 0; t < ins.T; ++t) {
            const double numerator = max2(0.0, max2(ins.ct[t] - ins.ts[k], ins.ctr[t] - ins.tsr[k]));
            const double denominator = min2(ins.tp[k], ins.tpr[k]);
            if (denominator <= 0.0) {
                throw std::runtime_error("Dont divide by zero in calc M");
            }
            M[k][t] = numerator / denominator;
        }
    }
    return M;
}

Schedule makeDummySchedule(const Instance& ins, int k, double penalty = 1.0e6) {
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

PricingResult solvePricingSubproblem( const Instance& ins, const std::vector<std::vector<double>>& M, int k, double sigma, const std::vector<double>& piProd, const std::vector<double>& piRem, double rcTolerance = -1e-6) {
    IloEnv env;
    PricingResult result;

    try {
        IloModel model(env);

        std::vector<IloNumVar> Q(ins.T), Qr(ins.T), Y(ins.T), Yr(ins.T);
        std::vector<IloNumVar> gamma(ins.T), gammar(ins.T);

        for (int t = 0; t < ins.T; ++t) {
            Q[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            Qr[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            Y[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            Yr[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
            gamma[t] = IloNumVar(env, 0.0, 1.0, ILOINT);
            gammar[t] = IloNumVar(env, 0.0, 1.0, ILOINT);
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
        cplex.setParam(IloCplex::Param::MIP::Display, 0);

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
                objective += schedulesPool[k][n].cost * delta[k].back(); // min Z = sum_k sum_n_from_S nc_kn * delta_kn gdzie S to zbior schedules
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
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
        //cplex.setParam(IloCplex::Param::LPMethod, IloCplex::Primal);

        if (!cplex.solve()) {
            throw std::runtime_error("Restricted master LP was not solved successfully.");
        }

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

std::vector<std::vector<double>> solveIntegerMasterProblemSPP(const Instance& datasetInstance,const std::vector<std::vector<Schedule>>& schedulesPool, double& objectiveValue) {
    IloEnv env;

    try {
        IloModel model(env);
        std::vector<std::vector<IloBoolVar>> delta(datasetInstance.K);
        IloExpr objective(env);

        for (int k = 0; k < datasetInstance.K; ++k) {
            delta[k].reserve(static_cast<std::size_t>(schedulesPool[k].size()));
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) {
                delta[k].push_back(IloBoolVar(env));
                objective += schedulesPool[k][n].cost * delta[k].back(); // min Z = sum_k sum_n_from_S nc_kn * delta_kn gdzie S to zbior schedules
            }
        }
        model.add(IloMinimize(env, objective));

        for (int k = 0; k < datasetInstance.K; ++k) {
            IloExpr lhs(env);
            for (std::size_t n = 0; n < schedulesPool[k].size(); ++n) lhs += delta[k][n];
            model.add(lhs == 1.0);
            lhs.end();
        }

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

void printScheduleBrief(const Schedule& s) {
    std::cout << "  cost=" << s.cost << (s.isDummy ? "  [dummy]" : "") << "\n";
    std::cout << "  Q   : ";
    for (double v : s.Q) std::cout << std::setw(8) << v;
    std::cout << "\n  Qr  : ";
    for (double v : s.Qr) std::cout << std::setw(8) << v;
    std::cout << "\n  gam : ";
    for (double v : s.gamma) std::cout << std::setw(8) << v;
    std::cout << "\n  gamr: ";
    for (double v : s.gammar) std::cout << std::setw(8) << v;
    std::cout << "\n";
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
            for (int k = 0; k < ins.K; ++k) std::cout << std::setw(6) << pool[k].size();
            std::cout << "\n\n";

            if (!addedAnyColumn) {
                std::cout << "No improving columns found. CG terminates.\n\n";
                break;
            }
        }

        if (iter == maxIterations) {
            std::cout << "Reached maxIterations = " << maxIterations << ".\n\n";
        }

        std::cout << "Lower bound (from relaxing problem) " << master.objective << "\n\n";

        double integerMasterObj = std::numeric_limits<double>::quiet_NaN();
        auto chosen = solveIntegerMasterProblemSPP(ins, pool, integerMasterObj);
        std::cout << "Upper Bound (from solving master problem over generated columns): " << integerMasterObj << "\n\n";

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
