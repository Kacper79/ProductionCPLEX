// ProductionBaseCPLEX.cpp : Defines the entry point for the application.
//

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

    bool useValidIneq = true;
};

Instance buildReferenceInstance() {
    Instance ins;

    // Data 1-1 from the article should yield expected result 9620
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

double max2(double a, double b) { return (a > b) ? a : b; }
double min2(double a, double b) { return (a < b) ? a : b; }

// Maximum possible production quantity given available capacity
std::vector<std::vector<double>> computeBigM(const Instance& ins) {
    std::vector<std::vector<double>> M(ins.K, std::vector<double>(ins.T, 0.0));

    for (int k = 0; k < ins.K; ++k) {
        for (int t = 0; t < ins.T; ++t) {
            // M[k][t] = max(0, max(ct[t] - ts[k], ctr[t] - tsr[k])) / min(tp[k], tpr[k])
            const double numerator = max2(0.0, max2(ins.ct[t] - ins.ts[k], ins.ctr[t] - ins.tsr[k]));
            const double denominator = min2(ins.tp[k], ins.tpr[k]);

            if (denominator <= 0.000001) {
                throw std::runtime_error("tp/tpr parameters must be > 0");
            }
            M[k][t] = numerator / denominator;
        }
    }
    return M;
}

int main() {
    IloEnv env;
    try {
        const Instance ins = buildReferenceInstance();
        const auto M = computeBigM(ins);

        IloModel model(env);

        // Decision variables
        std::vector<std::vector<IloNumVar>> Q(ins.K, std::vector<IloNumVar>(ins.T));
        std::vector<std::vector<IloNumVar>> Qr(ins.K, std::vector<IloNumVar>(ins.T));
        std::vector<std::vector<IloNumVar>> Y(ins.K, std::vector<IloNumVar>(ins.T));
        std::vector<std::vector<IloNumVar>> Yr(ins.K, std::vector<IloNumVar>(ins.T));
        std::vector<std::vector<IloBoolVar>> gamma(ins.K, std::vector<IloBoolVar>(ins.T));
        std::vector<std::vector<IloBoolVar>> gammar(ins.K, std::vector<IloBoolVar>(ins.T));

        for (int k = 0; k < ins.K; ++k) {
            for (int t = 0; t < ins.T; ++t) {
                Q[k][t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
                Qr[k][t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
                Y[k][t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
                Yr[k][t] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
                gamma[k][t] = IloBoolVar(env);
                gammar[k][t] = IloBoolVar(env);
            }
        }

        // Objective function (1)
        IloExpr objective(env);
        IloExpr holdingCost(env);
        IloExpr variableCost(env);
        IloExpr setupCost(env);

        for (int k = 0; k < ins.K; ++k) {
            for (int t = 0; t < ins.T; ++t) {
                holdingCost += ins.hc[k] * Y[k][t] + ins.hcr[k] * Yr[k][t];
                variableCost += ins.pc[k] * Q[k][t] + ins.pcr[k] * Qr[k][t];
                setupCost += ins.sc[k] * gamma[k][t] + ins.scr[k] * gammar[k][t];
            }
        }
        objective += holdingCost + variableCost + setupCost;
        model.add(IloMinimize(env, objective));

        // Constraints

        // (2) Balance serviceables
        for (int k = 0; k < ins.K; ++k) {
            // t = 0 (i.e. period 1 in article notation)
            {
                IloExpr lhs(env);
                lhs += ins.y0[k] + Q[k][0] + Qr[k][0];
                model.add(lhs == ins.d[k][0] + Y[k][0]);
                lhs.end();
            }
            // t = 1..T-1
            for (int t = 1; t < ins.T; ++t) {
                IloExpr lhs(env);
                lhs += Y[k][t - 1] + Q[k][t] + Qr[k][t];
                model.add(lhs == ins.d[k][t] + Y[k][t]);
                lhs.end();
            }
        }

        // (3) Balance recoverables
        for (int k = 0; k < ins.K; ++k) {
            {
                IloExpr lhs(env);
                lhs += ins.yr0[k] + ins.r[k][0];
                model.add(lhs == Qr[k][0] + Yr[k][0]);
                lhs.end();
            }
            for (int t = 1; t < ins.T; ++t) {
                IloExpr lhs(env);
                lhs += Yr[k][t - 1] + ins.r[k][t];
                model.add(lhs == Qr[k][t] + Yr[k][t]);
                lhs.end();
            }
        }

        // (4) Capacity production
        for (int t = 0; t < ins.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < ins.K; ++k) {
                lhs += ins.tp[k] * Q[k][t] + ins.ts[k] * gamma[k][t];
            }
            model.add(lhs <= ins.ct[t]);
            lhs.end();
        }

        // (5) Capacity remanufacturing
        for (int t = 0; t < ins.T; ++t) {
            IloExpr lhs(env);
            for (int k = 0; k < ins.K; ++k) {
                lhs += ins.tpr[k] * Qr[k][t] + ins.tsr[k] * gammar[k][t];
            }
            model.add(lhs <= ins.ctr[t]);
            lhs.end();
        }

        // (6) Linking for production
        for (int k = 0; k < ins.K; ++k) {
            for (int t = 0; t < ins.T; ++t) {
                model.add(Q[k][t] <= M[k][t] * gamma[k][t]);
            }
        }

        // (7) Linking for remanufacturing
        for (int k = 0; k < ins.K; ++k) {
            for (int t = 0; t < ins.T; ++t) {
                model.add(Qr[k][t] <= M[k][t] * gammar[k][t]);
            }
        }

        // (18)-(19) optional valid inequalities
        if (ins.useValidIneq) {
            // Ensures enough inventory for demand if no setups
            // (18) Y[k][t] >= sum_{s=t+1}^{t+p} d[k][s] - sum_{s=t+1}^{t+p} M[k][s] * (gamma[k][s] + gammar[k][s])
            for (int k = 0; k < ins.K; ++k) {
                for (int t = 0; t <= ins.T - 2; ++t) {
                    for (int p = 1; p <= ins.T - 1 - t; ++p) {
                        IloExpr rhs(env);
                        for (int s = t + 1; s <= t + p; ++s) {
                            rhs += ins.d[k][s];
                            rhs -= M[k][s] * gamma[k][s];
                            rhs -= M[k][s] * gammar[k][s];
                        }
                        model.add(Y[k][t] >= rhs);
                        rhs.end();
                    }
                }
            }

            // Ensures enough stock for returns if no reman setups
            // (19) Yr[k][t] >= sum_{s=t-p}^{t} r[k][s] - sum_{s=t-p}^{t} M[k][s] * gammar[k][s]
            for (int k = 0; k < ins.K; ++k) {
                for (int t = 1; t < ins.T; ++t) {
                    for (int p = 1; p <= t; ++p) {
                        IloExpr rhs(env);
                        for (int s = t - p; s <= t; ++s) {
                            rhs += ins.r[k][s];
                            rhs -= M[k][s] * gammar[k][s];
                        }
                        model.add(Yr[k][t] >= rhs);
                        rhs.end();
                    }
                }
            }
        }

        // Solve
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Threads, 1);

        if (!cplex.solve()) {
            std::cerr << "CPLEX has not found any solutions with status " << cplex.getStatus() << "\n";
            objective.end();
            holdingCost.end();
            variableCost.end();
            setupCost.end();
            env.end();
            return 1;
        }

        // Results report
        std::cout << std::fixed << std::setprecision(2) << "\n\n";
        std::cout << "-----------------------------------------------\n";
        std::cout << "Model CLSP-RM-SS (C++, no schedules)\n";
        std::cout << "Status          = " << cplex.getStatus() << "\n";
        std::cout << "Objective value = " << cplex.getObjValue() << "\n";
        std::cout << "Holding cost    = " << cplex.getValue(holdingCost) << "\n";
        std::cout << "Variable cost   = " << cplex.getValue(variableCost) << "\n";
        std::cout << "Setup cost      = " << cplex.getValue(setupCost) << "\n";
        std::cout << "-----------------------------------------------\n\n";

        // Helper lambdas for printing results
        auto printNumMatrix = [&](const std::string& name, const std::vector<std::vector<IloNumVar>>& x) {
            std::cout << name << ":\n";

            for (int k = 0; k < ins.K; ++k) {
                std::cout << "k=" << (k + 1) << ": ";

                for (int t = 0; t < ins.T; ++t) {
                    std::cout << std::setw(12) << cplex.getValue(x[k][t]);
                }

                std::cout << "\n";
            }
            std::cout << "\n";
        };

        auto printBoolMatrix = [&](const std::string& name, const std::vector<std::vector<IloBoolVar>>& x) {
            std::cout << name << ":\n";

            for (int k = 0; k < ins.K; ++k) {
                std::cout << "k=" << (k + 1) << ": ";

                for (int t = 0; t < ins.T; ++t) {
                    std::cout << std::setw(12) << cplex.getValue(x[k][t]);
                }

                std::cout << "\n";
            }
            std::cout << "\n";
        };

        // Print decision variables
        printNumMatrix("Q   (production)", Q);
        printNumMatrix("Qr  (remanufacturing)", Qr);
        printNumMatrix("Y   (serviceable inventory)", Y);
        printNumMatrix("Yr  (recoverable inventory)", Yr);
        printBoolMatrix("gamma  (prod setup)", gamma);
        printBoolMatrix("gammar (reman setup)", gammar);

        // Display Big-M - maximum possible production quantity per period (6) (7)
        std::cout << "M values:\n";

        for (int k = 0; k < ins.K; ++k) {
            std::cout << "k=" << (k + 1) << ": ";

            for (int t = 0; t < ins.T; ++t) {
                std::cout << std::setw(8) << M[k][t];
            }

            std::cout << "\n";
        }
        std::cout << "\n";

        // Capacity usage vs available (4) (5)
        std::cout << "Capacity usage:\n";

        for (int t = 0; t < ins.T; ++t) {
            double prodUsed = 0.0;
            double remUsed = 0.0;

            for (int k = 0; k < ins.K; ++k) {
                prodUsed += ins.tp[k] * cplex.getValue(Q[k][t]) + ins.ts[k] * cplex.getValue(gamma[k][t]);
                remUsed += ins.tpr[k] * cplex.getValue(Qr[k][t]) + ins.tsr[k] * cplex.getValue(gammar[k][t]);
            }

            std::cout << "t=" << (t + 1) << "  prod=" << prodUsed << "/" << ins.ct[t] << "  rem=" << remUsed << "/" << ins.ctr[t] << "\n";
        }

        // Free CPLEX resources
        objective.end();
        holdingCost.end();
        variableCost.end();
        setupCost.end();
        env.end();

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

    env.end();
    return 1;
}