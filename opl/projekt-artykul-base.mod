/*********************************************
 * OPL 22.1.2.0 Model
 * Author: Kacper79
 * Creation Date: 28 Mar 2026 at 14:20:49
 *********************************************/

// Indicies - ilosc uzywanych produktow i okresow czasowych
int K = ...;                    // number of products
int T = ...;                    // number of periods
range Products = 1..K;
range Periods  = 1..T;

// Parametry kosztu
float hc[Products]   = ...;     // holding cost of serviceables
float hcr[Products]  = ...;     // holding cost of recoverables
float pc[Products]   = ...;     // production cost per unit
float pcr[Products]  = ...;     // remanufacturing cost per unit
float sc[Products]   = ...;     // setup cost for production
float scr[Products]  = ...;     // setup cost for remanufacturing

// Parametry kosztu jednostkowego produkcji i capacity dostepnego
float tp[Products]   = ...;     // production time per unit
float tpr[Products]  = ...;     // remanufacturing time per unit
float ts[Products]   = ...;     // setup time for production
float tsr[Products]  = ...;     // setup time for remanufacturing

// Pojemnosci - capacity
float ct[Periods] = ...;   // manufacturing capacity c_t
float ctr[Periods] = ...;   // remanufacturing capacity cr_t

// Demand / returns
float d[Products][Periods] = ...;   // demand d_kt
float r[Products][Periods] = ...;   // returns r_kt

// Initial inventories (default 0)
float y0[Products]  = ...;      // initial serviceable inventory Y_{k,0}
float yr0[Products] = ...;      // initial recoverable inventory Y^r_{k,0}

// Optional switch for valid inequalities (18)-(19)
int useValidIneq = ...;         // 0 = off, 1 = on

// Zawezanie wartosci M
float M[k in Products][t in Periods] = maxl(0, maxl(ct[t] - ts[k], ctr[t] - tsr[k])) / minl(tp[k], tpr[k]);

// Zmienne decyzyjne
dvar float+ Q [Products][Periods];    // production quantity
dvar float+ Qr[Products][Periods];    // remanufacturing quantity
dvar float+ Y [Products][Periods];    // end inventory of serviceables
dvar float+ Yr[Products][Periods];    // end inventory of recoverables

// Binarne wartosci czy ponosimy koszt ustawiania (czy one sie dzieje)
dvar boolean gamma [Products][Periods];   // production setup
dvar boolean gammar[Products][Periods];   // remanufacturing setup

// Cele optymalizacji
dexpr float holdingCost = sum(k in Products, t in Periods) (hc[k]  * Y[k][t] + hcr[k] * Yr[k][t]);
dexpr float variableCost = sum(k in Products, t in Periods) (pc[k]  * Q[k][t] + pcr[k] * Qr[k][t]);
dexpr float setupCost = sum(k in Products, t in Periods) (sc[k]  * gamma[k][t] + scr[k] * gammar[k][t]);

minimize holdingCost + variableCost + setupCost; //koszt calkowity to suma kosztow

subject to {

  // (2) Balans zatowarowania (ilosc towarow)
  // wiaze stan magazynowy poprzedni [t-1] + produkcje -> zapotrzebowanie + resztki na stan magazynowy
  forall(k in Products)
    y0[k] + Q[k][1] + Qr[k][1] == d[k][1] + Y[k][1]; //przypadek bazowy

  forall(k in Products, t in 2..T)
    Y[k][t-1] + Q[k][t] + Qr[k][t] == d[k][t] + Y[k][t];

  // (3) Balans recyklingu
  // wiaze liczbe zwrotow produktow w danym czasie + stan magazynu zwroconych -> ilosc zrecyklingowanych rzeczy + reszta stanu magazynu zwroconych
  forall(k in Products)
    yr0[k] + r[k][1] == Qr[k][1] + Yr[k][1]; //przypadek bazowy

  forall(k in Products, t in 2..T)
    Yr[k][t-1] + r[k][t] == Qr[k][t] + Yr[k][t];

  // (4) Produkcja roznych produktow do limitu capacity
  // wiaze czas produkcji towaru i ilosc towarow + czas przezbrajania maszyn (nie zawsze - zmienna binarna) -> do limitu capacity produkcji
  forall(t in Periods)
    sum(k in Products) (tp[k] * Q[k][t] + ts[k] * gamma[k][t]) <= ct[t];

  // (5) Recykling zwroconych produktow do limit capacity recyklingu
  // wiaze czas recyklingu towaru i ilosc towarow + czas przezbrajania maszyn (nie zawsze - zmienna binarna) -> do limitu capacity produkcji
  forall(t in Periods)
    sum(k in Products) (tpr[k] * Qr[k][t] + tsr[k] * gammar[k][t]) <= ctr[t];

  // (6) Wiaze produkcje towarow z zmienna binarna czy konieczne jest przezbrajanie
  forall(k in Products, t in Periods)
    Q[k][t] <= M[k][t] * gamma[k][t];

  // (7) Wiaze recykling towarow z zmienna binarna czy konieczne jest przezbrajanie
  forall(k in Products, t in Periods)
    Qr[k][t] <= M[k][t] * gammar[k][t];

  // (18) Opcjonalne valid inequalities dla produkcji do przyspieszania obliczen
  forall(k in Products, t in 1..T-1, p in 1..T-t : useValidIneq == 1)
    Y[k][t] >= sum(s in t+1..t+p) d[k][s] - sum(s in t+1..t+p) M[k][s] * (gamma[k][s] + gammar[k][s]);

  // (19) Opcjonalne valid inequalities dla recyklingu do przyspieszania obliczen
  forall(k in Products, t in 2..T, p in 1..t-1 : useValidIneq == 1)
    Yr[k][t] >= sum(s in t-p..t) r[k][s] - sum(s in t-p..t) M[k][s] * gammar[k][s];
}

execute DISPLAY_SUMMARY {
  writeln("----------------------------------------------");
  writeln("Model CLSP-RM-SS (bez harmonogramowania)");
  writeln("Objective value = ", cplex.getObjValue());
  writeln("Holding cost    = ", holdingCost.solutionValue);
  writeln("Variable cost   = ", variableCost.solutionValue);
  writeln("Setup cost      = ", setupCost.solutionValue);
  writeln("----------------------------------------------");
}

execute DISPLAY_TABLES {
  writeln("\nQ (production):");
  for(var k in Products) {
    write("k=", k, ": ");
    for(var t in Periods) write(Q[k][t].solutionValue, " ");
    writeln();
  }

  writeln("\nQr (remanufacturing):");
  for(var k2 in Products) {
    write("k=", k2, ": ");
    for(var t2 in Periods) write(Qr[k2][t2].solutionValue, " ");
    writeln();
  }

  writeln("\nY (serviceable inventory):");
  for(var k3 in Products) {
    write("k=", k3, ": ");
    for(var t3 in Periods) write(Y[k3][t3].solutionValue, " ");
    writeln();
  }

  writeln("\nYr (recoverable inventory):");
  for(var k4 in Products) {
    write("k=", k4, ": ");
    for(var t4 in Periods) write(Yr[k4][t4].solutionValue, " ");
    writeln();
  }

  writeln("\ngamma (prod setup):");
  for(var k5 in Products) {
    write("k=", k5, ": ");
    for(var t5 in Periods) write(gamma[k5][t5].solutionValue, " ");
    writeln();
  }

  writeln("\ngammar (reman setup):");
  for(var k6 in Products) {
    write("k=", k6, ": ");
    for(var t6 in Periods) write(gammar[k6][t6].solutionValue, " ");
    writeln();
  }
}
 