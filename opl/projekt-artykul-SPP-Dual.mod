/*********************************************
 * OPL 22.1.2.0 Model
 * Author: Kacper79
 * Creation Date: 7 May 2026 at 12:03:20
 *********************************************/

int K = ...;
int T = ...;
int maxN = ...;

range Products = 1..K;
range Periods  = 1..T;
range Cols     = 1..maxN;

int N[Products] = ...;

float ct[Periods]  = ...;
float ctr[Periods] = ...;

float nc[Products][Cols] = ...;
float prodUse[Products][Cols][Periods] = ...;
float remUse[Products][Cols][Periods]  = ...;

// sigma is free because it corresponds to equality choose-one constraints.
dvar float sigma[Products];

// Capacity constraints are <= in a minimization primal, hence these duals are <= 0.
dvar float piProd[Periods];
dvar float piRem[Periods];

maximize
  sum(k in Products) sigma[k]
  + sum(t in Periods) ct[t]  * piProd[t]
  + sum(t in Periods) ctr[t] * piRem[t];

subject to {
  forall(t in Periods) {
    piProd[t] <= 0;
    piRem[t]  <= 0;
  }

  // One dual constraint for every active column of the primal master.
  forall(k in Products, n in Cols : n <= N[k])
    sigma[k]
    + sum(t in Periods) prodUse[k][n][t] * piProd[t]
    + sum(t in Periods) remUse[k][n][t]  * piRem[t]
    <= nc[k][n];
}
