/*********************************************
 * OPL 22.1.2.0 Model
 * Author: Kacper79
 * Creation Date: 7 May 2026 at 12:04:36
 *********************************************/

int K = ...;
int T = ...;
int k0 = ...;

range Products = 1..K;
range Periods  = 1..T;

float d[Products][Periods] = ...;
float r[Products][Periods] = ...;

float hc[Products]  = ...;
float hcr[Products] = ...;
float pc[Products]  = ...;
float pcr[Products] = ...;
float sc[Products]  = ...;
float scr[Products] = ...;

float tp[Products]  = ...;
float tpr[Products] = ...;
float ts[Products]  = ...;
float tsr[Products] = ...;

float y0[Products]  = ...;
float yr0[Products] = ...;

float M[Products][Periods] = ...;

float sigma = ...;
float piProd[Periods] = ...;
float piRem[Periods]  = ...;

// Integer quantities: convenient for the numerical example and close to your C++ variant.
dvar int+ Q[Periods];
dvar int+ Qr[Periods];
dvar int+ Y[Periods];
dvar int+ Yr[Periods];

dvar boolean gamma[Periods];
dvar boolean gammar[Periods];

dexpr float originalCost =
  sum(t in Periods)
    ( hc[k0]  * Y[t]
    + hcr[k0] * Yr[t]
    + pc[k0]  * Q[t]
    + pcr[k0] * Qr[t]
    + sc[k0]  * gamma[t]
    + scr[k0] * gammar[t] );

dexpr float reducedCost =
  originalCost
  - sigma
  - sum(t in Periods) piProd[t] * (tp[k0]  * Q[t]  + ts[k0]  * gamma[t])
  - sum(t in Periods) piRem[t]  * (tpr[k0] * Qr[t] + tsr[k0] * gammar[t]);

minimize reducedCost;

subject to {
  // Serviceables balance.
  y0[k0] + Q[1] + Qr[1] == d[k0][1] + Y[1];

  forall(t in 2..T)
    Y[t-1] + Q[t] + Qr[t] == d[k0][t] + Y[t];

  // Recoverables balance.
  yr0[k0] + r[k0][1] == Qr[1] + Yr[1];

  forall(t in 2..T)
    Yr[t-1] + r[k0][t] == Qr[t] + Yr[t];

  // Linking constraints.
  forall(t in Periods) {
    Q[t]  <= M[k0][t] * gamma[t];
    Qr[t] <= M[k0][t] * gammar[t];
  }

  // Valid inequalities (18), separate setups.
  forall(t in 1..T-1, p in 1..T-t)
    Y[t] >=
      sum(s in t+1..t+p) d[k0][s]
      - sum(s in t+1..t+p) M[k0][s] * (gamma[s] + gammar[s]);

  // Valid inequalities (19), separate setups.
  forall(t in 2..T, p in 1..t-1)
    Yr[t] >=
      sum(s in t-p..t) r[k0][s]
      - sum(s in t-p..t) M[k0][s] * gammar[s];
}