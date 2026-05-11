/*********************************************
 * OPL 22.1.2.0 Model
 * Author: Kacper79
 * Creation Date: 10 May 2026 at 10:04:43
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

dvar float+ delta[Products][Cols];

minimize
  sum(k in Products, n in Cols : n <= N[k])
    nc[k][n] * delta[k][n];

subject to {
  forall(k in Products)
    sum(n in Cols : n <= N[k]) delta[k][n] == 1;

  forall(t in Periods)
    sum(k in Products, n in Cols : n <= N[k])
      prodUse[k][n][t] * delta[k][n] <= ct[t];

  forall(t in Periods)
    sum(k in Products, n in Cols : n <= N[k])
      remUse[k][n][t] * delta[k][n] <= ctr[t];

  forall(k in Products, n in Cols : n > N[k])
    delta[k][n] == 0;
}