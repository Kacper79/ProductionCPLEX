/*********************************************
 * OPL 22.1.2.0 Model
 * Author: Kacper79
 * Creation Date: 7 May 2026 at 12:01:25
 *********************************************/

int K = ...;
int T = ...;
range Products = 1..K;
range Periods = 1..T;

float dummyCost = ...;
int maxIterations = ...;
float rcTolerance = ...;

float ct[Periods] = ...;
float ctr[Periods] = ...;

float d[Products][Periods] = ...;
float r[Products][Periods] = ...;

float hc[Products] = ...;
float hcr[Products] = ...;
float pc[Products] = ...;
float pcr[Products] = ...;
float sc[Products] = ...;
float scr[Products] = ...;

float tp[Products] = ...;
float tpr[Products] = ...;
float ts[Products] = ...;
float tsr[Products] = ...;

float y0[Products] = ...;
float yr0[Products] = ...;

float M[k in Products][t in Periods] =
   maxl(0, maxl(ct[t] - ts[k], ctr[t] - tsr[k])) / minl(tp[k], tpr[k]);

execute {
  writeln("Loaded CLSP-RM-SS instance: K=", K, ", T=", T);
}

main {
  thisOplModel.generate();

  var K = thisOplModel.K;
  var T = thisOplModel.T;

  var dummyCost = thisOplModel.dummyCost;
  var maxIterations = thisOplModel.maxIterations;
  var rcTolerance = thisOplModel.rcTolerance;

  var ct  = thisOplModel.ct;
  var ctr = thisOplModel.ctr;

  var d = thisOplModel.d;
  var r = thisOplModel.r;

  var hc  = thisOplModel.hc;
  var hcr = thisOplModel.hcr;
  var pc  = thisOplModel.pc;
  var pcr = thisOplModel.pcr;
  var sc  = thisOplModel.sc;
  var scr = thisOplModel.scr;

  var tp  = thisOplModel.tp;
  var tpr = thisOplModel.tpr;
  var ts  = thisOplModel.ts;
  var tsr = thisOplModel.tsr;

  var y0  = thisOplModel.y0;
  var yr0 = thisOplModel.yr0;
  var M   = thisOplModel.M;

  var dualSource    = new IloOplModelSource("projekt-artykul-SPP-Dual.mod");
  var primalSource  = new IloOplModelSource("projekt-artykul-SPP-Primal.mod");
  var pricingSource = new IloOplModelSource("projekt-artykul-pricing.mod");

  var dualDef    = new IloOplModelDefinition(dualSource);
  var primalDef  = new IloOplModelDefinition(primalSource);
  var pricingDef = new IloOplModelDefinition(pricingSource);

  var N = new Array(K + 1);
  var cost = new Array(K + 1);
  var prodUse = new Array(K + 1);
  var remUse = new Array(K + 1);
  var schedQ = new Array(K + 1);
  var schedQr = new Array(K + 1);
  var schedY = new Array(K + 1);
  var schedYr = new Array(K + 1);
  var schedGamma = new Array(K + 1);
  var schedGammar = new Array(K + 1);

  for (var k = 1; k <= K; k++) {
    N[k] = 1;
    cost[k] = new Array();
    prodUse[k] = new Array();
    remUse[k] = new Array();
    schedQ[k] = new Array();
    schedQr[k] = new Array();
    schedY[k] = new Array();
    schedYr[k] = new Array();
    schedGamma[k] = new Array();
    schedGammar[k] = new Array();

    cost[k][1] = dummyCost;
    prodUse[k][1] = new Array(T + 1);
    remUse[k][1] = new Array(T + 1);
    schedQ[k][1] = new Array(T + 1);
    schedQr[k][1] = new Array(T + 1);
    schedY[k][1] = new Array(T + 1);
    schedYr[k][1] = new Array(T + 1);
    schedGamma[k][1] = new Array(T + 1);
    schedGammar[k][1] = new Array(T + 1);

    for (var t = 1; t <= T; t++) {
      prodUse[k][1][t] = 0;
      remUse[k][1][t] = 0;
      schedQ[k][1][t] = 0;
      schedQr[k][1][t] = 0;
      schedY[k][1][t] = 0;
      schedYr[k][1][t] = 0;
      schedGamma[k][1][t] = 0;
      schedGammar[k][1][t] = 0;
    }
  }

  function writeVector(f, name, a, from, to) {
    f.write(name, " = [");
    for (var i = from; i <= to; i++) {
      if (i > from) f.write(" ");
      f.write(a[i]);
    }
    f.writeln("]; ");
  }

  function writeMatrix2D(f, name, a, dim1, dim2) {
    f.writeln(name, " = [");
    for (var i = 1; i <= dim1; i++) {
      f.write("  [");
      for (var j = 1; j <= dim2; j++) {
        if (j > 1) f.write(" ");
        f.write(a[i][j]);
      }
      f.writeln("]");
    }
    f.writeln("]; ");
  }

  function writeMatrix3D(f, name, a, dim1, N, dim3, maxN) {
    f.writeln(name, " = [");
    for (var k = 1; k <= dim1; k++) {
      f.writeln("  [");
      for (var n = 1; n <= maxN; n++) {
        f.write("    [");
        for (var t = 1; t <= dim3; t++) {
          if (t > 1) f.write(" ");
          if (n <= N[k]) f.write(a[k][n][t]);
          else f.write(0);
        }
        f.writeln("]");
      }
      f.writeln("  ]");
    }
    f.writeln("]; ");
  }

  function writeMasterDataFile(fileName, maxN) {
    var f = new IloOplOutputFile(fileName);
    f.writeln("K = ", K, ";");
    f.writeln("T = ", T, ";");
    f.writeln("maxN = ", maxN, ";");
    f.write("N = [");
    for (var k = 1; k <= K; k++) {
      if (k > 1) f.write(" ");
      f.write(N[k]);
    }
    f.writeln("]; ");

    writeVector(f, "ct", ct, 1, T);
    writeVector(f, "ctr", ctr, 1, T);

    f.writeln("nc = [");
    for (var k = 1; k <= K; k++) {
      f.write("  [");
      for (var n = 1; n <= maxN; n++) {
        if (n > 1) f.write(" ");
        if (n <= N[k]) f.write(cost[k][n]);
        else f.write(0);
      }
      f.writeln("]");
    }
    f.writeln("]; ");

    writeMatrix3D(f, "prodUse", prodUse, K, N, T, maxN);
    writeMatrix3D(f, "remUse",  remUse,  K, N, T, maxN);
    f.close();
  }

  function writePricingDataFile(fileName, k0, sigma, piProd, piRem) {
    var f = new IloOplOutputFile(fileName);
    f.writeln("K = ", K, ";");
    f.writeln("T = ", T, ";");
    f.writeln("k0 = ", k0, ";");

    writeMatrix2D(f, "d", d, K, T);
    writeMatrix2D(f, "r", r, K, T);
    writeVector(f, "hc", hc, 1, K);
    writeVector(f, "hcr", hcr, 1, K);
    writeVector(f, "pc", pc, 1, K);
    writeVector(f, "pcr", pcr, 1, K);
    writeVector(f, "sc", sc, 1, K);
    writeVector(f, "scr", scr, 1, K);
    writeVector(f, "tp", tp, 1, K);
    writeVector(f, "tpr", tpr, 1, K);
    writeVector(f, "ts", ts, 1, K);
    writeVector(f, "tsr", tsr, 1, K);
    writeVector(f, "y0", y0, 1, K);
    writeVector(f, "yr0", yr0, 1, K);
    writeMatrix2D(f, "M", M, K, T);
    f.writeln("sigma = ", sigma, ";");
    writeVector(f, "piProd", piProd, 1, T);
    writeVector(f, "piRem", piRem, 1, T);
    f.close();
  }

  function scheduleAlreadyExists(k, q, qr, y, yr, g, gr) {
    var eps = 1e-6;
    for (var n = 1; n <= N[k]; n++) {
      var same = true;
      for (var t = 1; t <= T; t++) {
        if (Math.abs(schedQ[k][n][t] - q[t]) > eps) same = false;
        if (Math.abs(schedQr[k][n][t] - qr[t]) > eps) same = false;
        if (Math.abs(schedY[k][n][t] - y[t]) > eps) same = false;
        if (Math.abs(schedYr[k][n][t] - yr[t]) > eps) same = false;
        if (Math.abs(schedGamma[k][n][t] - g[t]) > eps) same = false;
        if (Math.abs(schedGammar[k][n][t] - gr[t]) > eps) same = false;
      }
      if (same) return true;
    }
    return false;
  }

  var masterObj = 0;

  for (var iter = 1; iter <= maxIterations; iter++) {
    var maxN = 0;
    for (var kk = 1; kk <= K; kk++) if (N[kk] > maxN) maxN = N[kk];

    writeMasterDataFile("tmp_master.dat", maxN);

    var cplexDual = new IloCplex();
    var dualModel = new IloOplModel(dualDef, cplexDual);
    var dualData = new IloOplDataSource("tmp_master.dat");
    dualModel.addDataSource(dualData);
    dualModel.generate();

    if (!cplexDual.solve()) {
      writeln("Dual RMP was not solved.");
      break;
    }

    masterObj = cplexDual.getObjValue();

    var sigma = new Array(K + 1);
    var piProd = new Array(T + 1);
    var piRem = new Array(T + 1);

    for (var k = 1; k <= K; k++) {
      sigma[k] = dualModel.sigma[k].solutionValue;
    }
    for (var t = 1; t <= T; t++) {
//      piProd[t] = cplexDual.getValue(dualModel.piProd[t]);
//      piRem[t]  = cplexDual.getValue(dualModel.piRem[t]);
      piProd[t] = dualModel.piProd[t].solutionValue;
	  piRem[t]  = dualModel.piRem[t].solutionValue;
    }

    writeln();
    writeln("Iteration ", iter);
    writeln("  RMP dual objective = ", masterObj);
    write("  Columns per product before pricing: ");
    for (var k = 1; k <= K; k++) write(N[k], " ");
    writeln();

    dualModel.end();
    dualData.end();
    cplexDual.end();

    var addedAny = false;

    for (var k = 1; k <= K; k++) {
      writePricingDataFile("tmp_pricing.dat", k, sigma[k], piProd, piRem);

      var cplexPricing = new IloCplex();
      var pricing = new IloOplModel(pricingDef, cplexPricing);
      var pricingData = new IloOplDataSource("tmp_pricing.dat");
      pricing.addDataSource(pricingData);
      pricing.generate();

      if (!cplexPricing.solve()) {
        writeln("  Pricing product ", k, " not solved.");
        pricing.end();
        pricingData.end();
        cplexPricing.end();
        continue;
      }

      var rc = cplexPricing.getObjValue();
      writeln("  Product ", k, " pricing rc = ", rc);

      if (rc < -rcTolerance) {
        var q = new Array(T + 1);
        var qr = new Array(T + 1);
        var y = new Array(T + 1);
        var yr = new Array(T + 1);
        var g = new Array(T + 1);
        var gr = new Array(T + 1);

        var originalCost = 0;
        var prodCap = new Array(T + 1);
        var remCap = new Array(T + 1);

        for (var t = 1; t <= T; t++) {
          q[t]  = pricing.Q[t].solutionValue;
		  qr[t] = pricing.Qr[t].solutionValue;
		  y[t]  = pricing.Y[t].solutionValue;
		  yr[t] = pricing.Yr[t].solutionValue;
		  g[t]  = pricing.gamma[t].solutionValue;
		  gr[t] = pricing.gammar[t].solutionValue;

          originalCost += hc[k] * y[t] + hcr[k] * yr[t]
                        + pc[k] * q[t] + pcr[k] * qr[t]
                        + sc[k] * g[t] + scr[k] * gr[t];
		  

          prodCap[t] = tp[k] * q[t] + ts[k] * g[t];
          remCap[t]  = tpr[k] * qr[t] + tsr[k] * gr[t];
        }

        if (!scheduleAlreadyExists(k, q, qr, y, yr, g, gr)) {
          N[k] = N[k] + 1;
          var nnew = N[k];
          cost[k][nnew] = originalCost;
          prodUse[k][nnew] = new Array(T + 1);
          remUse[k][nnew] = new Array(T + 1);
          schedQ[k][nnew] = new Array(T + 1);
          schedQr[k][nnew] = new Array(T + 1);
          schedY[k][nnew] = new Array(T + 1);
          schedYr[k][nnew] = new Array(T + 1);
          schedGamma[k][nnew] = new Array(T + 1);
          schedGammar[k][nnew] = new Array(T + 1);

          for (var t = 1; t <= T; t++) {
            prodUse[k][nnew][t] = prodCap[t];
            remUse[k][nnew][t]  = remCap[t];
            schedQ[k][nnew][t] = q[t];
            schedQr[k][nnew][t] = qr[t];
            schedY[k][nnew][t] = y[t];
            schedYr[k][nnew][t] = yr[t];
            schedGamma[k][nnew][t] = g[t];
            schedGammar[k][nnew][t] = gr[t];
          }

          addedAny = true;
          writeln("    -> added column #", nnew, " original cost=", originalCost);
        } else {
          writeln("    -> duplicate ignored");
        }
      }

      pricing.end();
      pricingData.end();
      cplexPricing.end();
    }

    if (!addedAny) {
      writeln("No improving columns. CG stops.");
      break;
    }
  }

  writeln();
  writeln("Final CG lower bound from dual RMP = ", masterObj);
  writeln("Final columns per product:");
  for (var k = 1; k <= K; k++){
    writeln("  Product ", k, ": ", N[k], " columns");
  }
   

  dualDef.end();
  primalDef.end();
  pricingDef.end();
  dualSource.end();
  primalSource.end();
  pricingSource.end();
}  