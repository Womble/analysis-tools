#define MAXPE 200
#define MAXLEV 19
#define MAXQ 9
#define MAXNAL 2
#define MAXD 3
#ifdef DOUBLE
#define float double
#endif
struct join {
  int ijoin, bound, done;
  bool empty;
  float *pa0, *dpdt, *dpdx, *cor, xj[MAXD], area;
  struct cell *left, *right;
  struct join  *prev, *nxt, *par;
};
struct cell {
  int icell, ichild, numtimesteps;
  bool refined, empty, halo, gotdiff;
  float qsol, dqsoldt, qa[MAXQ], pa[MAXQ], dqdt[MAXQ], diff[MAXQ],
    xc[MAXD], vol, dudx[MAXD][MAXD], radR,radz;

  struct cell *child[1<<MAXD], *par, *nxt, *prev;
  struct join *ljoin[MAXD], *rjoin[MAXD];
  struct mcell *mcptr;
};

struct level {
  int nc[MAXD], nt, ndt;
  float dxbar[MAXD], xl[MAXD], xr[MAXD];
  struct cell *fcell, *lcell;
  struct join **fjoin, **ljoin;
};

struct cell_list {
  struct cell *cptr;
  struct cell_list *prev, *nxt;
};

struct join_list {
  struct join *jptr;
  struct join_list *nxt;
};
struct mcell {
  float qa, qt, mso, pso, diff, vol, xc[MAXD];
  struct mcell *child[1<<MAXD], *lcell[MAXD], *rcell[MAXD], *par, *nxt;
  struct cell *cptr;
};
