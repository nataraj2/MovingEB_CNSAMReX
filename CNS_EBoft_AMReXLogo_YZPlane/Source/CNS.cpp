
#include <CNS.H>
#include <CNS_F.H>

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>

#include <climits>

using namespace amrex;

constexpr int CNS::level_mask_interior;
constexpr int CNS::level_mask_covered;
constexpr int CNS::level_mask_notcovered;
constexpr int CNS::level_mask_physbnd;

constexpr int CNS::NUM_GROW;

BCRec     CNS::phys_bc;

int       CNS::verbose = 0;
IntVect   CNS::hydro_tile_size {AMREX_D_DECL(1024,16,16)};
Real      CNS::cfl       = 0.3;
int       CNS::do_reflux = 1;
int       CNS::refine_cutcells          = 1;
int       CNS::refine_max_dengrad_lev   = -1;
Real      CNS::refine_dengrad           = 1.0e10;
Vector<RealBox> CNS::refine_boxes;

CNS::CNS ()
{}

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if (do_reflux && level > 0) {
        flux_reg.define(bl, papa.boxArray(level-1),
                        dm, papa.DistributionMap(level-1),
                        level_geom, papa.Geom(level-1),
                        papa.refRatio(level-1), level, NUM_STATE);
    }

    buildMetrics();
}

CNS::~CNS ()
{}

void
CNS::init (AmrLevel& old)
{
    auto& oldlev = dynamic_cast<CNS&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);

    MultiFab& C_new = get_new_data(Cost_Type);
    FillPatch(old,C_new,0,cur_time,Cost_Type,0,1);
}

void
CNS::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);

    MultiFab& C_new = get_new_data(Cost_Type);
    FillCoarsePatch(C_new, 0, cur_time, Cost_Type, 0, 1);
}

void
CNS::initData ()
{
    BL_PROFILE("CNS::initData()");

    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      cns_initdata(&level, &cur_time,
                   BL_TO_FORTRAN_BOX(box),
                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                   dx, prob_lo);
    }

    MultiFab& C_new = get_new_data(Cost_Type);
    C_new.setVal(1.0);
}

void
CNS::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& ref_ratio,
                       Vector<Real>&          dt_level,
                       Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  if (level > 0) {
    return;
  }

  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
CNS::computeNewDt (int                    finest_level,
                   int                    sub_cycle,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& ref_ratio,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) {
        return;
    }

    for (int i = 0; i <= finest_level; i++)
    {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else 
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::post_regrid (int lbase, int new_finest)
{
}

void
CNS::post_timestep (int iteration)
{


    if (do_reflux && level < parent->finestLevel()) {
        CNS& fine_level = getLevel(level+1);
        MultiFab& S_crse = get_new_data(State_Type);
        MultiFab& S_fine = fine_level.get_new_data(State_Type);
        fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac);
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }
}

void
CNS::postCoarseTimeStep (Real time)
{
    // This only computes sum on level 0
    if (verbose >= 2) {
        printTotal();
    }
}

void
CNS::printTotal () const
{
    const MultiFab& S_new = get_new_data(State_Type);
    MultiFab mf(grids, dmap, 1, 0);
    std::array<Real,5> tot;
    for (int comp = 0; comp < 5; ++comp) {
        MultiFab::Copy(mf, S_new, comp, 0, 1, 0);
        MultiFab::Multiply(mf, *volfrac, 0, 0, 1, 0);
        tot[comp] = mf.sum(0,true) * geom.ProbSize();
    }
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(tot.data(), 5, ParallelDescriptor::IOProcessorNumber());
            amrex::Print().SetPrecision(17) << "\n[CNS] Total mass       is " << tot[0] << "\n"
                                            <<   "      Total x-momentum is " << tot[1] << "\n"
                                            <<   "      Total y-momentum is " << tot[2] << "\n"
                                            <<   "      Total z-momentum is " << tot[3] << "\n"
                                            <<   "      Total energy     is " << tot[4] << "\n";
#ifdef BL_LAZY
        });
#endif
}

void
CNS::post_init (Real)
{
    if (level > 0) return;
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    if (verbose >= 2) {
        printTotal();
    }
}

void
CNS::post_restart ()
{
}

void
CNS::errorEst (TagBoxArray& tags, int, int, Real time, int, int)
{
    BL_PROFILE("CNS::errorEst()");

    if (refine_cutcells) {
        const MultiFab& S_new = get_new_data(State_Type);
        amrex::TagCutCells(tags, S_new);
    }

    if (!refine_boxes.empty())
    {
        const Real* problo = geom.ProbLo();
        const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
            auto& fab = tags[mfi];
            const Box& bx = fab.box();
            for (BoxIterator bi(bx); bi.ok(); ++bi)
            {
                const IntVect& cell = bi();
                RealVect pos {AMREX_D_DECL((cell[0]+0.5)*dx[0]+problo[0],
                                           (cell[1]+0.5)*dx[1]+problo[1],
                                           (cell[2]+0.5)*dx[2]+problo[2])};
                for (const auto& rbx : refine_boxes) {
                    if (rbx.contains(pos)) {
                        fab(cell) = TagBox::SET;
                    }
                }
            }
        }
    }

    if (level < refine_max_dengrad_lev)
    {
        int ng = 1;
        const auto& rho = derive("density", time, ng);
        const MultiFab& S_new = get_new_data(State_Type);

        const char   tagval = TagBox::SET;
        const char clearval = TagBox::CLEAR;

        auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S_new.Factory());
        auto const& flags = fact.getMultiEBCellFlagFab();
        const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rho,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            const auto& sfab = S_new[mfi];
            const auto& flag = flags[mfi];

            const FabType typ = flag.getType(bx);
            if (typ != FabType::covered)
            {
                cns_tag_denerror(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(tags[mfi]),
                                 BL_TO_FORTRAN_ANYD((*rho)[mfi]),
                                 BL_TO_FORTRAN_ANYD(flag),
                                 &refine_dengrad, &tagval, &clearval, dx, &time);
            }
        }
    }
}

void
CNS::read_params ()
{
    ParmParse pp("cns");

    pp.query("v", verbose);
 
    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
	for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }
   
    pp.query("cfl", cfl);

    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    pp.query("do_reflux", do_reflux);

    pp.query("refine_cutcells", refine_cutcells);

    pp.query("refine_max_dengrad_lev", refine_max_dengrad_lev);
    pp.query("refine_dengrad", refine_dengrad);

    int irefbox = 0;
    Vector<Real> refboxlo, refboxhi;
    while (pp.queryarr(("refine_box_lo_"+std::to_string(irefbox)).c_str(), refboxlo))
    {
        pp.getarr(("refine_box_hi_"+std::to_string(irefbox)).c_str(), refboxhi);
        refine_boxes.emplace_back(refboxlo.data(), refboxhi.data());
        ++irefbox;
    }
}

void
CNS::avgDown ()
{
    BL_PROFILE("CNS::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level+1);

    MultiFab& S_crse =          get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
    volume.setVal(1.0);
    amrex::EB_average_down(S_fine, S_crse, volume, fine_lev.volFrac(),
                           0, S_fine.nComp(), fine_ratio);

    const int nghost = 0;
    computeTemp (S_crse, nghost);
}

void
CNS::buildMetrics ()
{
    BL_PROFILE("CNS::buildMetrics()");

    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    if (std::abs(dx[0]-dx[1]) > 1.e-12*dx[0] || std::abs(dx[0]-dx[2]) > 1.e-12*dx[0]) {
        amrex::Abort("CNS: must have dx == dy == dz\n");
    }

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
    
    volfrac = &(ebfactory.getVolFrac());
    bndrycent = &(ebfactory.getBndryCent());
    centroid = &(ebfactory.getCentroid());
    areafrac = ebfactory.getAreaFrac();
    facecent = ebfactory.getFaceCent();

    level_mask.clear();
    level_mask.define(grids,dmap,1,1);
    level_mask.BuildMask(geom.Domain(), geom.periodicity(), 
                         level_mask_covered,
                         level_mask_notcovered,
                         level_mask_physbnd,
                         level_mask_interior);
}

Real
CNS::estTimeStep ()
{
    BL_PROFILE("CNS::estTimeStep()");

    Real estdt = std::numeric_limits<Real>::max();

    const Real* dx = geom.CellSize();
    const MultiFab& S = get_new_data(State_Type);

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt)
#endif
    {
        Real dt = std::numeric_limits<Real>::max();
        for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();

            const auto& sfab = S[mfi];
            const auto& flag = flags[mfi];

            if (flag.getType(box) != FabType::covered) {
                cns_estdt(BL_TO_FORTRAN_BOX(box),
                          BL_TO_FORTRAN_ANYD(S[mfi]),
                          dx, &dt);
            }
        }
        estdt = std::min(estdt,dt);
    }

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);
    estdt=5.0e-5;
    return estdt;
}

Real
CNS::initialTimeStep ()
{
    return estTimeStep();
}

void
CNS::computeTemp (MultiFab& State, int ng)
{
    BL_PROFILE("CNS::computeTemp()");

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(State.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    // This will reset Eint and compute Temperature 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        const auto& sfab = State[mfi];
        const auto& flag = flags[mfi];

        if (flag.getType(bx) != FabType::covered) {
            cns_compute_temperature(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(State[mfi]));
        }
    }
}

void
CNS::computeVorticity (const MultiFab& State, MultiFab& Vort)
{
    BL_PROFILE("CNS::computeTemp()");

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(State.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
	const Real* dx = geom.CellSize();

    // This will reset Eint and compute Temperature 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(0);

        const auto& sfab = State[mfi];
        const auto& flag = flags[mfi];

        if (flag.getType(bx) != FabType::covered) {
            cns_compute_vorticity(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(State[mfi]), dx, BL_TO_FORTRAN_ANYD(Vort[mfi]));
        }
    }
}



