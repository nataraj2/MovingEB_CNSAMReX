#ifndef CNS_H_
#define CNS_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>

    void GlobalTest();

using namespace amrex;
class CNS
    :
    public amrex::AmrLevel
{
public:

    CNS ();
    CNS (amrex::Amr&            papa,
         int                    lev,
         const amrex::Geometry& level_geom,
         const amrex::BoxArray& bl,
         const amrex::DistributionMapping& dm,
         amrex::Real            time);
    virtual ~CNS ();

    CNS (const CNS& rhs) = delete;
    CNS& operator= (const CNS& rhs) = delete;
    
    // Restart from a checkpoint file.
    virtual void restart (amrex::Amr&     papa,
                          std::istream&   is,
                          bool            bReadSpecial = false) override;

    // Write checkpoint
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            amrex::VisMF::How  how = amrex::VisMF::NFiles,
                            bool               dump_old = true) override;

    virtual std::string thePlotFileType () const override {
        return {"HyperCLaw-V1.1"};
    }

    // Write a plotfile to specified directory.
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;

    // Initialize data on this level from another CNS (during regrid).
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;
    
    // Initialize grid data at problem start-up.
    virtual void initData () override;

    // Advance grids at this level in time.
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;

    virtual void computeInitialDt (int                                 finest_level,
                                   int                                 sub_cycle,
                                   amrex::Vector<int>&                  n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                         stop_time) override;

    virtual void computeNewDt (int                                 finest_level,
                               int                                 sub_cycle,
                               amrex::Vector<int>&                  n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                         stop_time,
                               int                                 post_regrid_flag) override;

    virtual void post_regrid (int lbase, int new_finest) override;

    // Do work after timestep().
    virtual void post_timestep (int iteration) override;

    // After a full time step
    virtual void postCoarseTimeStep (amrex::Real time) override;

    // Do work after init().
    virtual void post_init (amrex::Real stop_time) override;

    virtual void post_restart () override;

    // Error estimation for regridding.
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
                           int                 n_error_buf = 0,
                           int                 ngrow = 0) override;

    virtual int WorkEstType () override { return Cost_Type; }

    // Define data descriptors.
    static void variableSetUp ();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp ();

    static int numGrow() { return NUM_GROW; };
   
    void ZeroOutSolidWalls (const int* lev, const amrex::MultiFab& S);
    void ComputeUncoveringRedistribute (const amrex::MultiFab& volfracborder, amrex::MultiFab& Sborder, amrex::MultiFab& cons_adjust, amrex::MultiFab& donortag);
    void FinalUncoveringRedistribute (const amrex::MultiFab& volfracborder, amrex::MultiFab& Sborder, amrex::MultiFab& cons_adjust, amrex::MultiFab& donortag);

    //MultiFab volfracmfab;

protected:

    static void read_params ();

    const amrex::MultiFab& volFrac () const { return *volfrac; }

    CNS& getLevel (int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

    void avgDown ();

    void buildMetrics ();

    amrex::Real estTimeStep ();

    // Compute initial time step.
    amrex::Real initialTimeStep ();

    void computeTemp (amrex::MultiFab& State, int ng);
	amrex::MultiFab Vort;
    void computeVorticity (const amrex::MultiFab& State, amrex::MultiFab& Vort);

    void compute_dSdt (const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
                       amrex::EBFluxRegister* fr_as_crse, amrex::EBFluxRegister* fr_as_fine, Real Time, int RKstep);

    void printTotal () const;

    const amrex::MultiFab* volfrac;
    const amrex::MultiCutFab* bndrycent;
    const amrex::MultiCutFab* centroid;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> facecent;

    amrex::iMultiFab level_mask;

    amrex::EBFluxRegister flux_reg;

    static constexpr int level_mask_interior   = 0; // valid cells
    static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2; // ghost cells not covered
    static constexpr int level_mask_physbnd    = 3; // outside domain

    static constexpr int NUM_GROW = 5;
    enum StateVariable {
        Density = 0, Xmom, Ymom, Zmom, Eden, Eint, Temp, volfraction,
        NUM_STATE
    };

    enum StateDataType {
        State_Type = 0,
        Cost_Type
    };
    static int num_state_data_types;

    static amrex::BCRec phys_bc;

    // Parameters
    static int verbose;
    static amrex::IntVect hydro_tile_size;
    static amrex::Real cfl;

    static int do_reflux;

    static int refine_cutcells;

    static int refine_max_dengrad_lev;
    static amrex::Real refine_dengrad;

    static amrex::Vector<amrex::RealBox> refine_boxes;
};

#endif
