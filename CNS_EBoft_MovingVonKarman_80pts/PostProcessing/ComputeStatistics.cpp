/*
  A very simple example of reading a plotfile and calling a function to perform a pointwise transformation
  based on a set of components that are specified by name on the command line.  The transformation is done
  in the accompanying fortran routine.  No grow cells are used, so the transformation cannot involve a 
  stencil operation.

  The output is a new plotfile with a single component set to the output of the transform routine.  This 
  new plotfile has metadata (number of levels, boxarray, grid spacing, etc) that is identical to the original
  plotfile.
 */
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>

extern "C" {
  void transform (const int* lo, const int* hi,
                  const amrex_real* sIn, const int* sInlo, const int* sInhi, const int* ncIn, double dx[3], double problo[3], double volume_and_mass[3]);
}

using namespace amrex;

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}


struct AMReXMeshHierarchy
{
/*
  A minimal class to describe the AMR hierarchy for analysis routines
 */
public:
  AMReXMeshHierarchy() {}
  void define(const AmrData& ad) {
    finestLevel = ad.FinestLevel();
    int nlevs = finestLevel + 1;
    ba.resize(nlevs);
    probSize = ad.ProbSize();
    probDomain = ad.ProbDomain();
    refRatio = ad.RefRatio();
    for (int lev=0; lev<nlevs; ++lev) {
      ba[lev] = &ad.boxArray(lev);
    }
  }
  int FinestLevel() const {return finestLevel;}
  const BoxArray &boxArray(int level) const {return *ba[level];}
  const Vector<int> &RefRatio() const {return refRatio;}
  const Vector<Real> &ProbSize() const {return probSize;}
  const Vector<Box> &ProbDomain() const {return probDomain;}

protected:
  int finestLevel;
  std::vector<const BoxArray*> ba;
  Vector<int> refRatio;
  Vector<Real> probSize;
  Vector<Box> probDomain;
};


struct AMReXDataHierarchy
{
/*
  Data on a AMReXMeshHierarchy, currently pointing to MultiFabs of
  named variables managed by an AmrData object.
*/
public:
  AMReXDataHierarchy(AmrData& ad, const Vector<std::string>& varNames) {
    mesh.define(ad);
    const Vector<std::string>& plotVarNames = ad.PlotVarNames();
    int nComp = varNames.size();
    int nlevs = mesh.FinestLevel() + 1;
    for (int i=0; i<nComp; ++i) {
      int idx = -1;
      for (int j=0; j<plotVarNames.size() && idx<0; ++j) {
        if (plotVarNames[j] == varNames[i]) {idx = j;}
      }
      if (ParallelDescriptor::IOProcessor() && idx<0) {
        Abort("Cannot find variable="+varNames[i]+" in pltfile");
      }
      std::vector<MultiFab*> mfs(nlevs);
      for (int lev=0; lev<nlevs; ++lev) {
        mfs[lev] = &ad.GetGrids(lev,idx); // Note: This lazily triggers a MultiFab read in the AmrData
      }
      varMap[varNames[i]] = std::make_pair(0,mfs);
    }
  }

  MultiFab &GetGrids(int level, const std::string& name) {
    if (varMap.find(name) == varMap.end()) {
      Abort("Unknown component requested");
    }
    return *(varMap[name].second[level]);
  }

  const AMReXMeshHierarchy& Mesh() const {return mesh;}

protected:
  AMReXMeshHierarchy mesh;
  std::map<std::string,std::pair<int,std::vector<MultiFab*>>> varMap;
};


int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    std::string infile;
    std::string outfile;
    int counter = 0;
    int isize = 32;
    int jsize = 32;
    int ksize = 32;

    double dx[3];
    double problo[3];
    double xvelmean[ksize], xvelmean1[ksize], xvelmean2[ksize], time_1, time_2, time_start, time_end;
    double volume_and_mass[3];
    using std::ofstream;
    ofstream outdata;

    outdata.open("DensityVsTime.txt"); // opens the file

    for (int filenum=0; filenum < 17000; filenum=filenum+100){
	amrex::Print() << "filenum = "  << filenum << "\n";
    	if(filenum<10){
	    counter = counter+1;
	    infile = "plt0000" + std::to_string(filenum);
    	    outfile = "pressure0000" + std::to_string(filenum) + ".txt";
	}
	else if(filenum<100){
	    counter = counter+1;
	    infile = "plt000" + std::to_string(filenum);
    	    outfile = "pressure000" + std::to_string(filenum) + ".txt";
	}
	else if(filenum<1000){
	    counter = counter+1;
	    infile = "plt00" + std::to_string(filenum);
    	    outfile = "pressure00" + std::to_string(filenum) + ".txt";
	}
	else if(filenum<10000){
	    counter = counter+1;
	    infile = "plt0" + std::to_string(filenum);
    	    outfile = "pressure0" + std::to_string(filenum) + ".txt";
	}
	else if(filenum<100000){
	    counter = counter+1;
	    infile = "plt" + std::to_string(filenum);
    	    outfile = "pressure" + std::to_string(filenum) + ".txt";
	}

    //outdata.open(outfile); // opens the file
    amrex::Print() << infile << "\n";
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();
    PlotFileData plotfile(infile);

	 int nv = 3;
        Vector<std::string> varNames(nv);
        varNames[0] = "density";
        varNames[1] = "vfrac";
        varNames[2] = "pressure";

 // Make a data struct for just the variables needed
    AMReXDataHierarchy data(amrData,varNames);
    const AMReXMeshHierarchy& mesh = data.Mesh();

    	const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    	int nCompIn = varNames.size();
    	Vector<int> destFillComps(nCompIn);
        for (int i=0; i<nCompIn; ++i) {
     	    destFillComps[i] = i;
	    }

        for (int i=0; i<nCompIn; ++i) {
	    int ivar = -1;
      	     for (int j=0; j<plotVarNames.size(); ++j) {
	        if (plotVarNames[j] == varNames[i]) {ivar = j;}
	      }
	      if (ParallelDescriptor::IOProcessor() && ivar<0) {
        	Abort("Cannot find variable="+varNames[i]+" in pltfile");
	      }
        }

    	const int nCompOut = 1;
    	const int nGrow = 0;
    	const int nLev = amrData.FinestLevel() + 1;

    	Array<Real,AMREX_SPACEDIM> prob_lo = plotfile.probLo();
    	Array<Real,AMREX_SPACEDIM> dx0 = plotfile.cellSize(0);
   

    	for (int i=0; i<3; i++){
		dx[i] = dx0[i];
		problo[i] = prob_lo[i];
    	}

        
    
    const int finestLevel = mesh.FinestLevel();
    Vector<MultiFab*> stateOut(nLev);

    double mass = 0.0;
    double volume = 0.0;
    double pressure = 0.0;
    for (int lev=0; lev<nLev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);

      const DistributionMapping dmap(ba);
      MultiFab stateIn(ba,dmap,nCompIn,nGrow);

      // Load input data from pltfile 
      amrData.FillVar(stateIn,lev,varNames,destFillComps);

      double lev_double = lev*1.0;
      for (int i=0; i<3; i++){
	  dx[i]=dx0[i]/pow(2.0,lev_double); 
	}

      // Compute transformation
      for (MFIter mfi(stateIn); mfi.isValid(); ++mfi) {

        const FArrayBox& sIn = stateIn[mfi];
        const Box& box = mfi.validbox();
	
        transform(box.loVect(), box.hiVect(),
                  BL_TO_FORTRAN_ANYD(sIn),&nCompIn, dx, problo, volume_and_mass);

	volume = volume + volume_and_mass[0];
	mass = mass + volume_and_mass[1];
	pressure = pressure + volume_and_mass[2];

      }	
      }

	amrex::Print() << "Density = " << mass/volume << "\n";
	outdata << plotfile.time() << " "  << volume << " " << mass/volume << " " << pressure/volume << "\n";
	
     }
	     outdata.close();
   } 
  amrex::Finalize();
  return 0;
}
