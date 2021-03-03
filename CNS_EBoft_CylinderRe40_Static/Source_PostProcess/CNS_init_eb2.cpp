
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
                const int max_coarsening_level,Real time)
{
    BL_PROFILE("initializeEB2");

    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    if (geom_type == "combustor")
    {
        ParmParse pp("combustor");
        
        Real fwl; 
        pp.get("far_wall_loc",fwl);

        EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
                             {AMREX_D_DECL(1. ,0.,0.)});
        
        Vector<Real> pl1pt, pl2pt, pl2nm, pl3pt; 
        pp.getarr("ramp_plane1_point", pl1pt);
        pp.getarr("ramp_plane2_point", pl2pt);
        pp.getarr("ramp_plane2_normal", pl2nm);
        pp.getarr("ramp_plane3_point", pl3pt);

        auto ramp = EB2::makeIntersection(EB2::PlaneIF({pl1pt[0], pl1pt[1], 0.},
                                                       {      0.,      -1., 0.}),
                                          EB2::PlaneIF({pl2pt[0], pl2pt[1], 0.},
                                                       {pl2nm[0], pl2nm[1], 0.}),
                                          EB2::PlaneIF({pl3pt[0], pl3pt[1], 0.},
                                                       {      1.,       0., 0.}));

        Vector<Real> pipelo, pipehi; 
        pp.getarr("pipe_lo", pipelo);
        pp.getarr("pipe_hi", pipehi);
        
        EB2::BoxIF pipe({pipelo[0], pipelo[1], -1.}, {pipehi[0], pipehi[1], 1.}, false);

        // where does plane 1 and plane 2 intersect?
        Real k2 = std::abs(pl2nm[0]/pl2nm[1]);
        Real secty = pl2pt[1] + k2*(pl3pt[0]-pl2pt[0]);
        // How much do we cut?
        Real dx = geom.CellSize(0);
        Real dycut = 4.*(1.+max_coarsening_level)*std::min(dx, k2*dx);
        EB2::BoxIF flat_corner({pl3pt[0], 0., -1.}, {1.e10, secty+dycut, 1.}, false);
        
        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);

        Real lenx = geom.ProbLength(0);
        Real leny = geom.ProbLength(1);
        auto pr = EB2::translate(EB2::lathe(polys), {lenx*0.5, leny*0.5, 0.});
        
        auto gshop = EB2::makeShop(pr);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }
    else if (geom_type == "moving_plane")
    {
        RealArray point;
        point[0]=0.0+0.5*5e5*time*time;
        point[1]=0.0;
        point[2]=0.0;

        RealArray normal;
        normal[0]=-1.0;
        normal[1]=0.0;
        normal[2]=0.0;

        EB2::PlaneIF pf(point, normal);

        EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        EB2::Build(gshop, geom, max_coarsening_level,
                   max_coarsening_level, 4);

    }
    else if (geom_type == "piston_chamber_compression")
    {
	RealArray point, point2;
        RealArray normal, normal2;

        point[0]=0.0+10.0*time;
        point[1]=0.0;
        point[2]=0.0;

        normal[0]=-1.0;
        normal[1]=0.0;
        normal[2]=0.0;

        point2[0]=0.4;
        point2[1]=0.0;
        point2[2]=0.0;

        normal2[0]=1.0;
        normal2[1]=0.0;
        normal2[2]=0.0;


        EB2::PlaneIF pf(point, normal);
        EB2::PlaneIF pf2(point2, normal2);

        auto polys = EB2::makeUnion(pf, pf2);

	 auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }

    else if (geom_type == "piston_bowl_compression")
    {
	RealArray point;
        RealArray normal;

        RealArray center;
        center[0]=0.5+10*time;
        center[1]=0.0;
        center[2]=0.0;

        Real radius;
        radius=0.5;

        bool has_fluid_inside=1;
        EB2::SphereIF sf(radius, center, has_fluid_inside);

        point[0]=0.4;
        point[1]=0.0;
        point[2]=0.0;

        normal[0]=1.0;
        normal[1]=0.0;
        normal[2]=0.0;


        EB2::PlaneIF pf(point, normal);

        auto polys = EB2::makeUnion(pf, sf);

	 auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }
 
   else if (geom_type == "cylinder_piston_compression")
    {
        RealArray point;
        RealArray normal;

        RealArray center;
        center[0]=0.1+10*time;
        center[1]=0.0;
        center[2]=0.0;

        Real radius;
        radius=0.3;

        bool has_fluid_inside=1;
        EB2::SphereIF sf(radius, center, has_fluid_inside);

        point[0]=0.18;
        point[1]=0.0;
        point[2]=0.0;

        normal[0]=1.0;
        normal[1]=0.0;
        normal[2]=0.0;

        EB2::PlaneIF pf(point, normal);

        EB2::CylinderIF cf1(0.125, 0.8, 0, {-0.2,0.0,0.0}, true);
        //EB2::CylinderIF cf2(0.125, 1.6, 0, {-1.0+10*time,0.0,0.0}, false);

        //auto polys = EB2::makeUnion(cf1, cf2);//, sf);//,pf);//, cf3);
        auto polys = EB2::makeUnion(cf1, sf);//, sf);//,pf);//, cf3);

        auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }

   else if (geom_type == "ICE")
    {
        RealArray point;
        RealArray normal;

        RealArray center;
        center[0]=1.0+9*time;
        center[1]=0.0;
        center[2]=0.0;

        Real radius;
        radius=1.0;

        bool has_fluid_inside=1;
        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::CylinderIF cf1(0.082, 0.09, 0, {0.045,0.0,0.0}, true);

        auto polys = EB2::makeUnion(cf1, sf);

        auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }

   else if (geom_type == "ICE_PistonBowl")
    {
        RealArray point;
        RealArray normal;

        RealArray center;
        center[0]=0.04-0.0125+0.02*cos(3000.0/8.0*time);
        center[1]=0.0;
        center[2]=0.0;

        Real radius;
        radius=0.02;

        bool has_fluid_inside=0;
        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::CylinderIF cf1(0.04, 0.09, 0, {0.045,0.0,0.0}, true);
        EB2::CylinderIF cf2(0.04, 0.125, 0, {-0.0125+0.02*cos(3000.0/8.0*time),0.0,0.0}, false);
        EB2::CylinderIF cf3(0.03, 0.125, 0, {-0.0125+0.02*cos(3000.0/8.0*time),0.0,0.0}, false);
        auto pipe = EB2::makeDifference(cf2,cf3);
        EB2::CylinderIF cf4(0.03, 0.10, 0, {-0.0125+0.02*cos(3000.0/8.0*time),0.0,0.0}, false);

	RealArray center2;
        center2[0]=0.0;
        center2[1]=0.0;
        center2[2]=0.0;

        Real radius2;
        radius2=0.09;

        bool has_fluid_inside2=1;
        EB2::SphereIF sf2(radius2, center2, has_fluid_inside2);

        auto polys = EB2::makeUnion(cf1, pipe, cf4, sf, sf2);

        auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }
   
    else if(geom_type == "moving_cylinder")

    {
        //EB2::CylinderIF cf1(0.2, 10.0, 2, {4.0+0.0*time,0.0,0.0}, false);
        EB2::CylinderIF cf1(0.2, 10.0, 0, {0.0,0.0,7.0+0.0*time}, false);
        auto polys = EB2::makeUnion(cf1);
        auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4,false);
    }


    else if (geom_type == "inclined_duct")
    {
        RealArray point1,point2,point3;
        RealArray normal1,normal2,normal3;

        point1[0]=0.0;
        point1[1]=0.5;
        point1[2]=0.0;

	normal1[0]=0.5;
        normal1[1]=-0.866;
        normal1[2]=0.0;

        point2[0]=0.0;
        point2[1]=0.8;
        point2[2]=0.0;

	normal2[0]=-0.5;
        normal2[1]=0.866;
        normal2[2]=0.0;

	point3[0]=0.4+0.1*cos(30*3.1415/180.0)*time;
        point3[1]=0.5+0.1*sin(30*3.1415/180.0)*time;
        point3[2]=0.0;

        normal3[0]=-0.866;
        normal3[1]=-0.5;
        normal3[2]=0.0;
        

        EB2::PlaneIF pf1(point1, normal1);
        EB2::PlaneIF pf2(point2, normal2);
        EB2::PlaneIF pf3(point3, normal3);

        auto polys = EB2::makeUnion(pf1,pf2,pf3);

	 auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);

        //EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        //EB2::Build(gshop, geom, max_coarsening_level,
         //          max_coarsening_level, 4);

    }
 
    else if (geom_type == "shockpush")
    {
        RealArray point1,point2,point3,point4;
        RealArray normal1,normal2,normal3,normal4;

        point1[0]=0.0;
        point1[1]=1.0;
        point1[2]=0.0;

	normal1[0]=0.0;
        normal1[1]=1.0;
        normal1[2]=0.0;

        point2[0]=0.0;
        point2[1]=0.5;
        point2[2]=0.0;

	normal2[0]=0.0;
        normal2[1]=-1.0;
        normal2[2]=0.0;

	point3[0]=0.5;
        point3[1]=0.0;
        point3[2]=0.0;

        normal3[0]=-0.999999;
        normal3[1]=-0.0014142132;
        normal3[2]=0.0;

	point4[0]=0.1+0.1*time;
        point4[1]=0.0;
        point4[2]=0.0;

        normal4[0]=-1.0;
        normal4[1]=0.0;
        normal4[2]=0.0; 

        EB2::PlaneIF pf1(point1, normal1);
        EB2::PlaneIF pf2(point2, normal2);
        EB2::PlaneIF pf3(point3, normal3);
        EB2::PlaneIF pf4(point4, normal4);

	auto ramp = EB2::makeIntersection(pf2,pf3);
        auto polys = EB2::makeUnion(pf4,pf1,ramp);

	 auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);


        //EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        //EB2::Build(gshop, geom, max_coarsening_level,
         //          max_coarsening_level, 4);

    }
     else if (geom_type == "axishockpush")
    {
	RealArray point, normal;
       	EB2::CylinderIF cf1(0.2, 1.4, 0, {0.5,0.0,0.0}, true);
       	EB2::CylinderIF cf2(0.1, 1.4, 0, {0.5,0.0,0.0}, true);
       	EB2::CylinderIF cf3(0.1,(0.2+337.974*time)*2.0, 0, {0.0+337.974*time,0.0,0.0}, false);

	auto pipe = EB2::makeDifference(cf2,cf1);
	auto polys = EB2::makeUnion(pipe,cf3);

	 auto gshop = EB2::makeShop(polys);

	//EB2::GeometryShop<EB2::CylinderIF> gshop(cf);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }

    else if (geom_type == "moving_sphere")
    {
        RealArray center;
        //pp.get("sphere_center", center);
        center[0]=0.0;
        center[1]=0.0;
        center[2]=0.0;
        std::cout << center[0] << center[1] << center[2] << "\n";

        Real radius;
        //pp.get("sphere_radius", radius);
        radius=0.01+1e-6*cos(2.0*3.1415*1000.0*time);

        bool has_fluid_inside=0;
        //pp.get("sphere_has_fluid_inside", has_fluid_inside);

        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::SphereIF> gshop(sf);
        EB2::Build(gshop, geom, max_coarsening_level,
                   max_coarsening_level, 4);

    }

    else
    {
        EB2::Build(geom, max_coarsening_level, max_coarsening_level, 4);
    }
   }
void finalize_EB2()
   {
     EB2::Finalize();
   }

