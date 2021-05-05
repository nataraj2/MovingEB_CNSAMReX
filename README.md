# MovingEB_CNSAMReX
Implementation of moving embedded boundary algorithm for compressible flow in the CNS code in AMReX.
Currently only prescribed motion can be done. i.e. the simulations are not two-way coupled - the 
EB does not move in accordance with the force exerted by the fluid on it. 

## How to run a case 

1. The prescribed motion is specified as a case in one of the `if` loops in `CNS_init_eb2.cpp`.
  ```if(geom_type == "moving_cylinder")

    {
        //EB2::CylinderIF cf1(0.2, 10.0, 2, {4.0+0.0*time,0.0,0.0}, false);
        EB2::CylinderIF cf1(0.2, 10.0, 0, {0.0,0.2*0.4*cos(2.0*3.14159265359*24.375*1.0*time),7.0}, false);
        auto polys = EB2::makeUnion(cf1);
        auto gshop = EB2::makeShop(polys);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4,false);
    }
	```
2. 
3. `CNS_init_eb2.cpp`
