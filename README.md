# ww3-gridgen-mo
Python and fortran based grid generation code for WAVEWATCH III models

V0 repo for grid generation, WAVEWATCH III formatting, display and propagation tests for Spherical Multiple-Cell grids
Contents as follows:

  Directory gridgen - core python3 code to generate regularGrid and smcGrid python class data from netCDF bathymetry, save staging data as
                      netCDF files and write out WW3 compatible text files. Accompanying code to pre-process GEBCO data (reduces memory
                      requirements for larger grids) and generate cell face array linkages for SMC grids (bash and fortran).
                      
  Directory gridedit - python3 code to create google-earth files from grid data, remove cells from SMC grids based on kml files, check
                       and update boundary cells for edited SMC grids.
                       
  Directory smc_test - bash scripts, fortran and python3 code to run SMC grid propagation tests and visualize grids.
  
  Directory build_mo_global - example python3 code to run gridgen for build of Met Office global model configurations at GW1.0. These
                              scripts build grids from coarsest to finest cells, saving the parent grid and grid components in netCDF
                              format as they are generated.
                              
  Directory build_mo_atlantic - example python3 code to run gridgen for build of Met Office Atlantic model configuration at GW1.0. These
                                scripts build grids from coarsest to finest cells, saving the parent grid and grid components in netCDF
                                format as they are generated.
                              
