MEX = /Applications/MATLAB_R2014a.app/bin/mex
BOOSTFLAG = -I/usr/local/include
CXXFLAGS = -Ofast -std=c++11

CXX = clang++


default:
					# compile ecoli_bacterium.cpp and concentration_grid.cpp as one file
					# mex will fail to link otherwise
					cat ./ecoli_bacterium.cpp ./concentration_grid.cpp >> ./bac_conc_combined.cpp
					$(CXX) $(BOOSTFLAG) $(CXXFLAGS) -O -c ./bac_conc_combined.cpp
					rm ./bac_conc_combined.cpp
					$(CXX) $(BOOSTFLAG) $(CXXFLAGS) -O -c ./geometry.cpp
					$(CXX) $(BOOSTFLAG) $(CXXFLAGS) -O -c ./random.cpp
					$(CXX) $(BOOSTFLAG) $(CXXFLAGS) -O -c ./simulation_wrapper.cpp

					$(MEX) $(BOOSTFLAG) ./mx_alloc_simulation.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_dealloc_simulation.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_get_cell_positions.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_setup_border.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_setup_producers.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_update_simulation.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_insert_cells.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_get_concentration.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_set_diffusion_params.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_set_growth_params.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o
					$(MEX) $(BOOSTFLAG) ./mx_get_producer_conc.cpp ./bac_conc_combined.o ./geometry.o ./random.o ./simulation_wrapper.o

clean:
					rm *.o
					rm *.mexmaci64
