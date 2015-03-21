


/*
 * simulation.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_


#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "sources.h"
#include "materials.h"
#include "initialization_values.h"


using namespace std;





/*------------ simulation class ------------------ */


class Simulation{

public:

    //!Constructor. Initializes the fields; initializes all update coefficients to free space values.
    /*! number of magnetic nodes = number of electric nodes - 1;
     * in this way the grid ends at both sides with an electric node (so e.g., the boundary conditions can be PEC at both sides)
     *
     * The fields are initialized (i.e., set to certain values before the simulation starts) by making use of an initialization
     * function, passed as a pointer. "initialization_value.h" contains the predefined (pointers to the) functions standing_wave and all_zero.
     *
     * \param name: name of the simulation, used for generating output filenames
     * \param dimension: total number of spatial nodes
     * \param S_c : courant number, S_c = delta_t * c / delta_x; mostly set to 1
     * \param delta_t :  discrete unit of time; delta_x follows from choice of S_c and delta_t
     * \param type_of_BDC: type of boundary condition, apart from PMLs. "ABC" = absorbing boundary condition, "PEC" = perfect electric conductor
     *  */
    Simulation(string name_, int dimension_, double S_c_, double delta_t_, string type_of_BDC_, double (*init_func_)(int dim_grid, int x));

    //! adds Source objects and a vector of Static_Material objects (3D materials); sets update coefficients at material nodes to correct values
    /*!
     * \param source_e: Source object containing the information of the (analytical) source signal for the E field
     * \param source_h: Source object containing the information of the (analytical) source signal for the H field
     * \param static_layers: vector containing Static_Material objects, i.e., the 3D non time-dependent materials such as PML's
     */
    void static_setup(Source& source_e, Source& source_h, vector <Static_Material>& static_layers);

    //! adds a vector of Dispersive_Metasurface Objects to the simulation, containing the information about the metasurfaces with time-dependent conductivity
    /*!
     * \param metasurfaces: vector containing Dispersive_Metasurface objects, i.e. the 2D materials with time-dependent conductivity
     */
    void add_metasurfaces(vector <Dispersive_Metasurface>& metasurfaces);

    //! Performs the simulation, with different output options for the generated data.
    /*!
     * - Data can be saved in the form of snapshots (i.e., saving the values of all spatial nodes at certain timesteps)
     * - and/or the total signal in the time domain can be saved at certain nodes.
     *
     * This can be done for the electric and magnetic field values independently.
     *
     * For generating snapshots, a snapshotmodulo differing from -1 should be specified. A snapshotmodulo of e.g. 10
     * will generate 1 snapshot in 10 timesteps. The generated snapshots can be used to create so-called "waterfall plots".
     *
     * To examine the signal in the time domain for certain nodes, a non-empty vector of nodes (integers) should be specified.
     * This results in a single file per node and per field type.
     *
     * \param total_timesteps : the total number of timesteps simulated. Timesteps * delta_t = real time simulated.
     * \param output_path : directory where output files will be located
     * \param snapshotmodulo_electric : e.g., snapshotmodulo = 10 will generate 1 snapshot every 10 timesteps. Specify -1 if no electric snapshot files desired (speeds up simulation)
     * \param snapshotmodulo_magnetic : e.g., snapshotmodulo = 10 will generate 1 snapshot every 10 timesteps. Specify -1 if no magnetic snapshot files desired (speeds up simulation)
     * \param readout_nodes_electric : vector containing the number of the nodes for which the time domain signal is saved. Specify an empty vector if saving of time domain signal is not desired (speeds up simulation)
     * \param readout_nodes_magnetic : vector containing the number of the nodes for which the time domain signal is saved. Specify an empty vector if saving of time domain signal is not desired (speeds up simulation)
     */
    void simulate(int total_timesteps,string output_path, int snapshotmodulo_electric, int snapshotmodulo_magnetic, vector<int> readout_nodes_electric, vector<int> readout_nodes_magnetic );


private:
    string name;
    int dimension;
    double S_c;
    double delta_t;
    string type_of_BDC;

    Source source_e;
    Source source_h;

	vector<double> e; //vector containing the E field values at all spatial nodes at a certain timestep, updated at every step in the simulation
	vector<double> h; //vector containing the H field values at all spatial nodes at a certain timestep, updated at every step in the simulation

	vector <Static_Material> static_layers;
	vector <Dispersive_Metasurface> metasurfaces;

 	vector <double> ce_e, ce_h; //electric-field update coefficients for the e- and h-terms of the update equation
	vector <double> ch_e, ch_h; //magnetic-field update coefficients for the e- and h-terms of the update equation


	//methods for performing the simulation

	/*double calculate_convolution_term_e(int q, Dispersive_Metasurface& sheet);
	double calculate_convolution_term_h(int q, Dispersive_Metasurface& sheet);*/

	void apply_source(double& timestep, vector <double>& field, Source& source);

    void update_magnetic_once(double& timestep);
    void update_electric_once(double& timestep);
    void update_system_once(double& timestep);

    ofstream* file_Q_electric;
    ofstream* file_Q_magnetic;

    // methods to write the output to files

	string generate_filename(int number, string field);

	ofstream* generate_file(string basename);
	vector<ofstream*> generate_node_output_files(string output_path, vector<int>& readout_nodes, string field);
	ofstream* generate_snapshot_file(string output_path, int framenumber, string field);

	void value_to_file(double& value, ofstream* file);
	void data_to_snapshot_file(ofstream* file, string field);

};




#endif /* SIMULATION_H_ */
