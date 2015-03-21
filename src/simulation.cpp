/*
 * simulation.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#include "simulation.h"

#include <time.h>


double lower_value_limit = pow(10.0,-50); //smaller values are taken to be zero


//ÑÑÑÑÑ simulation constructor

Simulation::Simulation(string name_, int dimension_, double S_c_, double delta_t_, string type_of_BDC_, double (*init_func_)(int dim_grid, int x)):
																																				dimension(dimension_), name(name_), S_c(S_c_), type_of_BDC(type_of_BDC_), delta_t(delta_t_), source_e(standard_no_source), source_h(standard_no_source) {
	file_Q_magnetic = new ofstream("/Users/Hannah_Pinson/Documents/vub/FDTD_TimeDependent_Metasurface/progress/test_ADE_staggered/Q_magnetic.txt", std::ios::app);
	file_Q_electric = new ofstream("/Users/Hannah_Pinson/Documents/vub/FDTD_TimeDependent_Metasurface/progress/test_ADE_staggered/Q_electric.txt", std::ios::app);

	/*initialization of e and h*/
	for (int i = 0; i < dimension; i++){  /* number of magnetic nodes = number of electric nodes - 1; in this way the grid ends at both sides with an electric node (so e.g. the boundary conditions can be PEC at both sides)*/
		e.push_back ((*init_func_)(dimension,i));
	}
	for (int i = 0; i < dimension-1; i++){ //number of magnetic nodes = number of electric nodes - 1
		h.push_back((*init_func_)(dimension,i));
	}
	/*set electric field and magnetic update coefficients*/
	// no layers present (yet) -> simulation::static_setup adds layers and modifies update coefficients
	for (int i = 0; i < dimension; i++){
		ch_e.push_back(S_c/imp0);
		ch_h.push_back(1);
		ce_e.push_back(1);
		ce_h.push_back(imp0*S_c);
	}
};


//ÑÑÑÑÑÑ static setup


void Simulation::static_setup(Source& source_e_, Source& source_h_, vector <Static_Material>& static_layers_){
	source_e = source_e_;
	source_h = source_h_;
	static_layers = static_layers_;

	for (int layer_number = 0; layer_number < static_layers.size(); layer_number++){
		Static_Material layer = static_layers[layer_number];
		double loss_h = layer.sigma_h * delta_t / (2*layer.mu);
		double loss_e = layer.sigma_e * delta_t / (2*layer.eps);

		for (int i = layer.start_node; i < layer.end_node + 1; i++){
			ch_h[i] = (1.0 - loss_h) /(1.0 + loss_h);
			ch_e[i] = (S_c/(layer.mu*c)) / (1.0 + loss_h);
			ce_e[i] = (1.0 - loss_e) /(1.0 + loss_e);
			ce_h[i] = (S_c / (layer.eps*c))/ (1.0 + loss_e);
		}
	}
}

//ÑÑÑÑÑ time-dependent setup


void Simulation::add_metasurfaces(vector <Dispersive_Metasurface>& metasurfaces_){
	metasurfaces = metasurfaces_;
}


//ÑÑÑÑÑ simulation procedures


void Simulation::update_magnetic_once(double& timestep){

	double q = timestep; //+ 0.5;

	//first update Q according to ADE with current value of the magnetic field
	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){ //for all metasurfaces present
			int location = metasurfaces.at(sheet_number).node;
			metasurfaces.at(sheet_number).sigma_functor_h.calculate_next(h.at(location), delta_t, q);
		}
	}

	//update all nodes, without current terms
	for (int i = 0; i < dimension-1; i++){     //loop from 0 to size-1 : there exists no node[dimension] for the h-field (to make the grid end with an electric node, see: constructor of simulation)
		h[i] = ch_h[i] * h[i] + ch_e[i] * (e[i + 1] - e[i]) ;} // staggered scheme

	//add current terms
	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){ //for all metasurfaces present
			int location = metasurfaces.at(sheet_number).node;
			double J = metasurfaces.at(sheet_number).sigma_functor_h.calculate_J_term() /mu_0;
			value_to_file(J, file_Q_magnetic);
			h.at(location) -= J;
			cout << "subtracted from sheet * imp0: " << J * imp0 << endl;
			cout << "h at sheet * imp0: " << h.at(location) * imp0 << endl;
			cout << "_______________" << endl;
			//cout << metasurfaces.at(sheet_number).sigma_functor_h.calculate_J_term()/mu_0 << endl;
		}
	}
}

void Simulation::update_electric_once(double& timestep){

	double q = timestep;

	//first update Q according to ADE with current value of the electric field
	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){ //for all metasurfaces present
			int location = metasurfaces.at(sheet_number).node;
			metasurfaces.at(sheet_number).sigma_functor_e.calculate_next(e.at(location), delta_t, q);
		}
	}

	//update all nodes, without current terms
	for (int i = 1; i < dimension-1; i++){ // loop from 1 to size-1 : outer nodes are used for boundary conditions
		e[i] = (ce_e[i] * e[i]) +  (ce_h[i] * (h[i] - h[i - 1]));} //staggered scheme


	//add current terms
	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){ //for all metasurfaces present
			int location = metasurfaces.at(sheet_number).node;
			double J = metasurfaces.at(sheet_number).sigma_functor_e.calculate_J_term() / eps_0;
			value_to_file(J, file_Q_electric);
			e.at(location) -= J;
			cout << "subtracted from sheet: " << J << endl;
			cout << "e at sheet: " << e.at(location) << endl;
			cout << "_______________" << endl;
		}
	}
}

void Simulation::apply_source(double& timestep, vector <double>& field, Source& source){

	if (source.type_of_source == 'h'){ //hard source  --> =
		field[source.node] = (*(source.ptr_source_functor))(timestep);
	}
	else if (source.type_of_source == 'a'){ //additive source --> +=
		field[source.node] += (*(source.ptr_source_functor))(timestep);
	}
	else if (source.type_of_source == 's'){  //subtractive source --> - =
		field[source.node] -= (*(source.ptr_source_functor))(timestep);
	}
	else {
		string msg = name;
		msg.append(" :: apply_source:: source not of type 'h' (hard), 'a' (additive), nor 's' (subtractive) ");
		throw msg ;
	}
}

void Simulation::update_system_once(double& timestep){

	if (type_of_BDC == "PEC"){ //perfect electric conductor
		e[0] = 0;
		e[dimension] = 0;
	}
	else if (type_of_BDC == "ABC"){ //simple absorbing boundary conditions
		e[0] = e[1];
		e[dimension-1] = e[dimension-2];
	}
	else{
		string error_msg = name;
		error_msg.append("::type of boundary condition not recognized");
		throw error_msg;
	}

	apply_source(timestep, e, source_e);
	update_electric_once(timestep);
	//timestep += 0.5;
	apply_source(timestep, h, source_h);
	update_magnetic_once(timestep);
}


/*ÑÑÑÑÑÑÑ--------------- write & output procedures-----------------*/



template <class T>
string convert_to_string (T& t){
	ostringstream convert;   // stream used for the conversion
	convert << t;      // insert the textual representation of t (e.g. an integer) in the characters in the stream
	return convert.str();
}

string Simulation::generate_filename(int number, string field){
	//returns string "nameOfSimulation_number.txt" with e.g. number = number of node
	return name + "_" + field + "_" + convert_to_string(number) + ".txt";
}

vector<ofstream*> Simulation::generate_node_output_files(string output_path, vector<int>& readout_nodes, string field){
	vector<ofstream*> output_files;
	for(int i = 0; i < readout_nodes.size(); i++){ //fills the vector with (pointers to) files with a name corresponding to the readout node
		string total = output_path + generate_filename(readout_nodes[i], field);
		output_files.push_back(new ofstream(total.c_str(), std::ios::app));
	}
	return output_files;
}

ofstream* Simulation::generate_snapshot_file(string output_path, int framenumber, string field){
	string total = output_path + "w_" + generate_filename(framenumber, field);
	return new ofstream (total.c_str(), std::ios::app);
}

void Simulation::value_to_file(double& value, ofstream* outputfile){
	if (outputfile->is_open()){
		if (abs(value) > lower_value_limit){
			(*outputfile) << value << endl;}
		else
			(*outputfile) << 0 << endl;}
	else{
		string error_msg = name + "::value_to_file:: error opening file";
		throw error_msg;
	}
}

void Simulation::data_to_snapshot_file(ofstream* file, string field){
	if (field == "e") {
		for (int  i = 0; i < dimension; i++)
			value_to_file(e[i], file);
	}
	else if (field == "h") {
		for (int  i = 0; i < dimension; i++)
			value_to_file(h[i], file);
	}
	else{
		string error_msg = name + "::data_to_snapshot_file:: field type not recognized";
		throw error_msg;
	}
}


//ÑÑÑÑÑÑ total simulation procedure

void Simulation::simulate(int total_timesteps,string output_path, int snapshotmodulo_electric, int snapshotmodulo_magnetic, vector<int> readout_nodes_electric, vector<int> readout_nodes_magnetic ){


	vector<ofstream*> node_output_files_electric = generate_node_output_files(output_path, readout_nodes_electric, "e");
	vector<ofstream*> node_output_files_magnetic = generate_node_output_files(output_path, readout_nodes_magnetic, "h");

	int frame_electric = 0;
	int frame_magnetic = 0;

	for (int time = 0; time < total_timesteps; time++) {

		double timestep = (double) time;

		update_system_once(timestep);

		//write electric field values (time domain) to file
		//will not save data if size of vector of specified nodes = 0 (=> no time domain information for electric field saved)
		for(int i = 0; i < node_output_files_electric.size(); i++){
			double electric_field_value = e[readout_nodes_electric[i]]; //value of electric field at readout_node i of all readout nodes
			value_to_file(electric_field_value,node_output_files_electric[i]);
		}

		//write magnetic field values (time domain) to file
		//will not save data if size of vector of specified nodes = 0 (=> no time domain information for magnetic field saved)
		for(int i = 0; i < node_output_files_magnetic.size(); i++){
			double magnetic_field_value = h[readout_nodes_magnetic[i]]; //value of magnetic field at readout_node i of all readout nodes
			value_to_file(magnetic_field_value,node_output_files_magnetic[i]);
		}

		//write electric field values (spatial domain) to file if timestep = multiple of snapshotmodulo
		// will not save data if specified snapshotmodulo = -1
		if (snapshotmodulo_electric != -1 && time % snapshotmodulo_electric == 0) {
			ofstream* snapshotfile = generate_snapshot_file(output_path, frame_electric, "e");
			data_to_snapshot_file(snapshotfile, "e");
			frame_electric++;
		}

		//write magnetic field values (spatial domain) to file if timestep = multiple of snapshotmodulo
		// will not save data if specified snapshotmodulo = -1
		if (snapshotmodulo_magnetic != -1 && time % snapshotmodulo_magnetic == 0) {
			ofstream* snapshotfile = generate_snapshot_file(output_path, frame_magnetic, "h");
			data_to_snapshot_file(snapshotfile, "h");
			frame_magnetic++;
		}
	}
}


/*---------------- MAIN ------------------- */


int main(){




	/*! ------general simulation parameters------*/

	/*
	 * grid_scaling example:
	 *
	 *  0    0									0 0 0 0
	 *   			grid_scaling = 2 -> 		0 0 0 0
	 *         			    					0 0 0 0
	 *  0    0									0 0 0 0
	 */

	double grid_scaling = 1; //! - grid_scaling: accuracy of simulation for same real-world distance, time, frequency => delta_t/grid_scaling, delta_x/grid_scaling; #timesteps * grid_scaling, #nodes * grid_scaling
	double S_c = 1; //! - S_c: courant number = delta_t * c / delta_x
	double delta_t = pow(10.0, -12)/grid_scaling; //! - delta_t: discrete unit of time, in seconds
	double delta_x = c*delta_t/S_c; //! - delta_x: discrete unit of space, in meter
	int zones = 12;//15; //! - zones: number of zones in simulation, used for easy specification of positions
	int nodes_per_zone = 1000 * grid_scaling; //! - nodes_per_zone: nodes per zone, used for changing the number of nodes without changing relative positions of materials and sources
	int dimension = zones*nodes_per_zone; //! - dimension: total number of spatial nodes


	/*! ------source------*/

	double initial_frequency = 1 * pow(10.0, 9);  //1.3 * pow(10.0, 9); //! - initial_frequency: frequency of the source signal, in Hz
	double initial_wavelength = c/initial_frequency; //! - initial_wavelength: wavelength of the source signal, in m
	double N_lambda = initial_wavelength/delta_x; //! - N_lambda: number of nodes per wavelength (free space and source signal)
	double T = 1/initial_frequency; //! - T: period of the source signal
	double amplitude_factor = 1 ; //! - amplitude_factor: amplitude factor for the source signal. Multiply by 1/impedance for magnetic sources.
	double width = 2700;//dimension/5;
	double delay = 5000;//2*width;
	int source_position = 5.05*nodes_per_zone;//6*nodes_per_zone; //! -source_position : node where the source is located
	char source_type = 'a'; //'a' = additive

	Gaussian_packet_functor gaussian_packet(S_c, N_lambda, delay, width, amplitude_factor);
	Source gaussian_source_e(source_position, source_type, &gaussian_packet);
	Gaussian_packet_functor gaussian_packet_h(S_c, N_lambda, delay, width, amplitude_factor/imp0);
	Source gaussian_source_h(source_position, source_type, &gaussian_packet_h);
	//Source no_source_h = standard_no_source;





	/*! ------boundary: perfectly matched layers (PMLs)------*/

	double PML_eps = 1.5*eps_0; //! -PML_eps : permittivity of the PML
	double PML_mu = 1.5*mu_0; //! -PML_mu : permeability of the PML
	double PML_sigma_e = PML_eps*0.01/(delta_t*grid_scaling); //! -PML_sigma_e : electric conductivity of the PML
	double PML_sigma_h = PML_mu*0.01/(delta_t*grid_scaling); //! -PML_sigma_h : magnetic conductivity of the PML. Needed make the PML impedance-matched, allowed since the PML is a non-physical region.
	int startPML_left = 0; //! -startPML_left : node where the left PML starts, usually node 0.
	int endPML_left = 2*nodes_per_zone;
	int startPML_right = dimension - 2*nodes_per_zone;
	int endPML_right = dimension;
	Static_Material PML_left(PML_sigma_e, PML_sigma_h, PML_eps, PML_mu, startPML_left, endPML_left);
	Static_Material PML_right(PML_sigma_e, PML_sigma_h, PML_eps, PML_mu, startPML_right, endPML_right);


	//static material objects (here: PMLs)
	vector<Static_Material> static_layers;
	static_layers.push_back(PML_left);
	static_layers.push_back(PML_right);



	/*! ------metasurface with time-varying conductivity------*/



	//functors for time-varying conductivities

	int m = 5;//res[m_index]; //number of resonances
	double g = 0.09;//gamma[g_index];
	double a = 0.85;//a_array[a_index];

	double beta = T*2; //?????
	Linear_Time_Variation t0_functor(g, beta);

	double ksi_factor_e = 2/imp0 / delta_x * delta_t;
	double ksi_factor_h = 2*imp0 / delta_x * delta_t;

	Sigma_ADE sigma_e(m, a, ksi_factor_e, &t0_functor);
	Sigma_ADE sigma_h(m, a, ksi_factor_h, &t0_functor);

	vector<Dispersive_Metasurface> metasurfaces;
	int sheet_position = 5.10*nodes_per_zone;
	Dispersive_Metasurface sheet(sheet_position, sigma_e, sigma_h, eps_0, mu_0);
	metasurfaces.push_back(sheet);


	/*! ------simulation objects and setup------*/

	//constructing free space simulation object & setup
	Simulation free_space("free_space_test", dimension, S_c, delta_t, "ABC", &all_zero);
	free_space.static_setup(gaussian_source_e, gaussian_source_h, static_layers);

	//constructing metasurface simulation object & setup
	string name = "test";//"res_" + convert_to_string(m) + "_gamma_" + convert_to_string(g) + "_a_" + convert_to_string(a);
	Simulation metasurface(name, dimension, S_c, delta_t, "ABC", &all_zero);
	metasurface.static_setup(gaussian_source_e, gaussian_source_h, static_layers);
	metasurface.add_metasurfaces(metasurfaces);


	//metasurface.add_metasurfaces(metasurfaces);


	//ÑÑÑÑÑ simulate
	try{ //try --> if error is 'thrown', rest of program skipped until 'catch' block catches the error

		/* consecutive simulations of the same Simulation object influence one another */

		//output options
		vector<int> sample_nodes_e;
		int readout_node = 5.15*nodes_per_zone;//9*nodes_per_zone;
		sample_nodes_e.push_back(readout_node);

		vector<int> no_sample_nodes_h;

		int snapshotmodulo_electric = -1;//500*grid_scaling;
		int snapshotmodulo_magnetic = -1;//500*grid_scaling;
		string output_path = "/Users/Hannah_Pinson/Documents/vub/FDTD_TimeDependent_Metasurface/progress/test_ADE_staggered/";

		//string output_path = "/user/vginis/Hannah/";

		//simulation
		cout << "starting simulation " << endl;//<< name << endl;
		clock_t start;
		start = clock();

		double timesteps = 25000 * grid_scaling / S_c ;
		metasurface.simulate(timesteps, output_path, snapshotmodulo_electric, snapshotmodulo_magnetic, sample_nodes_e, no_sample_nodes_h);
		free_space.simulate(timesteps, output_path, snapshotmodulo_electric, snapshotmodulo_magnetic, sample_nodes_e, no_sample_nodes_h);
		double diff = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout << "-----finished in "<< diff/60 << " minutes." << endl;

	}


	catch (string error_message){
		cerr << error_message << endl;
	}

	return 0;
}
