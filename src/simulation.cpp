/*
 * simulation.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#include "simulation.h"
#include "output_methods.h"

using namespace std;


//ÑÑÑÑÑ simulation constructor

Simulation::Simulation(string name_, int dimension_, double S_c_, double delta_t_, string type_of_BDC_, Source source_e_, Source source_h_):
				dimension(dimension_), name(name_), S_c(S_c_), type_of_BDC(type_of_BDC_), delta_t(delta_t_), source_e(source_e_), source_h(source_h_){

	/*initialization of e and h*/

	for (int i = 0; i < dimension; i++){  /* number of magnetic nodes = number of electric nodes - 1; in this way the grid ends at both sides with an electric node (so e.g. the boundary conditions can be PEC at both sides)*/
		e.push_back (0);
	}
	for (int i = 0; i < dimension-1; i++){ //number of magnetic nodes = number of electric nodes - 1
		h.push_back(0);
	}

	/*set electric field and magnetic update coefficients*/
	// no layers present (yet) -> simulation::static_setup adds layers and modifies update coefficients
	for (int i = 0; i < dimension; i++){
		ch_e.push_back(S_c/imp0);
		ch_h.push_back(1);
		ce_h.push_back(S_c*imp0);
		ce_e.push_back(1);
	}
};


//ÑÑÑÑÑÑ static setup:


void Simulation::static_setup(vector <Static_Material>& static_layers_){

	static_layers = static_layers_;

	for (int layer_number = 0; layer_number < static_layers.size(); layer_number++){
		Static_Material layer = static_layers[layer_number];
		double loss_h = layer.sigma_h * delta_t / (2*layer.mu);
		double loss_e = layer.sigma_e * delta_t / (2*layer.eps);

		for (int i = layer.start_node; i < layer.end_node + 1; i++){
			ch_e[i] = (S_c/(layer.mu*c)) / (1.0 + loss_h);
			ch_h[i] = (1.0 - loss_h) /(1.0 + loss_h);
			ce_h[i] = (S_c / (layer.eps*c))/ (1.0 + loss_e);
			ce_e[i] = (1.0 - loss_e) /(1.0 + loss_e);
		}
	}
}

//ÑÑÑÑÑ time-dependent setup


void Simulation::add_metasurfaces(vector <Dispersive_Metasurface>& metasurfaces_){
	metasurfaces = metasurfaces_;
}


//ÑÑÑÑÑ simulation procedures


double Simulation::calculate_convolution_term_e(double q, Dispersive_Metasurface& sheet){

	double coefficient = - ( pow(delta_t,2) / eps0);
	double half_integer_time_correction = (sheet.sigma_functor_e(q*delta_t, q*delta_t) * sheet.saved_e.at(0))/2;
	double summation = 0;

	for (int i = 0; i < q; i++){
		summation += sheet.sigma_functor_e((q)*delta_t, (i)*delta_t) * sheet.saved_e.at(q-i);
	}
	double total = coefficient*(summation + half_integer_time_correction);
	return total;
}


double Simulation::calculate_convolution_term_h(double q, Dispersive_Metasurface& sheet){

	double coefficient = - ( pow(delta_t,2) / mu0);
	double integer_time_correction = (sheet.sigma_functor_h( q*delta_t, q*delta_t ) * sheet.saved_h.at(0))/2;
	double summation = 0;

	for (int i = 0; i < q; i++){ // i < q-1
		summation += sheet.sigma_functor_h((q)*delta_t, (i)*delta_t)  * sheet.saved_h.at(q-i);
	}
	double total = coefficient*(summation + integer_time_correction);
	return total;
}


void Simulation::update_magnetic_once(double& timestep){


	for (int i = 0; i < dimension-1; i++){     //loop from 0 to size-1 : there exists no node[dimension] for the h-field (to make the grid end with an electric node, see: constructor of simulation)
		h[i] = ch_h[i] * h[i] + ch_e[i] * (e[i+1] - e[i]);
	}


	//all nodes updated, without convolution terms
	//add convolution terms and save local field
	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){ //for all metasurfaces present
			int location = metasurfaces.at(sheet_number).node;
			h.at(location) += calculate_convolution_term_h(timestep, metasurfaces.at(sheet_number)); // minus sign in coefficient of calculation
			double local_field = (h.at(location-1) + h.at(location)) / 2;
			metasurfaces.at(sheet_number).saved_h.push_back(local_field); //save local field
		}
	}
}


void Simulation::update_electric_once(double& timestep){

	for (int i = 1; i < dimension-1; i++){ // loop from 1 to size-1 : outer nodes are used for boundary conditions
		e[i] = ce_e[i] * e[i] + ce_h[i] * (h[i] - h[i-1]);
	}

	//all nodes updated, without convolution terms
	//add convolution terms and save local field

	if (metasurfaces.size() != 0){
		for (int sheet_number = 0; sheet_number < metasurfaces.size(); sheet_number++){
			int location = metasurfaces.at(sheet_number).node;
			e.at(location) += calculate_convolution_term_e(timestep, metasurfaces.at(sheet_number));
			double local_field = (e.at(location) + e.at(location+1)) / 2;
			//double local_field = e.at(location);
			metasurfaces.at(sheet_number).saved_e.push_back(local_field); //save local field
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

void Simulation::update_system_once(int& timestep){

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

	double timestep_magnetic = timestep;
	double timestep_electric = timestep;
	//double timestep_electric = timestep + 0.5; // half integer time correction is accounted for in the calculation of the convolution integral

	apply_source(timestep_magnetic, h, source_h);
	update_magnetic_once(timestep_magnetic);
	apply_source(timestep_electric, e, source_e);
	update_electric_once(timestep_electric);


}



//ÑÑÑÑÑÑ total simulation procedure

void Simulation::simulate(int total_timesteps,string output_path, int snapshotmodulo_electric, int snapshotmodulo_magnetic, vector<int> readout_nodes_electric, vector<int> readout_nodes_magnetic ){


	vector<ofstream*> node_output_files_electric = generate_node_output_files(output_path, name, "e", readout_nodes_electric);
	vector<ofstream*> node_output_files_magnetic = generate_node_output_files(output_path, name, "h", readout_nodes_magnetic);

	int frame_electric = 0;
	int frame_magnetic = 0;

	for (int time = 1; time < total_timesteps; time++) {

		update_system_once(time);

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
			ofstream* snapshotfile = generate_snapshot_file(output_path, frame_electric, name, "e");
			data_to_snapshot_file(snapshotfile, e);
			frame_electric++;
		}

		//write magnetic field values (spatial domain) to file if timestep = multiple of snapshotmodulo
		// will not save data if specified snapshotmodulo = -1
		if (snapshotmodulo_magnetic != -1 && time % snapshotmodulo_magnetic == 0) {
			ofstream* snapshotfile = generate_snapshot_file(output_path, frame_magnetic, name, "h");
			data_to_snapshot_file(snapshotfile, h);
			frame_magnetic++;
		}
	}
}






