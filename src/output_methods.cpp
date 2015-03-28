/*
 * output_methods.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: Hannah_Pinson
 */


#include "output_methods.h"
#include <sstream>



/*ÑÑÑÑÑÑÑ--------------- write & output procedures-----------------*/


string convert_to_string (int number){
	return static_cast<ostringstream*>( &(ostringstream() << number) )->str();
}

string generate_filename(string name, int number, string field){
	//returns string "nameOfSimulation_number.txt" with e.g. number = number of node
	return name + "_" + field + "_" + convert_to_string(number) + ".txt";
}

vector<ofstream*> generate_node_output_files(string output_path, string name, string field, vector<int>& readout_nodes){
	vector<ofstream*> output_files;
	for(int i = 0; i < readout_nodes.size(); i++){ //fills the vector with (pointers to) files with a name corresponding to the readout node
		string total = output_path + generate_filename(name, readout_nodes[i], field);
		output_files.push_back(new ofstream(total.c_str(), std::ios::app));
	}
	return output_files;
}

ofstream* generate_snapshot_file(string output_path, int framenumber, string name, string field){
	string total = output_path + "w_" + generate_filename(name, framenumber, field);
	return new ofstream (total.c_str(), std::ios::app);
}

void value_to_file(double& value, ofstream* outputfile){
	if (outputfile->is_open()){
		if (abs(value) > lower_value_limit){
			(*outputfile) << value << endl;}
		else
			(*outputfile) << 0 << endl;}
	else{
		string error_msg = "::value_to_file:: error opening file";
		throw error_msg;
	}
}

void data_to_snapshot_file(ofstream* file, vector<double>& e_or_h ){
	for (int  i = 0; i < e_or_h.size(); i++)
		value_to_file(e_or_h[i], file);
}


