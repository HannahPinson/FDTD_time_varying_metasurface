/*
 * output_methods.h
 *
 *  Created on: Mar 27, 2015
 *      Author: Hannah_Pinson
 */

#ifndef OUTPUT_METHODS_H_
#define OUTPUT_METHODS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

const double lower_value_limit = pow(10.0,-30); //smaller values are taken to be zero

string convert_to_string(int number);

string generate_filename(string name, int number, string field);
ofstream* generate_file(string basename);
vector<ofstream*> generate_node_output_files(string output_path, string name, string field, vector<int>& readout_nodes);
ofstream* generate_snapshot_file(string output_path, int framenumber, string name, string field);

void value_to_file(double& value, ofstream* file);
void data_to_snapshot_file(ofstream* file, vector<double>& e_or_h);


#endif /* OUTPUT_METHODS_H_ */
