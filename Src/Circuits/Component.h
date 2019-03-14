#pragma once
#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

// Component Elements
struct Component
{
	string name;
	complex<double> voltage;
	complex<double> current;
	complex<double> Z;
	complex<double> Y;
	int node1, node2, VS_num;
	double factor;
};
