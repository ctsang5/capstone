#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <limits>
#include <optional>
#include <stdexcept>


// Define a Matrix struct that stores a 2D grid of double values. 
struct Matrix 
{
		int rows;
		int cols;
		std::vector<double> data;

		Matrix() : rows(0) , cols(0) {}
		Matrix(int row, int col, double fill = 0.0)
				: rows(row), cols(col), data(row * col, fill) {}

		double& operator()(int i, int j) { return data[i * cols + j]; }
		double operator()(int i, int j) const { return data[i * cols + j]; }
};

// Define a struct for simulate_AMS and AMS.
struct simResult
{
		Matrix S;
		Matrix V;
		std::string scheme;		// might remove l8ter because might not implement heston models
}

return List::create(Named("S") = Sfull);


struct AMSResult 
{
		double price;
		double standard_error;
}

