#include <vector>
#include <algorithm>

using namespace std;

struct NeuralNet {
	vector<int> layers; // Layer sizes of neural net including input and output layers
	vector<double> weights; // Per neuron weights to previous layer, plus one for unity weight

	// Construct NeuralNet
	NeuralNet(vector<int> layers, vector<double> weights) : layers(layers), weights(weights) {
		int total = 0;
		for(int i = 1; i < layers.size(); ++i)
			total += (layers[i-1] + 1) * layers[i];
		if (total != weights.size()) throw "NeuralNet layer size and weights mismatch";
	}

	// Activation function used for each neuron
	static double activation_function(double x) {
		if (x < -1.0) return -1.0;
		else if (x > 1.0) return 1.0;
		else return x;
	}

	// Given a vector of inputs (must be of same size as input layer), evaluates the net and returns the values of the output layer
	vector<double> evaluate(vector<double> inputs) {
		if (inputs.size() != layers.front()) throw "NeuralNet input size mismatch";
		vector<double> prev = inputs;
		auto weight = weights.begin();
		for (auto layer = ++layers.begin(); layer != layers.end(); ++layer) {
			for (auto el : prev)
				cerr << " " << 
			vector<double> curr(*layer, 0.0);
			// Calculate weighted input for each neuron in this layer
			for (double & neuron : curr) {
				for (double value : prev)
					neuron += *(weight++) * value;
				neuron += *(weight++);
			}
			// Apply activation function to curr and overwrite prev
			prev.resize(curr.size());
			transform(curr.begin(), curr.end(), prev.begin(), activation_function);
		}
		// prev currently holds the value of the output neurons, so return it
		return prev;
	}
};