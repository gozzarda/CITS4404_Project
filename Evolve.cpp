#include "Pong.cpp"
#include "NeuralNet.cpp"
#include <cmath>
#include <chrono>
#include <thread>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <map>
#include <limits>

using namespace std;

const vector<int> layers = {8, 16, 1};
const int gen_size = 32;
const int gen_limit = 64;

struct NeuroPlayer : PlayerController {
	vector<int> layers;
	vector<double> weights;
	NeuroPlayer(vector<int> layers, vector<double> weights) : layers(layers), weights(weights) {}
	vector<double> tick(vector<double> state) override {
		return evaluate_neural_net(layers, weights, state);
	}
};

void randomize_genomes(vector<vector<double>> & genomes) {
	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	uniform_real_distribution<double> distribution(-1.0, 1.0);
	auto rng = bind(distribution, generator);
	for (auto & genome : genomes) {
		generate(genome.begin(), genome.end(), rng);
	}
}

vector<vector<double>> fittest(int keep, const vector<vector<double>> & population) {
	map<int, pair<int, int>> scores;
	for (int li = 0; li < population.size(); ++li) {
		for (int ri = 0; ri < population.size(); ++ri) {
			if (li == ri) continue;
			NeuroPlayer left(layers, population[li]);
			NeuroPlayer right(layers, population[ri]);
			PongGame game(left, right);
			game.simulate();
			if (game.left_score > game.right_score) {
				++scores[li].second;
			} else {
				++scores[ri].second;
			}
			scores[li].first += game.left_returns;
			scores[ri].first += game.right_returns;
		}
	}
	vector<pair<pair<int, int>, int>> rankings;
	for (auto & kv : scores) {
		rankings.push_back(pair<pair<int, int>, int>(kv.second, kv.first));
	}
	partial_sort(rankings.begin(), rankings.begin() + keep, rankings.end(), greater<pair<pair<int, int>, int>>());
	vector<vector<double>> selected;
	for (int i = 0; i < keep; ++i) {
		selected.push_back(population[rankings[i].second]);
	}
	cerr << " Best score: <" << rankings.front().first.first << ", " << rankings.front().first.second << ">.";
	return selected;
}

/*vector<double> crossover(const vector<double> & lp, const vector<double> & rp) {
	if (lp.size() != rp.size()) throw "crossover: parent genome lengths do not match";
	vector<double> result(lp);
	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	uniform_int_distribution<int> selector(0, result.size());
	for (int i = selector(generator); i < result.size(); ++i) {
		result[i] = rp[i];
	}
	return result;
}*/

vector<double> crossover(const vector<double> & lp, const vector<double> & rp) {
	if (lp.size() != rp.size()) throw "crossover: parent genome lengths do not match";
	vector<double> result(lp);
	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	uniform_int_distribution<int> selector(0, 1);
	for (int i = 0; i < result.size(); ++i) {
		result[i] = (selector(generator) == 1) ? lp[i] : rp[i];
	}
	return result;
}

vector<double> mutation(const vector<double> & parent) {
	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	uniform_int_distribution<int> selector(0, parent.size() - 1);
	normal_distribution<double> distribution(0.0, 1.0);
	vector<double> result = parent;
	result[selector(generator)] += distribution(generator);
	return result;
}

int main() {
	const int keep = (int) sqrt(gen_size);

	vector<vector<double>> population(gen_size, vector<double>(layers_to_weights(layers)));
	randomize_genomes(population);

	for (int gen = 0; gen < gen_limit; ++gen) {
		cerr << "Evaulating generation " << gen << "...";

		population = fittest(keep, population);

		for (int i = 0; i < keep; ++i) {
			for (int j = i + 1; j < keep; ++j) {
				population.push_back(crossover(population[i], population[j]));
			}
		}

		default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
		uniform_int_distribution<int> distribution(0, keep-1);
		while (population.size() < gen_size) {
			population.push_back(mutation(population[distribution(generator)]));
		}

		cerr << " Done." << endl;
	}

	while (true) {
		cerr << "Press Enter to run simulation" << endl;
		cin.ignore();
		cerr << endl << endl << endl;
		NeuroPlayer left(layers, population[0]);
		NeuroPlayer right(layers, population[0]);
		PongGame pong(left, right);
		while (max(pong.left_score, pong.right_score) < pong.max_score) {
			cerr << " ";
			for (int i = 0; i < (int) pong.length / 10; ++i)
				cerr << "=";
			cerr << " " << endl;
			for (int i = (int) -pong.width / 20; i <= (int) pong.width / 20; ++i) {
				if (abs(pong.left_pos/10 - i) <= pong.paddle_width / 20)
					cerr << "|";
				else
					cerr << " ";
				for (int j = 0; j < (int) pong.length / 10; ++j) {
					if ((int) (pong.ball_pos.y / 10) == i && (int) (pong.ball_pos.x / 10 + pong.length / 20) == j) {
						cerr << "O";
					} else {
						cerr << " ";
					}
				}
				if (abs(pong.right_pos/10 - i) <= pong.paddle_width / 20)
					cerr << "|";
				else
					cerr << " ";
				cerr << endl;
			}
			cerr << " ";
			for (int i = 0; i < (int) pong.length / 10; ++i)
				cerr << "=";
			cerr << " " << endl;
			cerr << "ball_pos: " << pong.ball_pos << "\tball_vel: " << pong.ball_vel << endl;
			cerr << "left_pos: " << pong.left_pos << "\tleft_vel: " << pong.left_vel << endl;
			cerr << "right_pos: " << pong.right_pos << "\tright_vel: " << pong.right_vel << endl;
			cerr << "left_score: " << pong.left_score << "\tright_score: " << pong.right_score << endl;
			pong.tick();
			//cin.ignore();
			this_thread::sleep_for(chrono::milliseconds(1000/pong.tickrate));
		}
		cerr << "SCORES: " << pong.left_score << ", " << pong.right_score << endl;
	}
	return 0;
}