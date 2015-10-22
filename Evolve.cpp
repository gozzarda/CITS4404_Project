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

using namespace std;

const vector<int> layers = {8, 8, 4, 2, 1};
const int gen_size = 128;

struct NeuroPlayer : PlayerController {
	vector<int> layers;
	vector<double> weights;
	NeuroPlayer(vector<int> layers, vector<double> weights) : layers(layers), weights(weights) {}
	vector<double> tick(vector<double> state) override {
		return evaluate_neural_net(layers, weights, state);
	}
};

void randomize_genomes(vector<vector<double>> & genomes) {
	default_random_engine generator;
	uniform_real_distribution<double> distribution(-1.0, 1.0);
	auto rng = bind(distribution, generator);
	for (auto & genome : genomes) {
		generate(genome.begin(), genome.end(), rng);
	}
}

vector<int> selection(int keep, const vector<int> & population, const vector<vector<double>> & genomes) {
	map<int, pair<int, int>> scores;
	for (int lp : population) {
		for (int rp : population) {
			if (lp == rp) continue;
			NeuroPlayer left(layers, genomes[lp]);
			NeuroPlayer right(layers, genomes[rp]);
			PongGame game(left, right);
			auto result = game.simulate();
			if (result.first > result.second) {
				++scores[lp].first;
			} else {
				++scores[rp].first;
			}
			scores[lp].second += result.first - result.second;
			scores[rp].second += result.second - result.first;
		}
	}
	vector<pair<pair<int, int>, int>> rankings;
	for (auto & kv : scores) {
		rankings.push_back(pair<pair<int, int>, int>(kv.second, kv.first));
	}
	partial_sort(rankings.begin(), rankings.begin() + keep, rankings.end(), greater<pair<pair<int, int>, int>>());
	vector<int> selected;
	for (int i = 0; i < keep; ++i) {
		selected.push_back(rankings[i].second);
	}
	return selected;
}

vector<double> crossover(const vector<double> & lp, const vector<double> & rp) {
	if (lp.size() != rp.size()) throw "crossover: parent genome lengths do not match";
	vector<double> result(lp.size());
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < result.size(); ++i) {
		result[i] = lp[i] + distribution(generator) * (rp[i] - lp[i]);
	}
	return result;
}

vector<double> mutation(const vector<double> & parent) {
	default_random_engine generator;
	normal_distribution<double> distribution(0.0, 0.5);
	vector<double> result = parent;
	for (auto & gene : result) {
		gene += distribution(generator);
	}
	return result;
}

int main() {
	vector<vector<double>> genomes(gen_size, vector<double>(layers_to_weights(layers)));
	randomize_genomes(genomes);

	vector<vector<int>> generations(1);
	for (int i = 0; i < gen_size; ++i) {
		generations.back().push_back(i);
	}

	for (int gen = 0; gen < 128; ++gen) {
		cout << "Evaulating generation " << gen << "...";
		cout.flush();

		int keep = (int) sqrt(generations.back().size());
		generations.push_back(selection(keep, generations.back(), genomes));

		for (int i = 0; i < keep; ++i) {
			for (int j = i + 1; j < keep; ++j) {
				generations.back().push_back(genomes.size());
				genomes.push_back(crossover(genomes[generations.back()[i]], genomes[generations.back()[j]]));
			}
		}

		default_random_engine generator;
		uniform_int_distribution<int> distribution(0, keep-1);
		while (generations.back().size() < gen_size) {
			int r = distribution(generator);
			generations.back().push_back(genomes.size());
			genomes.push_back(mutation(genomes[generations.back()[r]]));
		}

		cout << " Done." << endl;
	}

	cout << "Press Enter to run simulation" << endl;
	cin.ignore();

	NeuroPlayer left(layers, genomes[generations.back().front()]);
	NeuroPlayer right(layers, genomes[generations.back().front()]);
	PongGame pong(left, right);
	while (max(pong.left_score, pong.right_score) < pong.max_score) {
		cout << " ";
		for (int i = 0; i < (int) pong.length / 10; ++i)
			cout << "=";
		cout << " " << endl;
		for (int i = (int) -pong.width / 20; i <= (int) pong.width / 20; ++i) {
			if (abs(pong.left_pos/10 - i) <= pong.paddle_width / 20)
				cout << "|";
			else
				cout << " ";
			for (int j = 0; j < (int) pong.length / 10; ++j) {
				if ((int) (pong.ball_pos.y / 10) == i && (int) (pong.ball_pos.x / 10 + pong.length / 20) == j) {
					cout << "O";
				} else {
					cout << " ";
				}
			}
			if (abs(pong.right_pos/10 - i) <= pong.paddle_width / 20)
				cout << "|";
			else
				cout << " ";
			cout << endl;
		}
		cout << " ";
		for (int i = 0; i < (int) pong.length / 10; ++i)
			cout << "=";
		cout << " " << endl;
		cout << "ball_pos: " << pong.ball_pos << "\tball_vel: " << pong.ball_vel << endl;
		cout << "left_pos: " << pong.left_pos << "\tleft_vel: " << pong.left_vel << endl;
		cout << "right_pos: " << pong.right_pos << "\tright_vel: " << pong.right_vel << endl;
		cout << "left_score: " << pong.left_score << "\tright_score: " << pong.right_score << endl;
		pong.tick();
		//cin.ignore();
		this_thread::sleep_for(chrono::milliseconds(1000/pong.tickrate));
	}
	cout << "SCORES: " << pong.left_score << ", " << pong.right_score << endl;
	return 0;
}