#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <utility>
#include <limits>

using namespace std;

struct Point {
	double x, y;
	Point() : x(0), y(0) {}
	Point(double x, double y) : x(x), y(y) {}
	Point operator+(const Point& rhs) const {
		return Point(x + rhs.x, y + rhs.y);
	}
	Point operator-(const Point& rhs) const {
		return Point(x - rhs.x, y - rhs.y);
	}
	Point operator*(const double rhs) const {
		return Point(x * rhs, y * rhs);
	}
	Point operator/(const double rhs) const {
		return Point(x / rhs, y / rhs);
	}
	double operator*(const Point& rhs) const {
		return (x * rhs.x) + (y * rhs.y);
	}
	double cross(const Point& rhs) const {
		return (x * rhs.y) - (y * rhs.x);
	}
	double length() const {
		return sqrt(x*x + y*y);
	}
	Point perp() const {
		return Point(-y, x);
	}
	Point norm() const {
		return Point(x, y) / length();
	}
	bool operator<(const Point& rhs) const {
		if (x == rhs.x) return y < rhs.y;
		else return x < rhs.x;
	}
	bool operator>(const Point& rhs) const {
		if (x == rhs.x) return y > rhs.y;
		else return x > rhs.x;
	}
	bool operator==(const Point& rhs) const {
		return x == rhs.x && y == rhs.y;
	}
	bool operator!=(const Point& rhs) const {
		return x != rhs.x || y != rhs.y;
	}
	friend ostream& operator<<(ostream& os, const Point& rhs) {
		os << "(" << rhs.x << ", " << rhs.y << ")";
		return os;
	}
};

enum inter_t {parallel = -3, collinear = -2, projective = -1, real = 0, overlap = 1};
struct Segment {
	Point s, e;
	Segment(Point s, Point e) : s(s), e(e) {}
	Segment(double sx, double sy, double ex, double ey) : s(Point(sx, sy)), e(Point(ex, ey)) {}
	pair<Point, inter_t> intersection(const Segment& rhs) const {
		Segment l(min(s,e), max(s,e)), r(min(rhs.s,rhs.e), max(rhs.s,rhs.e));
		Point lse = l.e - l.s, rse = r.e - r.s;
		Point diff = l.s - r.s;
		if (lse.x * rse.y == rse.x * lse.y) {
			if (diff.x * rse.y == rse.x * diff.y) {
				if (l.e.x >= r.s.x && r.e.x >= l.s.x) {
					return pair<Point, inter_t>(max(l.s, r.s), overlap);
				} else {
					return pair<Point, inter_t>(max(l.s, r.s), collinear);
				}
			} else {
				return pair<Point, inter_t>(Point(), parallel);
			}
		} else {
			double lt = (l.s-r.s).cross(rse)/rse.cross(lse);
			double rt = (l.s-r.s).cross(lse)/rse.cross(lse);
			Point inter = max(l.s + lse * lt, r.s + rse * rt);
			bool proj = lt * rt > min(lt, rt) || lt + rt < max(lt, rt);
			if (proj) {
				return pair<Point, inter_t>(inter, projective);
			} else {
				return pair<Point, inter_t>(inter, real);
			}
		}
	}
	double length() const {
		return (e - s).length();
	}
	bool operator==(const Segment& rhs) const {
		return s == rhs.s && e == rhs.e;
	}
	friend ostream& operator<<(ostream& os, const Segment& rhs) {
		os << rhs.s << " - " << rhs.e;
		return os;
	}
};

struct PlayerController {
	virtual vector<double> tick(vector<double> state) = 0;
};

struct PongGame {
	const int tickrate = 60;
	int max_score = 1;
	int left_score = 0, right_score = 0;
	int left_returns = 0, right_returns = 0;
	const double length = 400, width = 300, paddle_width = width/5, paddle_max_vel = width/tickrate;
	const Point ball_start_vel = Point(length/tickrate, length/tickrate);
	Point ball_pos = Point(0, 0), ball_vel = ball_start_vel;
	double left_pos = 0, left_vel = 0;
	double right_pos = 0, right_vel = 0;
	PlayerController & left, & right;
	PongGame(PlayerController & left, PlayerController & right) : left(left), right(right) {}
	void tick() {
		// Get left velocity from controller
		double left_cont = left.tick(vector<double>({
			2*ball_pos.x/length, 2*ball_pos.y/width,
			2*ball_vel.x/length, 2*ball_vel.y/width,
			2*left_pos/width, 2*left_vel/width,
			2*right_pos/width, 2*right_vel/width
		})).front() * paddle_max_vel;

		// Get right velocity from controller
		double right_cont = right.tick(vector<double>({
			-2*ball_pos.x/length, -2*ball_pos.y/width,
			-2*ball_vel.x/length, -2*ball_vel.y/width,
			-2*right_pos/width, -2*right_vel/width,
			-2*left_pos/width, -2*left_vel/width
		})).front() * -1 * paddle_max_vel;

		// Update paddle positions and velocities
		left_vel = left_cont;
		left_pos += left_vel;
		if (left_pos < paddle_width/2 - width/2) left_pos = paddle_width/2 - width/2;
		if (left_pos > width/2 - paddle_width/2) left_pos = width/2 - paddle_width/2;

		right_vel = right_cont;
		right_pos += right_vel;
		if (right_pos < paddle_width/2 - width/2) right_pos = paddle_width/2 - width/2;
		if (right_pos > width/2 - paddle_width/2) right_pos = width/2 - paddle_width/2;

		// Prepare geometry segments
		Segment mvmt(ball_pos, ball_pos + ball_vel);
		Segment upr_wall(Point(-length/2, -width/2), Point(length/2, -width/2));
		Segment lwr_wall(Point(-length/2, width/2), Point(length/2, width/2));
		Segment left_seg(Point(-length/2, left_pos - paddle_width/2), Point(-length/2, left_pos + paddle_width/2));
		Segment right_seg(Point(length/2, right_pos - paddle_width/2), Point(length/2, right_pos + paddle_width/2));

		vector<Segment> collidees = {upr_wall, lwr_wall, left_seg, right_seg};
		for (int hit = 0; hit >= 0;) {
			hit = -1;
			double best = numeric_limits<double>::max();
			for (int i = 0; i < collidees.size(); ++i) {
				auto inter = mvmt.intersection(collidees[i]);
				if (inter.second >= real && inter.first != mvmt.s) {
					if ((inter.first - mvmt.s).length() < best) {
						best = (inter.first - mvmt.s).length();
						hit = i;
					}
				}
			}
			if (hit >= 0) {
				Segment seg = collidees[hit];
				collidees.erase(collidees.begin() + hit);
				if (seg == left_seg) ++left_returns;
				if (seg == right_seg) ++right_returns;
				auto inter = mvmt.intersection(seg);
				mvmt.s = inter.first;
				Point perpv = (seg.e - seg.s).perp().norm();
				mvmt.e = mvmt.e - perpv * 2.0 * (perpv * (mvmt.e - mvmt.s));
				ball_vel = ball_vel - perpv * 2.0 * (perpv * ball_vel);
			}
		}

		ball_pos = mvmt.e;

		if (abs(ball_pos.x) > length/2) {
			if (ball_pos.x < 0) {
				++right_score;
				ball_vel = ball_start_vel * -1;
			} else {
				++left_score;
				ball_vel = ball_start_vel;
			}
			ball_pos = Point(0, 0);
			left_pos = 0;
			right_pos = 0;
		}

		assert(abs(ball_pos.y) <= width/2); // We should now still be inside bounds
	}
	pair<int, int> simulate() {
		int timelimit = 2 * 8 * max_score * length / abs(ball_start_vel.x);
		while (max(left_score, right_score) < max_score && timelimit > 0) {
			tick();
			--timelimit;
		}
		return pair<int, int>(left_score, right_score);
	}
};