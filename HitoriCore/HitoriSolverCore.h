#pragma once
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <queue>
#include <random>
#include <string>
#include <vector>
#include <set>
#include <assert.h>

namespace HitoriSolverCore {

enum class SearchType { GreedyBfs, AStar, HillClimbing, SimulatedAnnealing };

// Custom less and greater operators for State scores. Used in priority queues
// and elsewhere.
class StateLessComparator {
 public:
  bool operator()(const State& lhs, const State& rhs) {
    return lhs.score < rhs.score;
  }
};
class StateGreaterComparator {
 public:
  bool operator()(const State& lhs, const State& rhs) {
    return lhs.score > rhs.score;
  }
};

// Softmax filter that turns a vector of real valued, positive numbers (which
// are in this case heuristic values for successsors) into a probability
// distribution. The 'alpha' factor can increase or decrease the preference of a
// smaller valued heuristic (more attractive choice) over a larger valued one
// (less attractive choice)
void softmax(std::vector<double>& scores, double alpha = 1.0f) {
  double sum = 0.f;
  for (auto it = scores.begin(); it != scores.end(); ++it) {
    // Since we want higher scores to have less probability, we use
    // exp(-x) instead of exp(x - MaximumScore).
    *it = std::exp(-*it * alpha);
    sum += *it;
  }
  for (auto it = scores.begin(); it != scores.end(); ++it) {
    *it /= sum;
  }
};

class State {
 public:
  // @TODO: Rename 'score' into something more meaningful. Something with a
  // closer meaning to cost!
  double score;  // The score based on which the state will be placed in the
                 // priority queue. How the score is calculated is dependent on
                 // the search algorithm (searchType).
  double costSoFar;  // The actual cost to reach this state.

  State();
  State(int);
  State(const State&);
  State(State&&);
  ~State();

  State& operator=(const State&);
  State& operator=(State&&);

 private:
  int _dimension;
  bool** _blackedOut;
  int _level;

  // @TODO: Remove friend class declaration, only IsGoal remains a friend
  friend class HitoriSolver;
  friend bool HitoriSolver::IsGoal(State&);
};

class HitoriSolver {
 public:
  HitoriSolver(const char*);
  ~HitoriSolver();

  bool IsGoal(State&);
  std::vector<State>* Successor(State, SearchType, double (*)(State));

  bool IsFeasible(bool*, State);
  void NShadeGenerator(int, std::vector<State>&, State);
  void Shade(bool*, int, std::vector<State>&, State, int);
  void PreProccess();
  void PrintState(State);

  State GreedyBfs(State, double (*)(State));
  State AStar(State, double (*)(State));
  State SteepestAscentHillClimbing(State, double (*)(State));
  State StochasticHillClimbing(State, double (*)(State));
  State KStartSteepestAscentHillClimbing(State, double (*)(State));
  State KStartStochasticHillClimbing(State, double (*)(State));
  State SimulatedAnnealing(
      State, double (*)(State),
      std::vector<double>* (*)(int));  // @TODO: Get scheduler as a parameter
 private:
  int _dimension;
  int** _gameBoard;
  bool** _whitedOut;
};
}  // namespace HitoriSolverCore