#pragma once
#include <string>
#include <vector>
#include <math>
#include <random>
#include <omp.h>

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

class State {
 public:
  // @TODO: Rename 'score' into something more meaningful. Something with a
  // closer meaning to cost!
  double score;  // The score based on which the state will be placed in the
                // priority queue. How the score is calculated is dependent on
                // the search algorithm (searchType).
  double costSoFar;  // The actual cost to reach this state.

  State();
  ~State();

  bool IsGoal();
  std::vector<State>* Successor(SearchType searchType,
                                double (*heuristic)(State));
};

class HitoriSolver {
 public:
  HitoriSolver(const char* filePath);
  ~HitoriSolver();

  State GreedyBfs(State initialState, double (*heuristic)(State));
  State AStar(State initialState, double (*heuristic)(State));
  State SteepestAscentHillClimbing(State initialState, double (*heuristic)(State));
  State StochasticHillClimbing(State initialState, double (*heuristic)(State));
  State KStartSteepestAscentHillClimbing(State initialState,
                                         double (*heuristic)(State));
  State KStartStochasticHillClimbing(State initialState,
                                     double (*heuristic)(State));
  State SimulatedAnnealing(
      State initialState,
      double (*heuristic)(State));  // @TODO: Get scheduler as a parameter

 private:
  int** _gameBoard;
};
}  // namespace HitoriSolverCore