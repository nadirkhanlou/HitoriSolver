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
  State(State&&) noexcept;
  ~State();

  State& operator=(const State&);
  State& operator=(State&&) noexcept;

 private:
  int _dimension;
  bool** _blackedOut;
  int _level;

  // @TODO: Remove friend class declaration, only IsGoal remains a friend
  friend class HitoriSolver;
};

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

class HitoriSolver {
 public:
  HitoriSolver(const char*);
  ~HitoriSolver();

  State InitialState();
  bool IsGoal(State&);
  std::vector<State> Successor(State&, SearchType, double (*)(State));

  bool IsFeasible(bool*, const State&);
  void NShadeGenerator(int, std::vector<State>&, const State&);
  void Shade(bool*, int, std::vector<State>&, const State&, int);
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