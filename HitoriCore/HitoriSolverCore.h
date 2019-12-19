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
#include <iomanip>

#define DIGITS(num) (num > 0 ? (int)std::log10((double)num) + 1 : 1)

#define RANDOM static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX)

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
  std::vector<State> PermutationSuccessor(State&, SearchType,
                               double (*)(const State&, int**));
  std::vector<State> NextConfilctSuccessor(State&, SearchType,
                                          double (*)(const State&, int**));
  std::vector<State> Successor(State&, SearchType,
                                double (*)(const State&, int**));

  bool IsFeasible(bool*, const State&);
  bool IsFeasible(const State&);
  bool IsFeasibleAdj(const State&, int, int);
  bool IsFeasibleSurround(const State&, int, int);
  void NShadeGenerator(int, std::vector<State>&, const State&);
  void Shade(bool*, int, std::vector<State>&, const State&, int);
  void PreProccess();
  void PrintState(const State&);

  static void PrintState(const State&, int**);
  static double HeuristicFunction1(const State&, int**);
  static double HeuristicFunction2(const State&, int**);
  static double HeuristicFunction3(const State&, int**);

  State GreedyBfs(State, double (*)(const State&, int**));
  State AStar(State, double (*)(const State&, int**));
  State SteepestAscentHillClimbing(State, double (*)(const State&, int**));
  State StochasticHillClimbing(State, double (*)(const State&, int**));
  State KStartSteepestAscentHillClimbing(State,
                                         double (*)(const State&, int**));
  State KStartStochasticHillClimbing(State, double (*)(const State&, int**));
  State SimulatedAnnealing(
      State, double (*)(const State&, int**),
      std::vector<double>* (*)(int));  // @TODO: Get scheduler as a parameter
 private:
  int _dimension;
  int** _gameBoard;
  bool** _whitedOut;
  int* _rowConflicts;
  int* _rowMaxShade;
};
}  // namespace HitoriSolverCore