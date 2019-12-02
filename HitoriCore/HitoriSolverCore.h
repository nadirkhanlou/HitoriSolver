#pragma once
#include <string>
#include <vector>

namespace HitoriSolverCore {

enum class SearchType { GreedyBfs, AStar, HillClimbing, SimulatedAnnealing };

// Custom less and greater operators for  State scores. Used in priority queues.
auto stateLessCompartor = [](State lhs, State rhs) {
  return lhs.score < rhs.score;
};
auto stateGreaterCompartor = [](State lhs, State rhs) {
  return lhs.score < rhs.score;
};

class State {
 public:
  float score;  // The score based on which the state will be placed in the
                // priority queue. How the score is calculated is dependent on
                // the search algorithm (searchType).
  float costSoFar;  // The actual cost to reach this state.

  State();
  ~State();

  bool IsGoal();
  std::vector<State>* Successor(SearchType searchType,
                                float (*heuristic)(State));
};

class HitoriSolver {
 public:
  HitoriSolver(const char* filePath);
  ~HitoriSolver();

  State GreedyBfs(State initialState, float (*heuristic)(State));
  State AStar(State initialState, float (*heuristic)(State));
  State HillClimbing(State initialState, float (*heuristic)(State));
  State RandomizedHillClimbing(State initialState, float (*heuristic)(State));
  State KStartHillClimbing(State initialState, float (*heuristic)(State));
  State KStartRandomizedHillClimbing(State initialState,
                                     float (*heuristic)(State));
  State SimulatedAnnealing(
      State initialState,
      float (*heuristic)(State));  // @TODO: Get scheduler as a parameter

 private:
  int** _gameBoard;
};
}  // namespace HitoriSolverCore