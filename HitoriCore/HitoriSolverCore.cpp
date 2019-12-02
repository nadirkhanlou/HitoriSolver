#include <algorithm>
#include <queue>
#include "HitoriSolverCore.h"

namespace HitoriSolverCore {
State::State() {}

State::~State() {}

bool State::IsGoal() { return true; }

std::vector<State>* State::Successor(SearchType searchType,
                                     float (*heuristic)(State)) {
  return {};
}

HitoriSolver::HitoriSolver(const char* filePath) {}

HitoriSolver::~HitoriSolver() {}

State HitoriSolver::GreedyBfs(State initialState, float (*heuristic)(State)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score =
  // heuristic. States with lower scores are higher up in the heap.
  std::priority_queue<
      std::priority_queue<int, std::vector<State>, std::greater<State>>,
      std::vector<State>, std::greater<State>>
      priorityQueue;
  priorityQueue.push(initialState);

  // Search until all possible states are checked
  while (!priorityQueue.empty()) {
    // Pop the state with the lowest score. This will be the state that is to be
    // examined in this iteration.
    State current_state = priorityQueue.top();
    priorityQueue.pop();

    // If current_state is the goal, return it
    if (current_state.IsGoal()) {
      return current_state;
    }
    // current_state was not the goal. So we expand it in our search tree.
    std::vector<State>* successors =
        current_state.Successor(SearchType::GreedyBfs, heuristic);
    for (auto it = successors->begin(); it != successors->end(); ++it) {
      // @TODO: Add a visited states buffer
      if (/* If *it was already visited */ true) {
        priorityQueue.push(*it);
        // @TODO: Add *it to the visited buffer
      }
    }
    successors->clear();
    delete successors;
  }

  // The goal could not be found.
  return State();
}

State HitoriSolver::AStar(State initialState, float (*heuristic)(State)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score = (cost so
  // far + heuristic). States with lower scores are higher up in the heap.
  std::priority_queue<State, std::vector<State>, std::greater<State>>
      priorityQueue;
  priorityQueue.push(initialState);

  // Search until all possible states are checked
  while (!priorityQueue.empty()) {
    // Pop the state with the lowest score. This will be the state that is to be
    // examined in this iteration.
    State current_state = priorityQueue.top();
    priorityQueue.pop();

    // If current_state is the goal, return it
    if (current_state.IsGoal()) {
      return current_state;
    }
    // current_state was not the goal. So we expand it in our search tree.
    std::vector<State>* successors =
        current_state.Successor(SearchType::AStar, heuristic);
    for (auto it = successors->begin(); it != successors->end(); ++it) {
      // @TODO: Add a visited states buffer
      if (/* If *it was already visited */ true) {
        priorityQueue.push(*it);
        // @TODO: Add *it to the visited buffer
      }
    }
    successors->clear();
    delete successors;
  }

  // The goal could not be found.
  return State();
}

State HitoriSolver::HillClimbing(State initialState,
                                 float (*heuristic)(State)) {
  // Initialize a currentState variable with the initial state
  State currentState = initialState;
  while (true) {
    // Check whether the current state is a goal
    if (currentState.IsGoal()) {
      return currentState;
    }

    // Generate the successors of the current state
    std::vector<State>* successors =
        currentState.Successor(SearchType::HillClimbing, heuristic);

    // Find the state with the score (heuristic)
    auto bestSuccessor = std::min_element(
        successors->begin(), successors->end(), stateLessCompartor);

    // If there was a successor with a better score than the current tate,
    // continue the search with that successor. Otherwise
    if (bestSuccessor->score <= currentState.score)
      currentState = *bestSuccessor;
    else {
      break;
    }
  }

  // Return the current state (which is a local minimum and NOT a global
  // minimum)
  return currentState;
}
}  // namespace HitoriSolverCore