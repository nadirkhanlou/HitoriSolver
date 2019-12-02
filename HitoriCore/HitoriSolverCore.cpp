#include "HitoriSolverCore.h"
#include <algorithm>
#include <queue>

namespace HitoriSolverCore {
State::State() {}

State::~State() {}

bool State::IsGoal() { return true; }

std::vector<State>* State::Successor(SearchType searchType,
                                     double (*heuristic)(State)) {
  return {};
}

HitoriSolver::HitoriSolver(const char* filePath) {}

HitoriSolver::~HitoriSolver() {}

State HitoriSolver::GreedyBfs(State initialState, double (*heuristic)(State)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score =
  // heuristic. States with lower scores are higher up in the heap.
  std::priority_queue<
      std::priority_queue<int, std::vector<State>, std::greater<State>>,
      std::vector<State>, StateGreaterComparator>
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

State HitoriSolver::AStar(State initialState, double (*heuristic)(State)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score = (cost so
  // far + heuristic). States with lower scores are higher up in the heap.
  std::priority_queue<State, std::vector<State>, StateGreaterComparator>
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

State HitoriSolver::SteepestAscentHillClimbing(State initialState,
                                 double (*heuristic)(State)) {
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
        successors->begin(), successors->end(), StateLessComparator);

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

State HitoriSolver::StochasticHillClimbing(State initialState,
                                           double (*heuristic)(State)) {
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

    // Now each state will be assigned a probability according to its score
    // (heuristic). The higher the score is, the lower the probability will be.
    // This is done using a softmax function
    auto softmax = [](std::vector<double>& scores) {
      double sum = 0.f;
      for (auto it = scores.begin(); it != scores.end(); ++it) {
        // Since we want higher scores to have less probability, we use
        // exp(-x) instead of exp(x - MaximumScore).
        *it = std::exp(-*it);
        sum += *it;
      }
      for (auto it = scores.begin(); it != scores.end(); ++it) {
        *it /= sum;
      }
    };

    // Create a new vector 'probabilities', which initially contains the score
    // of each successor (which are NOT probabilities)
    std::vector<double> probabilities;
    for (auto it = successors->begin(); successors->end(); ++it) {
      probabilities.push_back(it->score)
    }
    // Apply softmax to 'probabilities' to get a probability distribution
    softmax(probabilities);

    // Randomly choose a state based on the probability distribution above.
    double randomNumber =
        static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

    State nextState;
    double low = 0;
    for (unsigned int i = 0; i < probabilities.size(); ++i) {
      if (randomNumber >= low && randomNumber < low + probabilities[i]) {
        // The i'th successor (with probability = probabilities[i]) was selected
        nextState = (*successors)[i];
        break;
      }
      low += probabilities[i];
    }

    // Check if 'nextState' is better than 'currentState'
    if (currentState.score < nextState.score) {
      break;
    } else {
      currentState = nextState;
    }
  }

  // Return the current state (which is a local minimum and NOT a global
  // minimum)
  return currentState;
}

State HitoriSolver::KStartSteepestAscentHillClimbing(
    State initialState, double (*heuristic)(State)) {
  return State();
}

State HitoriSolver::KStartStochasticHillClimbing(State initialState,
                                                 double (*heuristic)(State)) {
  // Use the maximum number of threads available
  omp_set_num_threads(omp_get_max_threads());

  // Create a vector the size of the number of threads available
  int numOfThreads = omp_get_num_threads();
  std::vector<State> results(numOfThreads);

#pragma omp parallel
  {
    // Get the thread number
    int threadNum = omp_get_thread_num();
    results[threadNum] = StochasticHillClimbing(initialState, heuristic);

#pragma omp barrier
    {}
  }

  // Check the results from each thread and return the first one that is a goal
  for (auto it = results.begin(); it != results.end(); ++it) {
    if (it->IsGoal()) {
      return *it;
    }
  }

  return State();
}
}  // namespace HitoriSolverCore