#include "HitoriSolverCore.h"
#include <fstream>
#include <iostream>

namespace HitoriSolverCore {
State::State() : _dimension(0), score(0), costSoFar(0) {}

State::State(int dimension) : _dimension(dimension), score(0), costSoFar(0) {
  _blackedOut = new bool*[_dimension];
  for (int i = 0; i < _dimension; ++i) {
    _blackedOut[i] = new bool[_dimension]{};  // Initializes all values to false
  }
}

State::State(const State& other)
    : score(other.score), costSoFar(other.costSoFar) {
  for (int i = 0; i < _dimension; ++i) {
    delete[] this->_blackedOut[i];
  }
  delete this->_blackedOut;
  this->_dimension = other._dimension;
  this->_blackedOut = new bool*[_dimension];
  for (int i = 0; i < _dimension; ++i) {
    this->_blackedOut[i] = new bool[_dimension];
    std::copy(other._blackedOut[i], other._blackedOut[i] + _dimension,
              this->_blackedOut[i]);
  }
}

State::State(State&& other)
    : _dimension(other._dimension),
      score(other.score),
      costSoFar(other.costSoFar) {
  this->_blackedOut = other._blackedOut;
  other._blackedOut = nullptr;
}

State::~State() {
  for (int i = 0; i < _dimension; ++i) {
    delete[] _blackedOut[i];
  }
  delete[] _blackedOut;
}

State& State::operator=(const State& other) {
  if (&other != this) {
    this->score = other.score;
    this->costSoFar = other.costSoFar;
    this->_level = other._level;
    for (int i = 0; i < _dimension; ++i) {
      delete[] this->_blackedOut[i];
    }
    delete this->_blackedOut;
    this->_dimension = other._dimension;
    this->_blackedOut = new bool*[_dimension];
    for (int i = 0; i < _dimension; ++i) {
      this->_blackedOut[i] = new bool[_dimension];
      std::copy(other._blackedOut[i], other._blackedOut[i] + _dimension,
                this->_blackedOut[i]);
    }
  }
  return *this;
}

State& State::operator=(State&& other) {
  if (&other != this) {
    this->score = other.score;
    this->costSoFar = other.costSoFar;
    for (int i = 0; i < _dimension; ++i) {
      delete[] this->_blackedOut[i];
    }
    delete this->_blackedOut;
    this->_dimension = other._dimension;
    this->_blackedOut = other._blackedOut;
    other._blackedOut = nullptr;
  }
  return *this;
}

HitoriSolver::HitoriSolver(const char* filePath) {
  std::ifstream input;
  input.open(filePath);

  input >> _dimension;
  _gameBoard = new int*[_dimension];

  for (int i = 0; i < _dimension; i++) {
    _gameBoard[i] = new int[_dimension];
    for (int j = 0; j < _dimension; j++) {
      input >> _gameBoard[i][j];
    }
  }

  input.close();
}

HitoriSolver::~HitoriSolver() {}

bool HitoriSolver::IsGoal(State& s) {
  assert(this->_dimension == s._dimension);

  std::set<int> set;
  int count;

  // Check all rows for duplicate numbers
  for (int i = 0; i < _dimension; ++i) {
    set.clear();
    count = 0;
    for (int j = 0; j < _dimension; ++j) {
      if (!s._blackedOut[i][j]) {
        set.insert(_gameBoard[i][j]);
        ++count;
      }
    }
    // Duplicate numbers found
    if (set.size() != count) {
      return false;
    }
  }
  // Check all columns for duplicate numbers
  for (int j = 0; j < _dimension; ++j) {
    set.clear();
    count = 0;
    for (int i = 0; i < _dimension; ++i) {
      if (!s._blackedOut[i][j]) {
        set.insert(_gameBoard[i][j]);
        ++count;
      }
    }
    // Duplicate numbers found
    if (set.size() != count) {
      return false;
    }
  }
  // Now columns or rows with duplicate numbers
  return true;
}

std::vector<State>* HitoriSolver::Successor(State currentState, SearchType searchType,
                                            double (*heuristic)(State)) {
  std::vector<State> successors;

  bool* cleanRow = new bool[_dimension];
  for (int i = 0; i < _dimension; i++) {
    cleanRow[i] = false;
  }

  if (IsFeasible(cleanRow, currentState)) {
    State cleanRowState = currentState;
    delete[] cleanRowState._blackedOut[currentState._level];
    cleanRowState._blackedOut[currentState._level] = cleanRow;
    cleanRowState._level++;
    successors.push_back(cleanRowState);
  }

  for (int i = 1; i <= (_dimension + 1) / 2; i++) {
    NShadeGenerator(i, successors, currentState);
  }

  return &successors;
}

bool HitoriSolver::IsFeasible(bool* shaded, State currentState) {
  for (int i = 1; i < _dimension - 1; i++) {
    if ((shaded[i] && shaded[i - 1]) || (shaded[i] && shaded[i + 1]))
      return false;
    if (shaded[i] && currentState._level > 0 &&
        (currentState._blackedOut[currentState._level - 1][i]))
      return false;
    if (shaded[i] && currentState._level == 1 &&
        currentState._blackedOut[currentState._level - 1][i + 1] &&
        currentState._blackedOut[currentState._level - 1][i - 1])
      return false;
    if (shaded[i] && currentState._level > 1 &&
        currentState._blackedOut[currentState._level - 2][i] &&
        currentState._blackedOut[currentState._level - 1][i + 1] &&
        currentState._blackedOut[currentState._level - 1][i - 1])
      return false;
  }

  if (shaded[0]) {
    if (currentState._level > 0 &&
        currentState._blackedOut[currentState._level - 1][0])
      return false;
    if (currentState._level == 1 &&
        currentState._blackedOut[currentState._level - 1][1])
      return false;
    if (currentState._level > 1 &&
        currentState._blackedOut[currentState._level - 1][1] &&
        currentState._blackedOut[currentState._level - 2][0])
      return false;
  }

  if (shaded[_dimension - 1]) {
    if (currentState._level > 0 &&
        currentState._blackedOut[currentState._level - 1][_dimension - 1])
      return false;
    if (currentState._level == 1 &&
        currentState._blackedOut[currentState._level - 1][_dimension - 2])
      return false;
    if (currentState._level > 1 &&
        currentState._blackedOut[currentState._level - 1][_dimension - 2] &&
        currentState._blackedOut[currentState._level - 2][_dimension - 1])
      return false;
  }

  for (int i = 0; i < _dimension; i++) {
    if (!shaded[i]) {
      for (int j = i + 1; j < _dimension; j++) {
        if (!shaded[j] && _gameBoard[currentState._level][i] ==
                              _gameBoard[currentState._level][j])
          return false;
      }
      for (int j = currentState._level - 1; j >= 0; j--) {
        if (!currentState._blackedOut[j][i] &&
            _gameBoard[j][i] == _gameBoard[currentState._level][i])
          return false;
      }
    }
  }

  return true;
}

void HitoriSolver::NShadeGenerator(int n, std::vector<State>& succcessorStates ,
                                   State currentState) {
  bool* tempShadedFlags = new bool[_dimension];
  for (int i = 0; i < _dimension; i++) {
    tempShadedFlags[i] = false;
  }
  for (int i = 0; i < _dimension; i++) {
    Shade(tempShadedFlags, i, succcessorStates, currentState, n - 1);
  }
}

void HitoriSolver::Shade(bool* shaded, int selectIndex, std::vector<State>& succcessorStates, State currentState, int recursiveLevel) {
  if (currentState._blackedOut[currentState._level][selectIndex]) return;

  shaded[selectIndex] = true;

  if (recursiveLevel < 1) {
    if (IsFeasible(shaded, currentState)) {
      State newState = currentState;
      for (int i = 0; i < _dimension; i++)
        newState._blackedOut[newState._level][i] = shaded[i];
      newState._level++;
      succcessorStates.push_back(newState);
    }
    shaded[selectIndex] = false;
    return;
  }
  for (int i = selectIndex; i < _dimension; i++) {
    if (!shaded[i])
      Shade(shaded, i, succcessorStates, currentState, recursiveLevel - 1);
  }
  shaded[selectIndex] = false;
}

void HitoriSolver::PreProccess() {
  _whitedOut = new bool*[_dimension];

  for (int i = 0; i < _dimension; i++) {
    _whitedOut[i] = new bool[_dimension];
    for (int j = 0; j < _dimension; j++) {
      _whitedOut[i][j] = true;

      for (int k = 0; k < _dimension; k++) {
        if (_gameBoard[i][j] == _gameBoard[i][k] ||
            _gameBoard[i][j] == _gameBoard[k][j])
          _whitedOut[i][j] = false;
      }
    }
  }

  for (int i = 0; i < _dimension; i++) {
    for (int j = 0; j < _dimension - 2; j++) {
      if (_gameBoard[i][j] == _gameBoard[i][j + 2]) _whitedOut[i][j + 1] = true;
      if (_gameBoard[j][i] == _gameBoard[j + 2][i]) _whitedOut[j + 1][i] = true;
    }
  }
}

void HitoriSolver::PrintState(State state) {
  std::cout << "\n";
  for (int i = 0; i < _dimension; i++) {
    for (int j = 0; j < _dimension; j++) {
      if (state._blackedOut[i][j]) {
        std::cout << "# ";
      } else {
        std::cout << _gameBoard[i][j] << " ";
      }
    }
    std::cout << "\n";
  }
}

State HitoriSolver::GreedyBfs(State initialState, double (*heuristic)(State)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score =
  // heuristic. States with lower scores are higher up in the heap.
  std::priority_queue<State, std::vector<State>, StateGreaterComparator>
      priorityQueue;
  priorityQueue.push(initialState);

  // Search until all possible states are checked
  while (!priorityQueue.empty()) {
    // Pop the state with the lowest score. This will be the state that is to be
    // examined in this iteration.
    State currentState = priorityQueue.top();
    priorityQueue.pop();

    // If currentState is the goal, return it
    if (IsGoal(currentState)) {
      return currentState;
    }
    // currentState was not the goal. So we expand it in our search tree.
    std::vector<State>* successors =
        Successor(currentState, SearchType::GreedyBfs, heuristic);
    for (auto it = successors->begin(); it != successors->end(); ++it) {
      // @TODO: Add a visited states buffer
      if (/* If *it was already visited */ true) {
        priorityQueue.push(*it);
        // @TODO: Add *it to the visited buffer
      }
    }
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
    State currentState = priorityQueue.top();
    priorityQueue.pop();

    // If currentState is the goal, return it
    if (IsGoal(currentState)) {
      return currentState;
    }
    // current_state was not the goal. So we expand it in our search tree.
    std::vector<State>* successors =
        Successor(currentState, SearchType::AStar, heuristic);
    for (auto it = successors->begin(); it != successors->end(); ++it) {
      // @TODO: Add a visited states buffer
      if (/* If *it was already visited */ true) {
        priorityQueue.push(*it);
        // @TODO: Add *it to the visited buffer
      }
    }
    successors->clear();
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
    if (IsGoal(currentState)) {
      return currentState;
    }

    // Generate the successors of the current state
    std::vector<State>* successors =
        Successor(currentState, SearchType::HillClimbing, heuristic);

    // Find the state with the score (heuristic)
    auto bestSuccessor = std::min_element(
        successors->begin(), successors->end(), StateLessComparator());

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
    if (IsGoal(currentState)) {
      return currentState;
    }

    // Generate the successors of the current state
    std::vector<State>* successors =
        Successor(currentState, SearchType::HillClimbing, heuristic);

    // Now each state will be assigned a probability according to its score
    // (heuristic). The higher the score is, the lower the probability will be.
    // This is done using a softmax function

    // Create a new vector 'probabilities', which initially contains the score
    // of each successor (which are NOT probabilities)
    std::vector<double> probabilities;
    for (auto it = successors->begin(); it != successors->end(); ++it) {
      probabilities.push_back(it->score);
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
    if (IsGoal(*it)) {
      return *it;
    }
  }

  return State();
}
State HitoriSolver::SimulatedAnnealing(State initialState,
                                       double (*heuristic)(State),
                                       std::vector<double>* (*scheduler)(int)) {
  return State();
}
}  // namespace HitoriSolverCore