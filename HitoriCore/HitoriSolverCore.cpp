#include "HitoriSolverCore.h"

#include <fstream>
#include <iostream>

//#define PERMUTATION_SUCCESSOR

namespace HitoriSolverCore {
// Softmax filter that turns a vector of real valued, positive numbers (which
// are in this case heuristic values for successsors) into a probability
// distribution. The 'alpha' factor can increase or decrease the preference of a
// smaller valued heuristic (more attractive choice) over a larger valued one
// (less attractive choice)
void Softmax(std::vector<double>& scores, double alpha = 1.0f) {
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
}

State::State()
    : _dimension(0), _level(0), score(0), costSoFar(0), _blackedOut(nullptr) {}

State::State(int dimension)
    : _dimension(dimension),
      _level(0),
      score(0),
      costSoFar(0),
      _blackedOut(new bool*[_dimension]) {
  for (int i = 0; i < _dimension; ++i) {
    _blackedOut[i] = new bool[_dimension];
    for (int j = 0; j < _dimension; ++j) {
      _blackedOut[i][j] = false;
    }
  }
}

State::State(const State& other)
    : _level(other._level), score(other.score), costSoFar(other.costSoFar) {
  this->_dimension = other._dimension;
  this->_blackedOut = new bool*[_dimension];
  for (int i = 0; i < _dimension; ++i) {
    this->_blackedOut[i] = new bool[_dimension];
    std::copy(other._blackedOut[i], other._blackedOut[i] + _dimension,
              this->_blackedOut[i]);
  }
}

State::State(State&& other) noexcept
    : _dimension(other._dimension),
      _level(other._level),
      score(other.score),
      costSoFar(other.costSoFar) {
  this->_blackedOut = other._blackedOut;
  other._blackedOut = nullptr;
  other._dimension = 0;
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
    delete[] this->_blackedOut;
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

State& State::operator=(State&& other) noexcept {
  if (&other != this) {
    this->score = other.score;
    this->costSoFar = other.costSoFar;
    this->_level = other._level;
    for (int i = 0; i < _dimension; ++i) {
      delete[] this->_blackedOut[i];
    }
    delete[] this->_blackedOut;
    this->_dimension = other._dimension;
    this->_blackedOut = other._blackedOut;
    other._blackedOut = nullptr;
    other._dimension = 0;
  }
  return *this;
}

HitoriSolver::HitoriSolver(const char* filePath) {
  std::ifstream input;
  input.open(filePath);

  input >> _dimension;
  _gameBoard = new int*[_dimension];
  _whitedOut = new bool*[_dimension];

  for (int i = 0; i < _dimension; i++) {
    _gameBoard[i] = new int[_dimension];
    _whitedOut[i] = new bool[_dimension];
    for (int j = 0; j < _dimension; j++) {
      input >> _gameBoard[i][j];
      _whitedOut[i][j] = false;
    }
  }
  input.close();
}

HitoriSolver::~HitoriSolver() {}

State HitoriSolver::InitialState() { return State(_dimension); }

bool HitoriSolver::IsGoal(State& s) {
  assert(this->_dimension == s._dimension);
  // We know that while using PermutationSuccessor, we will not approach an
  // answer without going through all rows
#ifdef PERMUTATION_SUCCESSOR
  if (s._level < _dimension) return false;
#endif

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

std::vector<State> HitoriSolver::PermutationSuccessor(
    State& currentState, SearchType searchType,
    double (*heuristic)(const State&, int**)) {

    //This successor creates all permutations of shading blocks in a row which are feasible

  std::vector<State> successors;
  if (currentState._level >= _dimension) return successors; // if the row is the last one it doesnt have any successor


  // clean row is a row which there is no need to shade any blocks in it
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

  //NShadeGenerator generates permutations of all numbers in a row taken N at a time
  
  //Here in this loop NShadeGenerator is called with parameter N = i
  
  //If PreProcess method is called at start a min and a max bound can be found for i
  //otherwise it will start from 1 to half the size of the row
  for (int i =
           std::max(_rowConflicts ? _rowConflicts[currentState._level] : 1, 1);
       i <=
       std::min((_dimension + 1) / 2,
                _rowMaxShade ? _rowMaxShade[currentState._level] : _dimension);
       i++) {
    NShadeGenerator(i, successors, currentState);
  }

  //Depends on search type the score will have different value
  //In GreedyBFS and HillClimbing it is equal to heuristic of each state but in other algorithms
  //it is needed to consider g(n) -which here is called costSoFar- as well
  
  if (searchType == SearchType::GreedyBfs ||
      searchType == SearchType::HillClimbing) {
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      (*it).score = heuristic(*it, _gameBoard);
    }
  } else {
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      (*it).score = ((*it).costSoFar) + heuristic(*it, _gameBoard);
    }
  }

  return successors;
}

std::vector<State> HitoriSolver::NextConfilctSuccessor(
    State& currentState, SearchType searchType,
    double (*heuristic)(const State&, int**)) {
  std::vector<State> successors;

  //This successor has 2 options as successor in each candidate block, if none of them be
  //feasible it will iterate to next candidate block

  //Candidate blocks are blocks which have one or more conflict with a block in same row or column
  
  //As it is mentioned earlier there are 2 options, one is shading the candidate block and other
  //option is shade all other blocks which have conflict with the candidate block
  
  //In this successor _level stores the last block checked in predecessor nodes 

  for (int I = currentState._level; I < _dimension * _dimension; I++) {
    div_t resultD = div(I, _dimension);
    int i = resultD.quot;
    int j = resultD.rem;

    //if the block is shaded in predecessors or it is not candidate to shade 
    if (currentState._blackedOut[i][j] || _whitedOut[i][j]) continue;

    //find conflicts
    std::vector<std::pair<int, int>> confilcts;
    for (int k = 0; k < _dimension; ++k) {
      if (k != j && !currentState._blackedOut[i][k] &&
          _gameBoard[i][j] == _gameBoard[i][k]) {
        confilcts.push_back(std::make_pair(i, k));
      }
      if (k != i && !currentState._blackedOut[k][j] &&
          _gameBoard[i][j] == _gameBoard[k][j]) {
        confilcts.push_back(std::make_pair(k, j));
      }
    }

    if (!confilcts.empty()) {
      State onlyOneShaded = currentState;
      onlyOneShaded._level = I;
      onlyOneShaded._blackedOut[i][j] = true;
      onlyOneShaded.costSoFar += 1;

      //check if candidate block is feasible to shade
      if (IsFeasible(onlyOneShaded)) {
        successors.push_back(onlyOneShaded);
      }

      State shadedConfilcts = currentState;
      shadedConfilcts._level = I;

      for (auto it = confilcts.begin(); it != confilcts.end(); ++it) {
        shadedConfilcts._blackedOut[(*it).first][(*it).second] = true;
        shadedConfilcts.costSoFar += 1;
      }

      //check if shading all other blocks in conflict with the candidate block is feasible
      if (IsFeasible(shadedConfilcts)) {
        shadedConfilcts._level++;
        successors.push_back(shadedConfilcts);
      }

      if (!successors.empty()) break;
    }
  }

  // Depends on search type the score will have different value
  // In GreedyBFS and HillClimbing it is equal to heuristic of each state but in other algorithms
  //it is needed to consider g(n) -which here is called costSoFar- as well
  if (searchType == SearchType::GreedyBfs ||
      searchType == SearchType::HillClimbing) {
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      (*it).score = heuristic(*it, _gameBoard);
    }
  } else {
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      (*it).score = ((*it).costSoFar) + heuristic(*it, _gameBoard);
    }
  }

  return successors;
}

std::vector<State> HitoriSolver::Successor(State& currentState,
                                           SearchType searchType,
                                           double (*heuristic)(const State&,
                                                               int**)) {
//Depending on PERMUTATION_SUCCESSOR definition the successor will decide which successor
//method should be called
#ifdef PERMUTATION_SUCCESSOR
  return PermutationSuccessor(currentState, searchType, heuristic);
#endif  // PERMUTATION_SUCCESSOR
#ifndef PERMUTATION_SUCCESSOR
  return NextConfilctSuccessor(currentState, searchType, heuristic);
#endif  // !PERMUTATION_SUCCESSOR
}

bool HitoriSolver::IsFeasible(const State& currentState) {

  //This feasiblity checker is used in NextConflictSuccessor

  div_t resultD = div(currentState._level, _dimension);
  int x = resultD.quot;
  int y = resultD.rem;

  //check feasibility of shading candidate block 
  if (currentState._blackedOut[x][y]) {
    if (!IsFeasibleAdj(currentState, x, y) ||
        !IsFeasibleSurround(currentState, x, y))
      return false;
    else
      return true;
  }

  //Check feasibility when it is desired to shade all other blocks in conficts with
  //candidate block but itself 
  for (int i = 0; i < _dimension; ++i) {
    if (currentState._blackedOut[i][y]) {
      if (!IsFeasibleAdj(currentState, i, y) ||
          !IsFeasibleSurround(currentState, i, y))
        return false;
    }
  }
  for (int i = 0; i < _dimension; ++i) {
    if (currentState._blackedOut[x][i]) {
      if (!IsFeasibleAdj(currentState, x, i) ||
          !IsFeasibleSurround(currentState, x, i))
        return false;
    }
  }

  return true;
}

bool HitoriSolver::IsFeasibleAdj(const State& currentState, int x, int y) {

  //It will retrun false if shading the block cause shaded blocks become adjacent

  if ((x < _dimension - 1 && currentState._blackedOut[x + 1][y]) ||
      (x > 0 && currentState._blackedOut[x - 1][y]) ||
      (y < _dimension - 1 && currentState._blackedOut[x][y + 1]) ||
      (y > 0 && currentState._blackedOut[x][y - 1])) {
    return false;
  }

  return true;
}

bool HitoriSolver::IsFeasibleSurround(const State& currentState, int x, int y) {

  //It will return false if shading the block cause a white block be surrounded by shaded blocks

  if (((x >= _dimension - 2 || currentState._blackedOut[x + 2][y]) &&
       (x < _dimension - 1 &&
        (y > _dimension - 1 || currentState._blackedOut[x + 1][y + 1])) &&
       (x < _dimension - 1 &&
        (y <= 0 || currentState._blackedOut[x + 1][y - 1]))) ||
      ((x <= 1 || currentState._blackedOut[x - 2][y]) &&
       (x > 0 &&
        (y >= _dimension - 1 || currentState._blackedOut[x - 1][y + 1])) &&
       (x > 0 && (y <= 0 || currentState._blackedOut[x - 1][y - 1]))) ||
      ((y >= _dimension - 2 || currentState._blackedOut[x][y + 2]) &&
       (y < _dimension - 1 &&
        (x >= _dimension - 1 || currentState._blackedOut[x + 1][y + 1])) &&
       (y < _dimension - 1 &&
        (x <= 0 || currentState._blackedOut[x - 1][y + 1]))) ||
      ((y <= 1 || currentState._blackedOut[x][y - 2]) &&
       (y > 0 &&
        (x >= _dimension - 1 || currentState._blackedOut[x + 1][y - 1])) &&
       (y > 0 && (x <= 0 || currentState._blackedOut[x - 1][y - 1]))))
    return false;

  return true;
}

bool HitoriSolver::IsFeasible(bool* shaded, const State& currentState) {

  //It will check feasiblity of shading blocks based on permutations in a row

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

  //Consider if conflicts are solved, if not the state is not feasible 
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

void HitoriSolver::NShadeGenerator(int n, std::vector<State>& succcessorStates,
                                   const State& currentState) {
  //NShadeGenerator generates permutations of all blocks in a row taken N at a time recursively
  
    bool* tempShadedFlags = new bool[_dimension];
  for (int i = 0; i < _dimension; i++) {
    tempShadedFlags[i] = false;
  }

  //It is not needed to generate permutations after _dimension - n + 1 because it is transparent
  //can not generate enough to get N
  for (int i = 0; i < _dimension - n + 1; i++) {
    Shade(tempShadedFlags, i, succcessorStates, currentState, n - 1);
  }
}

void HitoriSolver::Shade(bool* shaded, int selectIndex,
                         std::vector<State>& succcessorStates,
                         const State& currentState, int recursiveLevel) {
  if (currentState._blackedOut[currentState._level][selectIndex]) return;

  if (!IsFeasibleSurround(currentState, currentState._level, selectIndex) ||
      !IsFeasibleAdj(currentState, currentState._level, selectIndex))
    return;

  //Select block to shade temporary in scope of this heap, It will undo it before exit the scope
  shaded[selectIndex] = true;

  if (recursiveLevel < 1) {
    if (IsFeasible(shaded, currentState)) {
      State newState = currentState;
      for (int i = 0; i < _dimension; i++) {
        newState._blackedOut[newState._level][i] = shaded[i];
        if (shaded[i]) {
          bool isConflictDecreased = false;
          for (int j = 0; j < _dimension; j++) {
            if (j != currentState._level &&
                _gameBoard[currentState._level][i] == _gameBoard[j][i])
              isConflictDecreased = true;
          }
          if (!isConflictDecreased)
            newState.costSoFar += 1; 
        }
      }
      newState._level++;
      succcessorStates.push_back(newState);
    }
    shaded[selectIndex] = false;
    return;
  }

  //If a block is shaded next one can not be shaded so it is enough to start from
  //selectIndext + 2
  //Also it is enough to generate until _dimension - recursiveLevel + 1 as the same reason as
  //explained in NShadeGenerator method
  for (int i = selectIndex + 2; i < _dimension - recursiveLevel + 1; i++) {
    if (!shaded[i])
      Shade(shaded, i, succcessorStates, currentState, recursiveLevel - 1);
  }
  shaded[selectIndex] = false;
}

void HitoriSolver::PreProcess() {

  //Before running algorithm some blocks can determined as whitedOut, It means there is no need
  //to shaded them to get to the solution
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

  //_rowConfilcts and _rowMaxShade give a bound to generate permutations in PermutationSuccessor
  _rowConflicts = new int[_dimension];
  _rowMaxShade = new int[_dimension];
  for (int L = 0; L < _dimension; L++) {
    bool* conflictsCounter = new bool[_dimension];
    _rowConflicts[L] = 0;
    _rowMaxShade[L] = 0;
    for (int i = 0; i < _dimension; i++) {
      conflictsCounter[i] = false;
    }
    for (int i = 0; i < _dimension; i++) {
      if (!conflictsCounter[i]) {
        for (int j = i + 1; j < _dimension; j++) {
          if (_gameBoard[L][i] ==
              _gameBoard[L][j]) {
            conflictsCounter[j] = true;
            _rowConflicts[L]++;
            _rowMaxShade[L]++;
            if (!conflictsCounter[i]) {
              conflictsCounter[i] = true;
              _rowMaxShade[L]++;
            }
          }
        }
        if (conflictsCounter[i]) continue;

        for (int j = 0; j < _dimension; j++) {
          if (j != i) {
            if (_gameBoard[L][i] == _gameBoard[j][i]) {
              _rowMaxShade[L]++;
              break;
            }
          }
        }

      }
    }
    delete[] conflictsCounter;
  }
  
}

void HitoriSolver::PrintState(const State& currentState) {
  PrintState(currentState, _gameBoard);
}

void HitoriSolver::PrintState(const State& state, int** gameBoard) {
  std::cout << "\n";
  for (int i = 0; i < state._dimension; i++) {
    for (int j = 0; j < state._dimension; j++) {
      if (state._blackedOut[i][j]) {
        std::cout << std::setfill(' ') << std::setw(DIGITS(state._dimension))
                  << '#' << ' ';
      } else {
        std::cout << std::setfill(' ') << std::setw(DIGITS(state._dimension))
                  << gameBoard[i][j] << ' ';
      }
    }
    std::cout << "\n";
  }
}

double HitoriSolver::HeuristicFunction1(const State& currentState,
                                        int** gameBoard) {
  double h = 0;

  //This heuristic function is compatible with the PermutationSuccessor
  //It will find the blocks in next levels in conflicts with the blocks in the row remained
  //after blocks choosen by the successor in the row are shaded

  if (currentState._level < 1) {
    for (int i = 0; i < currentState._dimension; i++) {
      for (int j = currentState._level; j < currentState._dimension; j++) {
        if (!currentState._blackedOut[0][i])
          if (gameBoard[0][i] == gameBoard[j][i]) h++;
      }
    }
    return h;
  }

  for (int i = 0; i < currentState._dimension; i++) {
    for (int j = currentState._level; j < currentState._dimension; j++) {
      if (!currentState._blackedOut[currentState._level - 1][i])
        if (gameBoard[currentState._level - 1][i] == gameBoard[j][i]) h++;
    }
  }
  return h;
}

double HitoriSolver::HeuristicFunction2(const State& currentState,
                                        int** gameBoard) {
  std::set<int> set;

  int row_conflicts = 0;
  // Check all rows for duplicate numbers
  for (int i = 0; i < currentState._dimension; ++i) {
    set.clear();
    for (int j = 0; j < currentState._dimension; ++j) {
      if (!currentState._blackedOut[i][j]) {
        if (set.find(gameBoard[i][j]) != set.end()) {
          ++row_conflicts;
        } else {
          set.insert(gameBoard[i][j]);
        }
      }
    }
  }
  int col_conflicts = 0;
  // Check all columns for duplicate numbers
  for (int j = 0; j < currentState._dimension; ++j) {
    set.clear();
    for (int i = 0; i < currentState._dimension; ++i) {
      if (!currentState._blackedOut[i][j]) {
        if (set.find(gameBoard[i][j]) != set.end()) {
          ++col_conflicts;
        } else {
          set.insert(gameBoard[i][j]);
        }
      }
    }
  }

  return std::max(row_conflicts, col_conflicts);
}

std::vector<double> HitoriSolver::SCurveSchedule(int steps) {
  double alpha = 1.0;
  double beta = 3.0;
  double t;
  double increment;
  t = 0.01;
  increment = (1 - t) / double(steps);
  std::vector<double> temperatures;
  while (t <= 1) {
    double temp = (alpha) / (alpha + pow((t / (1 - t)), (beta)));
    temperatures.push_back(temp);
    t += increment;
  }
  return temperatures;
}

State HitoriSolver::GreedyBfs(State initialState,
                              double (*heuristic)(const State&, int**)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score =
  // heuristic. States with lower scores are higher up in the heap.
  std::priority_queue<State, std::vector<State>, StateGreaterComparator>
      priorityQueue;

  priorityQueue.push(initialState);

  int counter = 0;

  // Search until all possible states are checked
  while (!priorityQueue.empty()) {
    // Pop the state with the lowest score. This will be the state that is to be
    // examined in this iteration.
    State currentState = priorityQueue.top();
    priorityQueue.pop();
    ++counter;

	// If currentState is the goal, return it
    if (IsGoal(currentState)) {
      std::cout << "Nodes checked: " << counter << "\n";
      std::cout << "States in queue: " << priorityQueue.size() << "\n";
      return currentState;
    }
    // currentState was not the goal. So we expand it in our search tree.
    std::vector<State> successors =
        Successor(currentState, SearchType::GreedyBfs, heuristic);
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      priorityQueue.push(*it);
    }
  }

  // The goal could not be found.
  std::cout << "Nodes checked: " << counter << "\n";
  std::cout << "States in queue: " << priorityQueue.size() << "\n";
  return State();
}

State HitoriSolver::AStar(State initialState,
                          double (*heuristic)(const State&, int**)) {
  // Create a new minimum priority queue (min heap) and insert the initial state
  // into it. States with compared to each other based on their score = (cost so
  // far + heuristic). States with lower scores are higher up in the heap.
  std::priority_queue<State, std::vector<State>, StateGreaterComparator>
      priorityQueue;
  priorityQueue.push(initialState);

  int counter = 0;

  // Search until all possible states are checked
  while (!priorityQueue.empty()) {
    // Pop the state with the lowest score. This will be the state that is to be
    // examined in this iteration.
    State currentState = priorityQueue.top();
    priorityQueue.pop();
    ++counter;
    // If currentState is the goal, return it
    if (IsGoal(currentState)) {
      std::cout << "Nodes checked: " << counter << "\n";
      std::cout << "States in queue: " << priorityQueue.size() << "\n";
      return currentState;
    }
    // current_state was not the goal. So we expand it in our search tree.
    std::vector<State> successors =
        Successor(currentState, SearchType::AStar, heuristic);
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      priorityQueue.push(*it);
    }
    successors.clear();
  }
  // The goal could not be found.
  std::cout << "Nodes checked: " << counter << "\n";
  std::cout << "States in queue: " << priorityQueue.size() << "\n";
  return State();
}

State HitoriSolver::SteepestAscentHillClimbing(State initialState,
                                               double (*heuristic)(const State&,
                                                                   int**)) {
  // Initialize a currentState variable with the initial state
  State currentState = initialState;
  currentState.score = heuristic(currentState, _gameBoard);
  int counter = 0;
  while (true) {
    ++counter;
    // Check whether the current state is a goal
    if (IsGoal(currentState)) {
      std::cout << "Nodes checked: " << counter << "\n";
      return currentState;
    }

    // Generate the successors of the current state
    std::vector<State> successors =
        Successor(currentState, SearchType::HillClimbing, heuristic);

    // No successors
    if (successors.begin() == successors.end()) break;

    // Find the state with the best score (lowest heuristic)
    auto bestSuccessor = std::min_element(successors.begin(), successors.end(),
                                          StateLessComparator());
    // If there was a successor with a better score than the current state,
    // continue the search with that successor. Otherwise return the current
    // state as the final result
    if (bestSuccessor->score <= currentState.score)
      currentState = *bestSuccessor;
    else {
      break;
    }
  }

  // Return the current state (which is a local minimum and COULD be a
  // global minimum)
  std::cout << "Nodes checked: " << counter << "\n";
  return currentState;
}

State HitoriSolver::StochasticHillClimbing(State initialState,
                                           double (*heuristic)(const State&,
                                                               int**)) {
  // Initialize a currentState variable with the initial state
  State currentState = initialState;
  currentState.score = heuristic(currentState, _gameBoard);
  int counter = 0;
  while (true) {
    ++counter;
    // Check whether the current state is a goal
    if (IsGoal(currentState)) {
      std::cout << "Nodes checked: " << counter << "\n";
      return currentState;
    }

    // Generate the successors of the current state
    std::vector<State> successors =
        Successor(currentState, SearchType::HillClimbing, heuristic);

    // Now each state will be assigned a probability according to its score
    // (heuristic). The higher the score is, the lower the probability will be.
    // This is done using a softmax function

    // Create a new vector 'probabilities', which initially contains the score
    // of each successor (which are NOT probabilities)
    std::vector<double> probabilities;
    for (auto it = successors.begin(); it != successors.end(); ++it) {
      probabilities.push_back(it->score);
    }

	// No successors
	if (successors.begin() == successors.end()) break;
    
    // Apply softmax to 'probabilities' to get a probability distribution
    Softmax(probabilities);

    // Randomly choose a state based on the probability distribution above.
    double randomNumber = RANDOM;
    State nextState;
    double low = 0;
    for (unsigned int i = 0; i < probabilities.size(); ++i) {
      if (randomNumber >= low && randomNumber < low + probabilities[i]) {
        // The i'th successor (with probability = probabilities[i]) was selected
        nextState = (successors)[i];
        break;
      }
      low += probabilities[i];
    }
    // Check if 'nextState' is better than 'currentState'
    if (nextState.score < currentState.score) {
      currentState = nextState;
    } else {
      break;
    }
  }
  // Return the current state (which is a local minimum and COULD be a
  // global minimum)
  std::cout << "Nodes checked: " << counter << "\n";
  return currentState;
}

State HitoriSolver::KStartStochasticHillClimbing(
    State initialState, double (*heuristic)(const State&, int**)) {
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

  // The goal could not be found.
  return State();
}
State HitoriSolver::SimulatedAnnealing(State initialState,
                                       double (*heuristic)(const State&, int**),
                                       std::vector<double> (*scheduler)(int)) {
  // Get the mapping of time (step) to temperature
  std::vector<double> schedule = scheduler(1000);
  double delta;
  double randomIndex;
  State currentState = initialState;
  currentState.score = heuristic(currentState, _gameBoard);
  // Start the simulation (1000 iterations)
  int counter = 0;
  for (auto T = schedule.begin(); T != schedule.end(); ++T) {
    ++counter;
    std::vector<State> successors =
        Successor(currentState, SearchType::SimulatedAnnealing, heuristic);
	// Check if there are any neighbours for the current state
    if (successors.empty()) {
      break;
    }
	// Randomly choose one of the neighbours
    randomIndex = std::rand() % successors.size();
    delta = successors[randomIndex].score - currentState.score;
	// If the selected neighbour is better, choose it
    if (delta < 0) {
      currentState = successors[randomIndex];
	// If the selected neighbour is not better, choose it with probability p = exp(-(delta / *T))
    } else {
      if (RANDOM < std::exp(-(delta / *T))) {
        currentState = successors[randomIndex];
      } else {
      }
    }
  }
  schedule.clear();

  std::cout << "Nodes checked: " << counter << "\n";
  return currentState;
}
}  // namespace HitoriSolverCore