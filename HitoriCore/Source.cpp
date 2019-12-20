#include <iostream>
#include <chrono>
#include "HitoriSolverCore.h"

#define BFS
//#define GREEDY_BFS
//#define A_STAR
//#define HILL_CLIMBING
//#define S_HILL_CLIMBING
//#define SIMULATED_ANNEALING

int main() {
  HitoriSolverCore::HitoriSolver* solver;
  const char* path = "5x5.txt";
  auto heuristic = HitoriSolverCore::HitoriSolver::HeuristicFunction1;

#ifdef BFS
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result = solver->GreedyBfs(
      solver->InitialState(),
      [](const HitoriSolverCore::State&, int**) { return double(0); });
  solver->PrintState(result);
  delete solver;
#endif

#ifdef GREEDY_BFS
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result =
      solver->GreedyBfs(solver->InitialState(), heuristic);
  solver->PrintState(result);
  delete solver;
#endif

#ifdef A_STAR
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result =
      solver->AStar(solver->InitialState(), heuristic);
  solver->PrintState(result);
  delete solver;
#endif

#ifdef HILL_CLIMBING
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result = solver->SteepestAscentHillClimbing(
      solver->InitialState(), heuristic);
  solver->PrintState(result);
  delete solver;
#endif

#ifdef S_HILL_CLIMBING
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result = solver->StochasticHillClimbing(
      solver->InitialState(), heuristic);
  solver->PrintState(result);
  delete solver;
#endif

#ifdef SIMULATED_ANNEALING
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  HitoriSolverCore::State result = solver->SimulatedAnnealing(
      solver->InitialState(), heuristic,
      HitoriSolverCore::HitoriSolver::SCurveSchedule);
  solver->PrintState(result);
  delete solver;
#endif

  return 0;
}