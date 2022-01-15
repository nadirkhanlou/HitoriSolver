#include <iostream>
#include <chrono>
#include "HitoriSolverCore.h"

//#define BFS
//#define GREEDY_BFS
//#define A_STAR
//#define HILL_CLIMBING
//#define S_HILL_CLIMBING
//#define K_START_S_HILL_CLIMBING
//#define SIMULATED_ANNEALING

int main() {
  HitoriSolverCore::HitoriSolver* solver;
  const char* path = "12x12.txt";
  auto heuristic = HitoriSolverCore::HitoriSolver::HeuristicFunction1;
  HitoriSolverCore::State result;
  std::chrono::high_resolution_clock::time_point t1, t2;
  int durationPast;

#ifdef BFS
  std::cout << "========================================\n";
  std::cout << "BFS\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->GreedyBfs(
      solver->InitialState(),
      [](const HitoriSolverCore::State&, int**) { return double(0); });
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

#ifdef GREEDY_BFS
  std::cout << "========================================\n";
  std::cout << "Greedy BFS\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->GreedyBfs(solver->InitialState(), heuristic);
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

#ifdef A_STAR
  std::cout << "========================================\n";
  std::cout << "A*\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->AStar(solver->InitialState(), heuristic);
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

#ifdef HILL_CLIMBING
  std::cout << "========================================\n";
  std::cout << "Steepest Ascent Hill Climbing\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result =
      solver->SteepestAscentHillClimbing(solver->InitialState(), heuristic);
  if (solver->IsGoal(result))
    std::cout << "Goal\n";
  else
    std::cout << "Not Goal\n";
  solver->PrintState(result);
  t2 =
      std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

#ifdef K_START_S_HILL_CLIMBING
  std::cout << "========================================\n";
  std::cout << "K Start Stochastic Hill Climbing\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->KStartStochasticHillClimbing(solver->InitialState(), heuristic);
  if (solver->IsGoal(result))
    std::cout << "Goal\n";
  else
    std::cout << "Not Goal\n";
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

  #ifdef S_HILL_CLIMBING
  std::cout << "========================================\n";
  std::cout << "Stochastic Hill Climbing\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->StochasticHillClimbing(solver->InitialState(), heuristic);
  if (solver->IsGoal(result))
    std::cout << "Goal\n";
  else
    std::cout << "Not Goal\n";
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

#ifdef SIMULATED_ANNEALING
  std::cout << "========================================\n";
  std::cout << "Simulated Annealing\n";
  t1 = std::chrono::high_resolution_clock::now();
  solver = new HitoriSolverCore::HitoriSolver(path);
  solver->PreProcess();
  result = solver->SimulatedAnnealing(
      solver->InitialState(), heuristic,
      HitoriSolverCore::HitoriSolver::SCurveSchedule);
  if (solver->IsGoal(result))
    std::cout << "Goal\n";
  else
    std::cout << "Not Goal\n";
  solver->PrintState(result);
  t2 = std::chrono::high_resolution_clock::now();
  durationPast =
      (int)(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
                .count() *
            1000);
  std::cout << "Duration: " << durationPast << "ms\n";
  std::cout << "========================================\n";
  delete solver;
#endif

  return 0;
}
