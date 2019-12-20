#include <iostream>
#include "HitoriSolverCore.h"

#define BFS
//#define GREEDY_BFS
//#define A_STAR
//#define HILL_CLIMBING
//#define S_HILL_CLIMBING
//#define SIMULATED_ANNEALING

double h(const HitoriSolverCore::State&, int**) { return 0; }

std::vector<double>* SCurveSchedule(int steps) {
  double t;
  double increment;
  t = 0.01;
  increment = (1 - t) / double(steps);
  std::vector<double>* temperatures = new std::vector<double>;
  while (t <= 1) {
    double temp = (1) / (1 + pow((t / (1 - t)), (3.0)));
    temperatures->push_back(temp);
    t += increment;
  }
  return temperatures;
}

int main() {
  HitoriSolverCore::HitoriSolver solver =
      HitoriSolverCore::HitoriSolver("20x20.txt");

  solver.PreProcess();

#ifdef BFS
  HitoriSolverCore::State result = solver.GreedyBfs(solver.InitialState(), h);
#endif

#ifdef GREEDY_BFS
  HitoriSolverCore::State result =
      solver.GreedyBfs(solver.InitialState(),
                       HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

#ifdef A_STAR
  HitoriSolverCore::State result =
      solver.AStar(solver.InitialState(),
                   HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

#ifdef HILL_CLIMBING
  HitoriSolverCore::State result = solver.SteepestAscentHillClimbing(
      solver.InitialState(),
      HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

#ifdef S_HILL_CLIMBING
  HitoriSolverCore::State result = solver.StochasticHillClimbing(
      solver.InitialState(),
      HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

#ifdef SIMULATED_ANNEALING
  HitoriSolverCore::State result = solver.SimulatedAnnealing(
      solver.InitialState(), HitoriSolverCore::HitoriSolver::HeuristicFunction2,
      SCurveSchedule);
#endif

  solver.PrintState(result);

  return 0;
}