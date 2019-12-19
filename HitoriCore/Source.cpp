#include <iostream>
#include "HitoriSolverCore.h"


#define BFS
//#define GREEDY_BFS
//#define A_STAR
//#define HILL_CLIMBING
//#define S_HILL_CLIMBING

double h(const HitoriSolverCore::State&, int**) { return 0; }

int main() {
  HitoriSolverCore::HitoriSolver solver =
      HitoriSolverCore::HitoriSolver("15x15.txt");

  solver.PreProccess();


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
  HitoriSolverCore::State result =
      solver.SteepestAscentHillClimbing(solver.InitialState(),
                   HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

#ifdef S_HILL_CLIMBING
  HitoriSolverCore::State result =
      solver.StochasticHillClimbing(
          solver.InitialState(),
          HitoriSolverCore::HitoriSolver::HeuristicFunction1);
#endif

  solver.PrintState(result);

  return 0;
}