//
// Created by Julien Rouzot on 20/05/24.
//

#include <iostream>
#include <memory>
#include <string>

#include "ortools/linear_solver/linear_solver.h"

#include "FlowModel.h"
#include "../Logger/Logger.h"

namespace optimization {

    FlowModel::FlowModel(dataflow::Instance &instance, logger::Logger &logger) :
    mInstance(instance), mLogger(logger) {}

    dataflow::FillRateSolution FlowModel::Solve() {
        dataflow::FillRateSolution solution({});
        std::unique_ptr <operations_research::MPSolver> solverPtr(operations_research::MPSolver::CreateSolver("SCIP"));// Create the solver instance
        operations_research::MPSolver &solver = *solverPtr;                                                            // Get the solver object
        if (solverPtr == nullptr) {                                                                                    // Check if solver is available
            mLogger.Error("SCIP solver is unavailable!");
        } else {
            mLogger.Debug("SCIP solver successfully initialized!");
        }
        const double infinity = solver.infinity();                                                                      // Largest double

        // Create a fill rate variable for each instrument and for each event
        std::vector< std::vector <operations_research::MPVariable*>> fillRateVars;
        fillRateVars.resize(mInstance.getNumInstruments());
        for(int i(0); i < mInstance.getNumInstruments(); ++i) {
            fillRateVars[i].resize(mInstance.getNumTimepoints());
            for(int j(0); j < mInstance.getNumTimepoints(); ++j) {
                auto name = "fillRate_{}_{}" + std::to_string(i) + std::to_string(j);
                operations_research::MPVariable* const fillRateVar = solver.MakeNumVar(
                        mInstance.getInstrument(i).mMinRate,
                        mInstance.getInstrument(i).mMaxRate,
                        name);
                fillRateVars[i][j] = fillRateVar;
            }
        }

        operations_research::MPVariable* const x = solver.MakeNumVar(0.0, infinity, "x");
        operations_research::MPVariable* const y = solver.MakeNumVar(0.0, infinity, "y");

        mLogger.Debug("Number of variables: ", solver.NumVariables());

        operations_research::MPConstraint* const c0 = solver.MakeRowConstraint(-infinity, 14.0);
        c0->SetCoefficient(x, 1.0);
        c0->SetCoefficient(y, 2.0);

        operations_research::MPConstraint* const c1 = solver.MakeRowConstraint(0.0, infinity);
        c1->SetCoefficient(x, 3.0);
        c1->SetCoefficient(y, -1.0);

        operations_research::MPConstraint* const c2 = solver.MakeRowConstraint(-infinity, 2.0);
        c2->SetCoefficient(x, 1.0);
        c2->SetCoefficient(y, -1.0);

        mLogger.Debug("Number of constraints: ", solver.NumConstraints());

        operations_research::MPObjective* const objective = solver.MutableObjective();
        objective->SetCoefficient(x, 3.0);
        objective->SetCoefficient(y, 4.0);
        objective->SetMaximization();

        const operations_research::MPSolver::ResultStatus resultStatus = solver.Solve();
        if (resultStatus != operations_research::MPSolver::OPTIMAL) {                                                  // Check that the problem has an optimal solution.
            mLogger.Warn("No solution found!");
        }

        mLogger.Debug("Optimal objective value = ", objective->Value());
        mLogger.Debug(x->name(), " = ", x->solution_value());
        mLogger.Debug(y->name(), " = ", y->solution_value());
        mLogger.Debug("Solution:");

        return solution;
    }
}

