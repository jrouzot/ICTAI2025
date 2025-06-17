//
// Created by jrouzot on 10/07/24.
//

#include "ortools/linear_solver/linear_solver.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"

#include "../Solution/Solution.h"
#include "FlowDataTransferModel.h"

namespace optimization {

    FlowDataTransferModel::FlowDataTransferModel(dataflow::Instance &instance, logger::Logger &logger) :
    mInstance(instance), mLogger(logger) {}

    /*
     * Build the LP model to solve the flow problem from @windowStart
     */
    dataflow::TransferSolution FlowDataTransferModel::Solve(
            int windowStart,
            std::vector<double> &memoryState) {

        // Start time must be 0 or the end of a previous downlink window
        assert(windowStart < mInstance.getNumDownlinks());

        mLogger.Debug("Building Data Transfer Model...");

        std::string solverName = "SCIP";
        // std::string solverName = "GLOP";

        // Create the solver instance
        std::unique_ptr <operations_research::MPSolver> solverPtr(operations_research::MPSolver::CreateSolver(solverName));

        // Check if solver is available
        if (solverPtr == nullptr) {
            mLogger.Error("Solver is unavailable!");
            exit(1);
        } else {
            mLogger.Debug("Solver successfully initialized!");
        }

        operations_research::MPSolver &solver = *solverPtr;

        // Setup the precision parameters
        std::string scipParameters = "numerics/feastol = 1e-4\nlimits/gap = 1e-6";
        solver.SetSolverSpecificParametersAsString(scipParameters);

        /***********************
         *                     *
         * -    VARIABLES    - *
         *                     *
         ***********************/

        /*
         * The transfer rate variables keeps tracks of the new transfer rate at each inflexion point, for each buffer
         */
        std::vector< std::vector <operations_research::MPVariable*>> transferRatesVars({});
        transferRatesVars.resize(mInstance.getNumInstruments());
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            transferRatesVars[b].resize(mInstance.getNumTimepoints());
            for(auto t(0); t < mInstance.getNumTimepoints(); ++t) {
                auto name = "transferRate_" + std::to_string(b) + "_" + std::to_string(t);
                // Get the maximum value at this timepoint (if outside downlink window -> 0, otherwise downlink window rate)
                double value = mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]);
                if((value > 0 or t == mInstance.getNumTimepoints()-2)
                and mInstance.getTimepoints()[t] >= mInstance.getDownlink(windowStart).start) {
                    operations_research::MPVariable* const transferRateVar = solver.MakeNumVar(
                            0.0,
                            value,
                            name);
                    transferRatesVars[b][t] = transferRateVar;
                } else {
                    transferRatesVars[b][t] = nullptr;
                }
            }
        }

        /*
         * The memory variables keep track of the memory state at each inflexion point for each buffer
         */
        std::vector< std::vector <operations_research::MPVariable*>> memoryVariables({});
        memoryVariables.resize(mInstance.getNumInstruments());
        for(auto b(0); b < mInstance.getNumInstruments(); ++b) {
            memoryVariables[b].resize(mInstance.getNumTimepoints());
            for (auto t(0); t < mInstance.getNumTimepoints(); ++t) {
                auto name = "memory_" + std::to_string(b) + "_" + std::to_string(t);
                double value = mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]);
                if(mInstance.getTimepoints()[t] >= mInstance.getDownlink(windowStart).start
                and (value > 0 or t == mInstance.getNumTimepoints()-2
                or std::find(mInstance.GetHandoverIndexes().begin(), mInstance.GetHandoverIndexes().end(), t) != mInstance.GetHandoverIndexes().end())) {
                    operations_research::MPVariable *const memoryVar = solver.MakeNumVar(
                            0.0,
                            mInstance.getInstrument(b).getCapacity(),
                            name);
                    memoryVariables[b][t] = memoryVar;
                } else {
                    memoryVariables[b][t] = nullptr;
                }
            }
        }

        /*
         * Auxiliary variable for the rmax
         */
        operations_research::MPVariable *alpha;
        alpha = solver.MakeNumVar(0.0, 1.0, "alpha");


        /*************************
         *                       *
         * -     OBJECTIVE     - *
         *                       *
         *************************/

        operations_research::MPObjective* const objective = solver.MutableObjective();
        objective->SetCoefficient(alpha, 1.0);
        objective->SetMinimization();



        /*************************
         *                       *
         * -    CONSTRAINTS    - *
         *                       *
         *************************/

        /*
         * Memory start
         */
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            auto name = "memoryStartCt_" + std::to_string(b);
            double memoryStart = memoryState[b];
            operations_research::MPConstraint* const memStartConstraint = solver.MakeRowConstraint(memoryStart, memoryStart, name);
            memStartConstraint->SetCoefficient(memoryVariables[b][windowStart], 1.0);
        }

        /*
         * Bandwidth share
         */
        for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
            if(transferRatesVars[0][t] != nullptr) {
                double bandwidth = mInstance.getDumpRateAtTime(mInstance.getTimepoints()[t]);
                auto name = "bandwidthShareCt_" + std::to_string(t);
                // sum(i) tr(i,j) = bandwidth --- for all j
                operations_research::MPConstraint* const bandwidthCt = solver.MakeRowConstraint(0.0, bandwidth, name);
                for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                    bandwidthCt->SetCoefficient(transferRatesVars[b][t], 1.0);
                }
            }
        }

        /*
         * Memory consistency (same window)
         */
        for(auto t{1}; t < mInstance.getNumTimepoints(); ++t) {
            if(transferRatesVars[0][t-1] != nullptr) {
                double timeLag = mInstance.getTimepoints()[t] - mInstance.getTimepoints()[t - 1];
                assert(timeLag > 0.0);
                for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                    auto name = "memConsistencyCt_" + std::to_string(b) + "_" + std::to_string(t - 1);
                    double fillRate = mInstance.getFillRateAt(b, mInstance.getTimepoints()[t - 1]);
                    // mem(i,j) = mem(i,j-1) + tr(i-1,j) * timeLag
                    operations_research::MPConstraint *const memConsistencyCt = solver.MakeRowConstraint(
                            fillRate * timeLag, fillRate * timeLag, name);
                    memConsistencyCt->SetCoefficient(memoryVariables[b][t], 1.0);
                    memConsistencyCt->SetCoefficient(memoryVariables[b][t - 1], -1.0);
                    memConsistencyCt->SetCoefficient(transferRatesVars[b][t - 1], timeLag);
                }
            }
        }

        /*
         * Memory consistency (non visibility)
         */
        for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            int prevIndex;
            int windowIndex;
            for (int w{1}; w < mInstance.GetStartWindowIndexes().size(); ++w) {
                if(w > windowStart) {
                    prevIndex = mInstance.GetHandoverIndexes()[w-1];
                    windowIndex = mInstance.GetStartWindowIndexes()[w];
                    auto name = "windowConsistencyCt_" + std::to_string(b) + "_" + std::to_string(w);
                    double sum = mInstance.getInstrument(b).getDataOffset(w);
                    operations_research::MPConstraint *const windowConsistencyCt = solver.MakeRowConstraint(
                            sum, sum, name);
                    windowConsistencyCt->SetCoefficient(memoryVariables[b][windowIndex], 1.0);
                    windowConsistencyCt->SetCoefficient(memoryVariables[b][prevIndex], -1.0);
                }
            }
        }

        /*
         * Handover consistency (the sum of handovers is constant so all the bandwidth is used)
         */
        int currentWindowIndex = windowStart;
        for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
            // Add the constraints only for the handover timepoints
            if(mInstance.GetHandoverIndexes()[currentWindowIndex] == t) {
                auto name = "handoverConsistencyCt_" + std::to_string(currentWindowIndex);
                operations_research::MPConstraint *const handoverCt = solver.MakeRowConstraint(
                        mInstance.getSumHandover(currentWindowIndex),
                        mInstance.getSumHandover(currentWindowIndex),
                        name);
                for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                    handoverCt->SetCoefficient(memoryVariables[b][t], 1.0);
                }
                ++currentWindowIndex;
            }
        }


        /*
         * Auxiliary variable for the objective
         */
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
                if(memoryVariables[0][t] != nullptr) {
                    auto name = "alpha_" + std::to_string(b) + "_" + std::to_string(t);
                    // alpha >= mem(i,j)
                    operations_research::MPConstraint *const alphaCt = solver.MakeRowConstraint(
                            0.0, solver.infinity(), name);
                    alphaCt->SetCoefficient(alpha, mInstance.getInstrument(b).getCapacity());
                    alphaCt->SetCoefficient(memoryVariables[b][t], -1.0);
                }
            }
        }


        /*************************
         *                       *
         * -    PRINT MODEL    - *
         *                       *
         *************************/

/*        std::string modelString;
        solver.ExportModelAsLpFormat(false, &modelString);
        std::cout << modelString << std::endl;*/

        /*************************
         *                       *
         * -       SOLVE       - *
         *                       *
         *************************/

        // Save current time
        absl::Time start_time = absl::Now();
        // Solve
        const operations_research::MPSolver::ResultStatus result_status = solver.Solve();
        // Calculate the solve time
        absl::Duration duration = absl::Now() - start_time;

        // Check that the problem has an optimal solution.
        if (result_status != operations_research::MPSolver::OPTIMAL) {
            mLogger.Debug("The problem does not have a solution!");
            return dataflow::TransferSolution();
        }

        std::vector<std::vector<double>> handover;
        handover.resize(mInstance.getNumDownlinks());
        currentWindowIndex = windowStart;
        double rmax = 0.0;
        for(auto t{0}; t < mInstance.getNumTimepoints(); ++t) {
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                if(memoryVariables[b][t]->solution_value() / mInstance.getInstrument(b).getCapacity() > rmax) {
                    rmax = memoryVariables[b][t]->solution_value() / mInstance.getInstrument(b).getCapacity();
                }
            }
            if(currentWindowIndex >= mInstance.getNumDownlinks()) {
                break;
            }
            if(mInstance.GetHandoverIndexes()[currentWindowIndex] == t) {
                handover[currentWindowIndex].resize(mInstance.getNumInstruments());
                for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                    handover[currentWindowIndex][b] = memoryVariables[b][t]->solution_value();
                }
                ++currentWindowIndex;
            }
        }

        dataflow::TransferSolution solution(rmax, handover, absl::ToDoubleMilliseconds(duration));
        return solution;
    }

} // optimization