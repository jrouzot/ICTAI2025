#include "ortools/linear_solver/linear_solver.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"

#include "FlowModelRelaxed.h"

namespace optimization {

    FlowModelRelaxed::FlowModelRelaxed(dataflow::Instance &instance, logger::Logger &logger) :
            mInstance(instance), mLogger(logger) {}

    dataflow::TransferSolution FlowModelRelaxed::Solve(
            std::vector<std::vector<double>> &handoverConstraints,
            std::vector<std::vector<double>> &handoverFixed) {

        mLogger.Debug("Building Relaxed Data Transfer Model...");

        std::string solverName = "SCIP";

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
         * The transfer rate variables keeps tracks of the new transfer rate at each window, for each buffer
         */
        std::vector< std::vector <operations_research::MPVariable*>> transferRatesVars({});
        transferRatesVars.resize(mInstance.getNumInstruments());
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            transferRatesVars[b].resize(mInstance.getNumDownlinks());
            for(auto w(0); w < mInstance.getNumDownlinks(); ++w) {
                double timeLag = mInstance.getDownlink(w).end - mInstance.getDownlink(w).start;
                assert(timeLag > 0.0);
                auto name = "transferRate_" + std::to_string(b) + "_" + std::to_string(w);
                // Get the maximum value at this timepoint (if outside downlink window -> 0, otherwise downlink window rate)
                double value = mInstance.getDownlink(w).rate;
                operations_research::MPVariable* const transferRateVar = solver.MakeNumVar(
                        0.0,
                        value * timeLag,
                        name);
                transferRatesVars[b][w] = transferRateVar;
            }
        }

        /*
         * The memory variables keep track of the memory state at each window for each buffer
         */
        std::vector< std::vector <operations_research::MPVariable*>> memoryVariables({});
        memoryVariables.resize(mInstance.getNumInstruments());
        for(auto b(0); b < mInstance.getNumInstruments(); ++b) {
            memoryVariables[b].resize(mInstance.getNumDownlinks());
            for (auto w(0); w < mInstance.getNumDownlinks(); ++w) {
                auto name = "memory_" + std::to_string(b) + "_" + std::to_string(w);
                if(w == 0) {    // Add the initial memory state
                    operations_research::MPVariable *const memoryVar = solver.MakeNumVar(
                            mInstance.getInstrument(b).getDataOffset(w) + mInstance.getInstrument(b).getInitialMemoryState(),
                            mInstance.getInstrument(b).getDataOffset(w) + mInstance.getInstrument(b).getInitialMemoryState(),
                            name);
                    memoryVariables[b][w] = memoryVar;
                } else if(handoverFixed[w][b] != -1) {
                    // When a fixed memory state is specified, fix the corresponding memory var
                    operations_research::MPVariable *const memoryVar = solver.MakeNumVar(
                            handoverFixed[w][b],
                            handoverFixed[w][b],
                            name);
                    memoryVariables[b][w] = memoryVar;
                } else {
                    // The upper bound for the memory var is either the capacity or a given handover
                    operations_research::MPVariable *const memoryVar = solver.MakeNumVar(
                            mInstance.getInstrument(b).getDataOffset(w),
                            std::min(mInstance.getInstrument(b).getCapacity(), handoverConstraints[w-1][b] + mInstance.getInstrument(b).getDataOffset(w)),
                            name);
                    memoryVariables[b][w] = memoryVar;
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
         * Bandwidth share
         */
        for(auto w{0}; w < mInstance.getNumDownlinks(); ++w) {
            double bandwidth = mInstance.getDownlink(w).rate * (mInstance.getDownlink(w).end - mInstance.getDownlink(w).start);
            auto name = "bandwidthShareCt_" + std::to_string(w);
            // sum(i) tr(i,j) <= bandwidth --- for all j
            operations_research::MPConstraint* const bandwidthCt = solver.MakeRowConstraint(0.0, bandwidth, name);
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                bandwidthCt->SetCoefficient(transferRatesVars[b][w], 1.0);
            }
        }

        /*
         * Memory consistency
         */
        for(auto w{1}; w < mInstance.getNumDownlinks(); ++w) {
            for (auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                auto name = "memConsistencyCt_" + std::to_string(b) + "_" + std::to_string(w);
                double totalFill = mInstance.getInstrument(b).getDataToDumpAt(mInstance.getDownlink(w-1).start);
                // mem(i,j) = mem(i,j-1) + [f(i,j-1) - tr(i-1,j)] * timeLag + data offset (during non visibility)
                operations_research::MPConstraint *const memConsistencyCt = solver.MakeRowConstraint(
                        totalFill + mInstance.getInstrument(b).getDataOffset(w),
                        totalFill + mInstance.getInstrument(b).getDataOffset(w), name);
                memConsistencyCt->SetCoefficient(memoryVariables[b][w], 1.0);
                memConsistencyCt->SetCoefficient(memoryVariables[b][w-1], -1.0);
                memConsistencyCt->SetCoefficient(transferRatesVars[b][w-1], 1.0);
            }
        }

        /*
         * Handover consistency (the sum of handovers is constant so all available bandwidth is used)
         */
        for(auto w{1}; w < mInstance.getNumDownlinks(); ++w) {
            auto name = "handoverConsistencyCt_" + std::to_string(w);
            double sumNonVisibility = 0.0;
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                sumNonVisibility += mInstance.getInstrument(b).getDataOffset(w);
            }
            operations_research::MPConstraint *const handoverCt = solver.MakeRowConstraint(
                mInstance.getSumHandover(w-1) + sumNonVisibility,
                mInstance.getSumHandover(w-1) + sumNonVisibility,
                name);
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                handoverCt->SetCoefficient(memoryVariables[b][w], 1.0);
            }
        }

        /*
         * Auxiliary variable for the objective
         */
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            for(auto w{0}; w < mInstance.getNumDownlinks(); ++w) {
                auto name = "alpha_" + std::to_string(b) + "_" + std::to_string(w);
                // alpha >= mem(i,j)
                operations_research::MPConstraint *const alphaCt = solver.MakeRowConstraint(
                        0.0, solver.infinity(), name);
                alphaCt->SetCoefficient(alpha, mInstance.getInstrument(b).getCapacity());
                alphaCt->SetCoefficient(memoryVariables[b][w], -1.0);
            }
        }

        std::string modelString;
        solver.ExportModelAsLpFormat(false, &modelString);
        std::cout << modelString << std::endl;

        /*************************
         *                       *
         * -       SOLVE       - *
         *                       *
         *************************/

        // Save current time
        absl::Time start_time = absl::Now();
        // Adding a 10 seconds time limit
        solver.SetTimeLimit(absl::Seconds(10));
        // Solve
        const operations_research::MPSolver::ResultStatus result_status = solver.Solve();
        // Calculate the solve time
        absl::Duration duration = absl::Now() - start_time;

        // Check that the problem has an optimal solution.
        if (result_status == operations_research::MPSolver::INFEASIBLE) {
            mLogger.Debug("The problem does not have a solution!");
            return dataflow::TransferSolution();
        }

        double rmax = alpha->solution_value();
        std::vector<std::vector<double>> handover;
        handover.resize(mInstance.getNumDownlinks());

        for(auto w{0}; w < mInstance.getNumDownlinks()-1; ++w) {
            for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
                handover[w].resize(mInstance.getNumInstruments());
                handover[w][b] = memoryVariables[b][w+1]->solution_value() - mInstance.getInstrument(b).getDataOffset(w+1);
            }
        }
        // Last window is a dummy window with no production, no transfer
        handover[mInstance.getNumDownlinks()-1].resize(mInstance.getNumInstruments());
        for(auto b{0}; b < mInstance.getNumInstruments(); ++b) {
            handover[mInstance.getNumDownlinks()-1][b] =
                    handover[mInstance.getNumDownlinks()-2][b] +
                    mInstance.getInstrument(b).getDataOffset(mInstance.getNumDownlinks()-1);
        }

        dataflow::TransferSolution solution(rmax, handover, absl::ToDoubleMilliseconds(duration));
        return solution;
    }
}

