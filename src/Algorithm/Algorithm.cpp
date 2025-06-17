#include "Algorithm.h"

namespace simulation {

    Algorithm::Algorithm(Instance &instance, logger::Logger &logger)
            : instance(instance), mLogger(logger), simulator(instance, logger) {

        rmax = 1;
        rmaxLowerBound = 0;

        noStop.resize(instance.getNumInstruments(), false);

        priorityTmp.resize(instance.getNumDownlinks());
        prioritySolution.resize(instance.getNumDownlinks());
        handoverSolution.resize(instance.getNumDownlinks());

        overflowingBuffers.reserve(instance.getNumInstruments());
        overflowingBuffersMinViolation.reserve(instance.getNumInstruments());

        buffers.reserve(instance.getNumInstruments());
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            buffers.emplace_back(b);
        }

        dump.resize(instance.getNumDownlinks());
        for(auto w{0}; w < instance.getNumDownlinks(); ++w) {
            dump[w] = instance.getDownlink(w).rate;
        }

        window.resize(instance.getNumDownlinks());

        target.resize(instance.getNumInstruments());
        swTarget.resize(instance.getNumInstruments());
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            target[b] = instance.getInstrument(b).getCapacity();
            swTarget[b] = instance.getInstrument(b).getCapacity();
        }

        stop.resize(instance.getNumDownlinks());
        for (int w{0}; w < instance.getNumDownlinks(); ++w) {
            stop[w].resize(instance.getNumInstruments());
            for (int b{0}; b < instance.getNumInstruments(); ++b) {
                stop[w][b] = instance.getDownlink(w).end;
            }
        }

        pointer.resize(instance.getNumDownlinks());
        pointer[0].resize(instance.getNumInstruments(), 0);
        usage.resize(instance.getNumDownlinks());
        usage[0].resize(instance.getNumInstruments());
        maxUsage.resize(instance.getNumDownlinks());
        maxUsage[0].resize(instance.getNumInstruments(), 0);
        handover.resize(instance.getNumDownlinks());
        handover[0].resize(instance.getNumInstruments(), 0);

        // Initialise usage for the first window to the initialMemoryState
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            usage[0][b] = instance.getInstrument(b).getInitialMemoryState();
        }
        // Initialise usage for other windows
        for (auto w{1}; w < instance.getNumDownlinks(); ++w) {
            usage[w].resize(instance.getNumInstruments(), 0);
            maxUsage[w].resize(instance.getNumInstruments(), 0);
            pointer[w].resize(instance.getNumInstruments(), 0); // There is one pointer to the downlink
            handover[w].resize(instance.getNumInstruments(), 0);
        }

        dummyPriorityWindow = {};
        dummyPriorityWindow.reserve(instance.getNumInstruments());
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            dummyPriorityWindow.push_back({});
            dummyPriorityWindow[b].reserve(instance.getNumInstruments());
        }
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            dummyPriorityWindow[0].emplace_back(b);
        }

        dummyHandovers.resize(instance.getNumDownlinks());
        for(auto w{0}; w < instance.getNumDownlinks(); ++w) {
            dummyHandovers[w].resize(instance.getNumInstruments());
            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                dummyHandovers[w][b] = instance.getInstrument(b).getCapacity();
            }
        }

        flatPriority.resize(instance.getNumInstruments(), 0);

        prioritySolution.resize(instance.getNumDownlinks());
        for (auto w{0}; w < instance.getNumDownlinks(); ++w) {
            prioritySolution[w] = dummyPriorityWindow;
        }

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            for (auto e{instance.getInstrument(b).getEvents().begin()}; e < instance.getInstrument(b).getEvents().end(); ++e) {
                events.push_back({e->time, b});
            }
        }
        sort(events.begin(), events.end());

        int currentWindowIndex = 0;
        auto lastEvent = events.begin();
        for(auto e{events.begin()}; e < events.end(); ++e) {
            if(currentWindowIndex < instance.getNumDownlinks() and e->time >= instance.getDownlink(currentWindowIndex).start) {
                window[currentWindowIndex] = e;
                lastEvent = e;
                ++currentWindowIndex;
            }
        }
        while(currentWindowIndex < instance.getNumDownlinks()) {
            window[currentWindowIndex] = lastEvent;
            ++ currentWindowIndex;
        }
        window.emplace_back(events.end());

        // Initialise sumHandover
        sumHandover.resize(instance.getNumDownlinks());

        // Compute an upper bound for the current downlink window
        rmaxUpperBound = std::numeric_limits<double>::max();
        playSolution(prioritySolution);
        instance.setSumHandover(sumHandover);

        // Compute a upper bound for each handover
        computeHandoverUb();

        // Play solution with relaxed bandwidth to compute LB for handovers
        handoverLb.resize(instance.getNumDownlinks());
        for(int w{0}; w < instance.getNumDownlinks(); ++w) {
            handoverLb[w].resize(instance.getNumInstruments());
        }
        playSolutionFullBandwidth();
    }

    /*
     * Set the target
     */
    void Algorithm::setTargetRmax(double newRmax) {
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            target[b] = instance.getInstrument(b).getCapacity() * newRmax;
        }
    }

    /*
     * Set the swTarget
     */
    void Algorithm::setSwTargetRmax(double newRmax) {
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            swTarget[b] = instance.getInstrument(b).getCapacity() * newRmax;
        }
    }

    /*
     * Play a downlink without stop.
     * This is useful to compute the global handover at each downlink window
     */
    void Algorithm::playDownlinkWithoutStop(const int i, std::vector<std::vector<int>> &priority) {

        simulator.initialise(pointer[i], usage[i]);

        simulator.run(
                window[i], window[i+1], dump[i],
                instance.getDownlink(i).start, instance.getDownlink(i).end,
                priority);

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            // Update the pointers for the next window
            if(i+1 < instance.getNumDownlinks()) {
                // Update the pointers for the next window
                pointer[i+1][b] = simulator.getPointer(b)+1;
                // Update the usage for the next window
                usage[i+1][b] = simulator.getUsage(b);
            }
        }

        sumHandover[i] = 0;
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            // Update max usages for the current window
            maxUsage[i][b] = simulator.getMaxUsage(b);
            if (maxUsage[i][b]/instance.getInstrument(b).getCapacity() > rmaxUpperBound) {
                rmaxUpperBound = maxUsage[i][b]/instance.getInstrument(b).getCapacity();
            }
            sumHandover[i] += simulator.getUsage(b);
        }
    }

    /*
     * Play a single downlink window @i with a given priority and target handovers but only check the violations
     */
    double Algorithm::playDownlinkViolation(const int i, std::vector<std::vector<int>> &priority,
                                            std::vector<double> &handoverTargets, std::vector<double> &handoverConstraints) {

        simulator.initialise(pointer[i], usage[i]);

        simulator.runWithEarliestStopViolation(
                window[i], window[i+1], dump[i],
                instance.getDownlink(i).start, instance.getDownlink(i).end,
                priority, stop[i], target, handoverTargets);

        double maxViolation = -MAXFLOAT;
        overflowingBuffers = {};

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            // Update max usages for the current window
            maxUsage[i][b] = simulator.getMaxUsage(b);
            if(maxUsage[i][b] > target[b]) {
                overflowingBuffers.emplace_back(b);
                if(maxUsage[i][b] - target[b] > maxViolation) {
                    maxViolation = maxUsage[i][b] - target[b];
                }
            } else if(simulator.getUsage(b) > handoverConstraints[b]) {
                overflowingBuffers.emplace_back(b);
                // Check handover CONSTRAINTS overflow (the target are not considered as an overflow if exceeded
                if(simulator.getUsage(b) - handoverConstraints[b] > maxViolation) {
                    maxViolation = simulator.getUsage(b) - handoverConstraints[b];
                }
            }
            if(i+1 < instance.getNumDownlinks()) {
                // Update the pointers for the next window
                pointer[i+1][b] = simulator.getPointer(b)+1;
                // Update the usage for the next window
                usage[i+1][b] = simulator.getUsage(b);
            }
        }

        assert(maxViolation > 0);

        return maxViolation;
    }


    /*
     * Play a solution for the whole instance
     */
    void Algorithm::playSolution(std::vector<std::vector<std::vector<int>>> &priority) {

        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            usage[0][b] = instance.getInstrument(b).getInitialMemoryState() + instance.getInstrument(b).getDataOffset(0);
        }

        for(auto i{0}; i < instance.getNumDownlinks(); ++i) {
            playDownlinkWithoutStop(i, priority[i]);
            // Add the non visibility data offset for next window
            if(i < instance.getNumDownlinks() - 1) {
                for(int b{0}; b < instance.getNumInstruments(); ++b) {
                    usage[i+1][b] += instance.getInstrument(b).getDataOffset(i+1);
                }
            }
        }
    }

    /*
     * Check if current best solution and handovers match given @rmax
     */
    bool Algorithm::checkSolution(double currentRmax) {

        setTargetRmax(currentRmax);
        rmax = 0.0;     // We reset rmax

        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            usage[0][b] = instance.getInstrument(b).getInitialMemoryState() + instance.getInstrument(b).getDataOffset(0);
        }

        for(auto i{0}; i < instance.getNumDownlinks(); ++i) {
            playDownlinkEarliestStop(i, prioritySolution[i], handoverSolution[i], handoverUb[i]);
            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(maxUsage[i][b] >= target[b] + 0.01) {
                    mLogger.Error("Exceeded target rmax at ", i, " for buffer ", b);
                    return false;
                }
                if(maxUsage[i][b] / instance.getInstrument(b).getCapacity() > rmax + precision) {
                    rmax = maxUsage[i][b] / instance.getInstrument(b).getCapacity();
                }
            }

            // Add the non visibility data offset for next window
            if(i < instance.getNumDownlinks() - 1) {
                for(int b{0}; b < instance.getNumInstruments(); ++b) {
                    usage[i+1][b] = handoverSolution[i][b] + instance.getInstrument(b).getDataOffset(i+1);
                }
            }
        }
        return true;
    }


    /*
     * Play a downlink window with earliest stops according to a priority and handovers
     */
    void Algorithm::playDownlinkEarliestStop(const int i, std::vector<std::vector<int>> &priority,
                                             std::vector<double> &handoverTargets, std::vector<double> &handoverConstraints) {

        simulator.initialise(pointer[i], usage[i]);
/*
        simulator.setNoStop(noStop);
*/

        int res = simulator.runWithEarliestStop(
                window[i], window[i+1], dump[i],
                instance.getDownlink(i).start, instance.getDownlink(i).end,
                priority, stop[i], target, handoverTargets, handoverConstraints);

        if(res != -1)
            mLogger.Debug(res);

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            // Update max usages for the current window
            maxUsage[i][b] = simulator.getMaxUsage(b);
            assert(maxUsage[i][b] <= target[b] + precision);
            if (i + 1 < instance.getNumDownlinks()) {
                // Update the pointers for the next window
                pointer[i + 1][b] = simulator.getPointer(b) + 1;
                // Update the usage for the next window
                usage[i + 1][b] = simulator.getUsage(b);
            }
        }
    }


    /*
     * Play a relaxed downlink (full bandwidth is allocated to every buffer)
     */
    void Algorithm::playDownlinkFullBandwidth(const int i) {

        simulator.initialise(pointer[i], usage[i]);

        simulator.runFullBandwidth(
                window[i], window[i+1], dump[i],
                instance.getDownlink(i).start, instance.getDownlink(i).end);

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            // Update max usages for the current window
            maxUsage[i][b] = simulator.getMaxUsage(b);
            if(i+1 < instance.getNumDownlinks()) {
                // Update the pointers for the next window
                pointer[i+1][b] = simulator.getPointer(b)+1;
                // Update the usage for the next window
                usage[i+1][b] = simulator.getUsage(b);
                // Update the handoverLB
                handoverLb[i][b] = simulator.getUsage(b);
            }
        }
    }


    /*
     * Play a relaxed solution (full bandwidth is allocated to every buffer)
     */
    void Algorithm::playSolutionFullBandwidth() {

        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            usage[0][b] = instance.getInstrument(b).getInitialMemoryState() + instance.getInstrument(b).getDataOffset(0);
        }

        for(auto i{0}; i < instance.getNumDownlinks(); ++i) {
            playDownlinkFullBandwidth(i);
            // Add the non visibility data offset for next window
            if(i < instance.getNumDownlinks() - 1) {
                for (int b{0}; b < instance.getNumInstruments(); ++b) {
                    usage[i + 1][b] += instance.getInstrument(b).getDataOffset(i + 1);
                }
            }
        }
    }

    /*
     * Play the simulation without data dumping. This is useful to compute handover upper bounds
     */
    void Algorithm::computeHandoverUb() {
        double memoryState;
        handoverUb.resize(instance.getNumDownlinks());
        for(auto w{0}; w < instance.getNumDownlinks(); ++w) {
            handoverUb[w].resize(instance.getNumInstruments());
        }
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            memoryState = instance.getInstrument(b).getInitialMemoryState();
            for(auto w{0}; w < instance.getNumDownlinks(); ++w) {
                memoryState += instance.getInstrument(b).getDataOffset(w);
                int timepoint = instance.getTimepoints()[instance.GetStartWindowIndexes()[w]];
                memoryState += instance.getInstrument(b).getDataToDumpAt(timepoint);
                if(memoryState >= instance.getInstrument(b).getCapacity()) {
                    memoryState = instance.getInstrument(b).getCapacity();
                }
                handoverUb[w][b] = memoryState;
            }
        }
    }


    /*
     *  SingleWindow with interruptions.
     *  This algorithm seeks a priority assignment for a single window and a given margin
     */
    bool Algorithm::singleWindow(int windowIndex, double targetRmax, std::vector<double> &currentHandoverTargets, std::vector<double> &currentHandoverConstraints) {

        ++nbSingleWindow;
        // Reset priority for this window
        bestPriority = dummyPriorityWindow;
        // Keep tracks of the positions of the buffers in bestPriority
        flatPriority.resize(instance.getNumInstruments());
        std::fill(flatPriority.begin(), flatPriority.end(), 0);
        // Least violation
        singleWindowMinViolation = MAXFLOAT;
        // Indicate if we try to reach the targets
        bool relax{false};
        // current constraints
        std::vector<double> myHandoverConstraints = currentHandoverTargets;

        // Try to build the priority assignment
        while (true) {

            // Initialise the simulator
            simulator.initialise(pointer[windowIndex], usage[windowIndex]);

            int overflowIndex = simulator.runWithEarliestStop(
                    window[windowIndex], window[windowIndex+1], dump[windowIndex],
                    instance.getDownlink(windowIndex).start, instance.getDownlink(windowIndex).end,
                    bestPriority, stop[windowIndex], target, currentHandoverTargets, myHandoverConstraints);

            ++nbSimulation;

            // Raise the priority of the first buffer that overflows
            if (overflowIndex >= 0) {

                // Compute the violation
                if (relax) {
                    double violation = playDownlinkViolation(windowIndex, bestPriority, currentHandoverTargets, currentHandoverConstraints);
                    if(violation < singleWindowMinViolation) {
                        singleWindowMinViolation = violation;
                        overflowingBuffersMinViolation = overflowingBuffers;
                        assert(overflowingBuffersMinViolation.size() >= 1);
                    }
                }

                if (flatPriority[overflowIndex] == instance.getNumInstruments() - 1) {
                    // buffer already had the best priority
                    if(relax) {
                        return false;
                    } else {
                        relax = true;
                        // Reset priority for this window
                        bestPriority = dummyPriorityWindow;
                        // Keep tracks of the positions of the buffers in bestPriority
                        flatPriority.resize(instance.getNumInstruments());
                        std::fill(flatPriority.begin(), flatPriority.end(), 0);
                        // Least violation
                        singleWindowMinViolation = MAXFLOAT;
                        // Update handover constraints
                        myHandoverConstraints = currentHandoverConstraints;
                    }
                } else {
                    auto it = std::find(bestPriority[flatPriority[overflowIndex]].begin(),
                                        bestPriority[flatPriority[overflowIndex]].end(), overflowIndex);

                    // Remove buffer from current priority pool (swap it with first element for more efficiency, order doesnt matter)
                    std::iter_swap(it, bestPriority[flatPriority[overflowIndex]].end() - 1);
                    bestPriority[flatPriority[overflowIndex]].pop_back();
                    // Push buffer in next pool
                    bestPriority[flatPriority[overflowIndex] + 1].push_back(overflowIndex);
                }

                // Check for inconsistency
                if (bestPriority[flatPriority[overflowIndex]].empty()) {
                    if(relax) {
                        return false;
                    } else {
                        relax = true;
                        // Reset priority for this window
                        bestPriority = dummyPriorityWindow;
                        // Keep tracks of the positions of the buffers in bestPriority
                        flatPriority.resize(instance.getNumInstruments());
                        std::fill(flatPriority.begin(), flatPriority.end(), 0);
                        // Least violation
                        singleWindowMinViolation = MAXFLOAT;
                        // Update handover constraints
                        myHandoverConstraints = currentHandoverConstraints;
                    }
                } else {
                    flatPriority[overflowIndex]++;
                }

            } else {
                if(windowIndex+1 < instance.getNumDownlinks()) {
                    for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                        usage[windowIndex+1][b] = simulator.getUsage(b);
                        // Update the pointers for the next window
                        pointer[windowIndex+1][b] = simulator.getPointer(b) + 1;
                    }
                }
                for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                    maxUsage[windowIndex][b] = simulator.getMaxUsage(b);
                }
                return true;
            }
        }
    }


    /*
     * Run downlinkCount heuristic to find target handover and rmax
     */
    void Algorithm::runDownlinkCount(std::vector<double> &initialMemoryState, int nbRuns) {

        double start = std::clock();

        // Best rmax
        rmax = MAXFLOAT;
        // Best handover
        std::vector<std::vector<double>> bestHandover = handover;

        std::vector<std::vector<std::vector<int>>> bestPrioritySolution;
        std::vector<std::vector<std::vector<int>>> currentPrioritySolution;

        for(int i{0}; i < nbRuns; ++i) {
            // Reset current Priority
            currentPrioritySolution = {};
            currentPrioritySolution.reserve(instance.getNumDownlinks());

            // Reset rmaxUpperBound
            rmaxUpperBound = 0.0;

            DownlinkCount dcHeuristic(instance, rmax);
            usage[0] = initialMemoryState;

            for(int w{0}; w < instance.getNumDownlinks(); ++w) {
                bestPriority = dcHeuristic.run(w, usage[w]);
                currentPrioritySolution.push_back(bestPriority);
                playDownlinkWithoutStop(w, bestPriority);
                // Add the non visibility data offset for next window
                for(int b{0}; b < instance.getNumInstruments(); ++b) {
                    if(w < instance.getNumDownlinks() - 1) {
                        usage[w+1][b] += instance.getInstrument(b).getDataOffset(w+1);
                        handover[w][b] = usage[w+1][b];
                    }
                }
            }

            if(rmaxUpperBound < rmax) {
                rmax = rmaxUpperBound;
                handoverSolution = handover;
                prioritySolution = currentPrioritySolution;
            }
        }
        rmaxUpperBound = rmax;
        handover = bestHandover;
    }


    /*
     * Post handover constraint
     */
    bool Algorithm::postHandoverConstraint(int windowIndex) {

        std::random_device rd;  // Seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

        // Select a random buffer among the overflowing buffers
        std::uniform_int_distribution<> distrib;
        int randomIndex = -1;
        int bufferIndex;
        double offset;

        // If offset is too small get another buffer
        do {
            // Remove buffer from selection
            if(randomIndex > -1) {
                overflowingBuffersMinViolation.erase(overflowingBuffersMinViolation.begin() + randomIndex);
            }
            // Check if overflowing buffer set is empty
            if(overflowingBuffersMinViolation.size()) {
                distrib = std::uniform_int_distribution<>(0, overflowingBuffersMinViolation.size() - 1);
                randomIndex = distrib(gen);
                bufferIndex = overflowingBuffersMinViolation[randomIndex];
                offset = std::min(singleWindowMinViolation, handover[windowIndex-1][bufferIndex] - handoverLb[windowIndex-1][bufferIndex]);
            } else {
                return false;
            }
        } while(offset < precision);

        handoverUb[windowIndex-1][bufferIndex] = handover[windowIndex-1][bufferIndex] - offset;
        handover[windowIndex-1][bufferIndex] = handoverUb[windowIndex-1][bufferIndex];
        return true;
    }

/*    *//*
     * Relax handover and then try to match the PL handovers as much as possible
     *//*
    bool Algorithm::relaxHandovers(int windowIndex) {

        assert(windowIndex < instance.getNumDownlinks()-1);

        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            handover[windowIndex][b] = std::min(handoverUb[windowIndex][b], target[b] - instance.getInstrument(b).getDataOffset(windowIndex+1));
        }

        return singleWindow(windowIndex, rmax, handover[windowIndex]);
    }*/

    void Algorithm::initHandoverConstraints() {
        for(auto w{0}; w < instance.getNumDownlinks(); ++w) {
            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                handoverUb[w][b] = target[b];
                if(handover[w][b] > handoverUb[w][b]) {
                    handover[w][b] = handoverUb[w][b];
                }
            }
        }
    }

    /*
     * Update handovers to usage of last solution
     */
    void Algorithm::updateHandoverTargets() {
        for(auto w{0}; w < instance.getNumDownlinks() - 1; ++w) {
            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                handover[w][b] = handoverSolution[w][b];
            }
        }
    }

    /*
     * Choose which buffer not to stop (this is the one that is stopped the lastest among the buffers stopped)
     */
    void Algorithm::setNoStop(int windowIndex) {
        double lastStop = 0.0;
        int index = -1;
        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            if(stop[windowIndex][b] < instance.getDownlink(windowIndex).end and stop[windowIndex][b] > lastStop) {
                lastStop = stop[windowIndex][b];
                index = b;
            }
        }
        noStop[index] = true;
    }


/*     *//*
      * Repair Solution to avoid bandwidth loss
      *//*
     void Algorithm::repair() {
         while()
     }*/


    /*
     * Main heuristic: Iteratively solve the PL and try to fit the handovers with SingleWindow
     */
    double Algorithm::solve(double timeLimit) {

        double start = std::clock();

        std::vector<double> initialMemoryState;
        initialMemoryState.resize(instance.getNumInstruments());

        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            initialMemoryState[b] = instance.getInstrument(b).getInitialMemoryState() + instance.getInstrument(b).getDataOffset(0);
        }

        // Init rmaxUB
        rmaxUpperBound = 1.0;
        // Run DownlinkCount heuristic
        runDownlinkCount(initialMemoryState, 100);
        // Run the PL
        // runPL(initialMemoryState, -1);
        double startPl = std::clock();

        std::cout << "time=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
        std::cout << "rmaxUB=" << rmaxUpperBound << std::endl;

        CplexFlowDataTransferModel model(instance, mLogger);
        dataflow::TransferSolution solution = model.Solve();
        if(solution.isSolution()) {
            rmaxLowerBound = solution.getRmax();
            handover = solution.getHandovers();
            handoverSolution = handover;
            std::cout << "timePL=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
            std::cout << "rmaxPL=" << rmaxLowerBound << std::endl;
        } else {
            mLogger.Debug("No solution with these handover constraints");
        }

        mLogger.Debug("rmax lower bound: ", rmaxLowerBound);

        timePL += std::clock() - startPl;

        // Improvement step size
        double epsilon = 0.0001;
        double epsilonScaleFactor;
        rmax = rmaxUpperBound - epsilon;

        if(rmaxUpperBound <= rmaxLowerBound + epsilon) {
            std::cout << "time=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
            std::cout << "rmax=" << rmax << std::endl;
            std::cout << "status=OPTIMAL" << std::endl;
            return rmaxUpperBound;
        }


        while((std::clock() - start) / CLOCKS_PER_SEC < timeLimit) {

            // Reset target
            rmax = rmaxUpperBound - epsilon;
            epsilonScaleFactor = 1;

            int currentDownlinkIndex = 0;\
            setTargetRmax(rmax);
            // Update handover targets to integrate the new handover constraints
            updateHandoverTargets();
            // Reset the handover constraints
            initHandoverConstraints();

            // Main loop, successively test current objective on each window,
            // local search for the unsat windows
            while(true) {

                if((std::clock() - start) / CLOCKS_PER_SEC >= timeLimit) {
                    break;
                }

                // Init memory state for the first window
                if(currentDownlinkIndex == 0) {
                    for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                        usage[0][b] = instance.getInstrument(b).getInitialMemoryState();
                    }
                }

                // Add the non-visibility values
                for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                    usage[currentDownlinkIndex][b] += instance.getInstrument(b).getDataOffset(currentDownlinkIndex);
                }

                // Single window checks if current target is sat according to the memory and handover state
                if(singleWindow(currentDownlinkIndex, rmax, handover[currentDownlinkIndex], handoverUb[currentDownlinkIndex])) {

                    if(currentDownlinkIndex == instance.getNumDownlinks() - 1) {
                        priorityTmp[currentDownlinkIndex] = bestPriority;
                        prioritySolution = priorityTmp;
                        handoverSolution = handover;

                        ++nbLoops;
                        assert(checkSolution(rmax));
                        rmaxUpperBound = rmax;
                        if (rmaxUpperBound < rmaxLowerBound + epsilon) {
                            timeTotal = std::clock() - start;
                            std::cout << "time=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
                            std::cout << "rmax=" << rmax << std::endl;
                            std::cout << "status=OPTIMAL" << std::endl;
                            return rmaxUpperBound;
                        } else {
                            mLogger.Debug("New solution with rmax: ", rmax);
                            std::cout << "time=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
                            std::cout << "rmax=" << rmax << std::endl;

                            // Check usage sum
                            double sum = 0.0;
                            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                                sum += usage[currentDownlinkIndex][b];
                            }
                            std::cout << "expected=" << instance.getMemState(instance.getNumTimepoints() - 1) << std::endl;
                            std::cout << "actual=" << sum << std::endl;

                            // Go back to the first window and reduce rmax target
                            rmax -= epsilon * epsilonScaleFactor;
                            epsilonScaleFactor *= 2;
                            setTargetRmax(rmax);
                            // Update handover targets to integrate the new handover constraints
                            updateHandoverTargets();
                            // Reset the handover constraints
                            initHandoverConstraints();
                            currentDownlinkIndex = 0;
                        }
                    } else {
                        // Save solution for current window
                        priorityTmp[currentDownlinkIndex] = bestPriority;
                        handover[currentDownlinkIndex] = usage[currentDownlinkIndex+1];
                        ++currentDownlinkIndex;
                    }

                } else {

                    if(currentDownlinkIndex >= 1) {
                        if (postHandoverConstraint(currentDownlinkIndex)) {
                            --currentDownlinkIndex;
                            // Remove the non-visibility values of previous window
                            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                                usage[currentDownlinkIndex][b] -= instance.getInstrument(b).getDataOffset(currentDownlinkIndex);
                            }
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        mLogger.Debug("status: FEASIBLE");
        std::cout << "time=" << (std::clock() - start) / CLOCKS_PER_SEC << std::endl;
        std::cout << "rmax=" << rmax << std::endl;
        std::cout << "status=FEASIBLE" << std::endl;
        return rmaxUpperBound;
    }


} // namespace