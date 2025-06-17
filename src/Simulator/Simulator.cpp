//
// Created by Julien Rouzot on 20/05/24.
//

#include <algorithm>
#include <iostream>
#include <vector>

#include "Simulator.h"

using namespace dataflow;

namespace simulation {

    bool Event::operator<(const Event &e) const {
        return this->time < e.time;
    }

    bool Event::operator==(const Event &e) const {
        return this->time == e.time;
    }

    bool Event::operator<=(const Event &e) const {
        return this->time <= e.time;
    }

    bool Event::operator>(const Event &e) const {
        return this->time > e.time;
    }

    bool Event::operator>=(const Event &e) const {
        return this->time >= e.time;
    }


    Simulator::Simulator(const Instance &i, logger::Logger &logger) : instance(i), mLogger(logger), now(0), end(-1) {

        timeOverflow = -1;
        bufferOverflow = -1;
        bufferIndex = -1;
        dumpRate = -1;
        overflowTimes.resize(instance.getNumInstruments(), never);
        maxUsage.resize(instance.getNumInstruments(), 0);
        current.resize(instance.getNumInstruments(), {0, 0, 0});
        noStop.resize(instance.getNumInstruments(), false);
        deactivatedBuffers.reserve(instance.getNumInstruments());
        // Init capacities
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            capacity.push_back(instance.getInstrument(b).getCapacity());
        }
    }

    // Initialise all variables
    void Simulator::initialise(std::vector<int> initialPointer, std::vector<double> initialUsage) {

        bufferIndex = -1;
        overflowTimes.clear();
        timeOverflow = never;
        bufferOverflow = -1;
        std::fill(noStop.begin(), noStop.end(), false);
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            current[b].pointer = initialPointer[b];
            current[b].usage = initialUsage[b];
            maxUsage[b] = initialUsage[b];
            current[b].transferRate = 0;
        }
        for(int i{0}; i < instance.getNumInstruments(); ++i) {
            deactivatedBuffers[i] = false;
        }
    }

    const double Simulator::never = std::numeric_limits<double>::max();

    double Simulator::getMaxUsage(const int b) const {
        return maxUsage[b];
    }

    double Simulator::getFillRate(const int b) {
        return instance.getInstrument(b).getEvents()[current[b].pointer].value;
    }

    double Simulator::getTransferRate(const int b) {
        return current[b].transferRate;
    }

    double Simulator::getUsage(const int b) const {
        return current[b].usage;
    }

    int Simulator::getPointer(const int b) const {
        return current[b].pointer;
    }

    void Simulator::setNoStop(std::vector<bool> &pNoStop) {
        noStop = pNoStop;
    }


    // Return the time when buffer b will be empty at the current rate and usage
    double Simulator::getEmptyTime(const int b) {

        // If buffer is filled faster than its transfer, or if it is deactivated, it will never get empty
        if (deactivatedBuffers[b] or getFillRate(b) >= getTransferRate(b)) {
            return Simulator::never;
        }
        if (getUsage(b) < precision) {
            mLogger.Debug("Buffer: ", b);
            mLogger.Debug("Usage: ", getUsage(b));
            mLogger.Debug("Usage < precision");
            exit(1);
        }
        // Otherwise, the buffer with be empty at usage / actual transfer rate (dump-fill)
        double t = now + getUsage(b) / (getTransferRate(b) - getFillRate(b));
        if (t <= now) {
            return Simulator::never;
        }
        return t;
    }

    // Return the earliest stop time for buffer @b to respect handover constraint
    double Simulator::getEarliestStopTime(int currentBufferIndex, double next) {

        assert(next <= end);
        assert(next >= 1);

        // Get the amount of data we have to dump at next
        double toDump(instance.getInstrument(currentBufferIndex).getDataToDumpAt(next));

/*
        assert(handover[currentBufferIndex] <= capacity[currentBufferIndex] + precision);
*/

        if(toDump >= handover[currentBufferIndex] + precision) {
            return never;
        }

        if(getTransferRate(currentBufferIndex) > 0) {
            double stopTime = now + (getUsage(currentBufferIndex) - (handover[currentBufferIndex] - toDump) + getFillRate(currentBufferIndex)*(next - now)) / getTransferRate(currentBufferIndex);
            return stopTime;
        // Edge case when buffer can be stopped right away
        } else if(getUsage(currentBufferIndex) + toDump + getFillRate(currentBufferIndex)*(next - now) <= handover[currentBufferIndex]) {
            return now;
        } else {
            return never;
        }
    }


    // Compute the rates at the current timepoint according to the priority assignment
    void Simulator::computeRates(std::vector<std::vector<int>> &priority) {

        double bandwidth = dumpRate;
        double bandwidthAlloc;
        auto batch{priority.rbegin()};

        while (batch != priority.rend()) {

            if (bandwidth >= precision) {
                // Sort buffers within the same priority pool by usage (or by fill rate if they are empty)
                std::sort(batch->begin(), batch->end(), [&](const int a, const int b) {

                    auto ua{getUsage(a)};
                    auto ub{getUsage(b)};
                    auto fa(getFillRate(a));
                    auto fb(getFillRate(b));
                    // Deactivated buffers must go first
                    if(deactivatedBuffers[a] && !deactivatedBuffers[b]) {
                        return true;
                    }
                    // Deactivated buffers must go first
                    if(!deactivatedBuffers[a] && deactivatedBuffers[b]) {
                        return false;
                    }
                    if (ua < precision and ub < precision)
                        return fa < fb;
                    return ua < ub;
                });
                // The bandwidth allocated to the buffers is the remaining bandwidth divided among the buffers with same priority
                bandwidthAlloc = (bandwidth / batch->size());

            } else {
                bandwidth = 0;
                bandwidthAlloc = 0;
            }

            for (auto b{batch->begin()}; b != batch->end(); ++b) {

                // Bandwidth available is none
                if (bandwidth < precision) {
                    current[*b].transferRate = 0;
                    // Buffer is deactivated -> transferRate is null
                } else if(deactivatedBuffers[*b]) {
                    current[*b].transferRate = 0;
                    auto d{batch->end() - b - 1};
                    if (d > 0)
                        bandwidthAlloc = bandwidth / d;
                    // If buffer b is non-empty or its bandwidth share is less than its fill rate
                } else if (getUsage(*b) > precision or bandwidthAlloc < getFillRate(*b)) {
                    current[*b].transferRate = bandwidthAlloc;         // buffer b takes all its bandwidth
                    bandwidth -= getTransferRate(*b);                  // we remove it from the overall bandwidth
                    // b is empty and bandwidth is share is greater than the fill rate
                } else {
                    current[*b].transferRate = getFillRate(*b);        // buffer b takes only the bandwidth necessary to compensate its fill rate
                    bandwidth -= getTransferRate(*b);                  // we remove it from the overall bandwidth
                    auto d{batch->end() - b - 1};
                    if (d > 0)
                        bandwidthAlloc = bandwidth / d;
                }
            }
            ++batch;
        }
/*
        // DEBUG
        double sumTransfer = 0.0;
        bool allEmpty = true;
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            sumTransfer += current[b].transferRate;
            if(current[b].usage > precision) {
                allEmpty = false;
            }
        }

        // RELAX
        if(not(allEmpty) and sumTransfer + precision < dumpRate) {
            mLogger.Debug("[SIMULATION] Bandwidth loss");
            for(auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(current[b].usage > precision) {
                    current[b].transferRate = sumTransfer;
                }
            }
        }*/
    }


    // Go to the next timepoint and update memory usage (and max memory usage).
    // If a buffer overflows, add it to the excessBuffer set.
    void Simulator::advanceToStop(const double t) {

        for (auto b{0}; b < instance.getNumInstruments(); ++b) {

            // Check the overflow time
            if(getFillRate(b) > getTransferRate(b)) {
                overflowTimes[b] = (capacity[b] - current[b].usage) / (getFillRate(b) - getTransferRate(b)) + now;
            } else {
                overflowTimes[b] = never;
            }
            // Calculate the memory state for the next time point
            current[b].usage += (getFillRate(b) - getTransferRate(b)) * (t - now);

            // Usage can't be negative
            if (current[b].usage < precision) {
                current[b].usage = 0;
            }
            // Update max usage
            maxUsage[b] = std::max(maxUsage[b], getUsage(b));
        }

        // Update the first overflow
        timeOverflow = never;
        bufferOverflow = -1;
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            if(overflowTimes[b] < timeOverflow) {
                timeOverflow = overflowTimes[b];
                bufferOverflow = b;
            }
        }
        now = t;
    }


    // Run without stop
    // Rmax and handovers are not constrained
    void Simulator::run(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd,
             double currentDumpRate, double startTime, double endTime, std::vector<std::vector<int>> &priority) {

        // Update dumpRate
        dumpRate = currentDumpRate;
        // Now variable keeps track of the current timepoint
        now = startTime;
        // End time for the window
        end = endTime;
        // Compute the initial rates
        computeRates(priority);

        // This is the POINTER to the current event in the event list
        auto cur{eventBegin};

        while(cur->time == now and cur < eventEnd) {
            ++cur;
        }

        // DEBUG
/*        if (std::binary_search(instance.getTimepoints().begin(), instance.getTimepoints().end(), now)) {
            double curMemState = 0.0;
            for(int b{0}; b < instance.getNumInstruments(); ++b) {
                curMemState += current[b].usage;
            }
            auto index = std::distance(instance.getTimepoints().begin(), std::lower_bound(instance.getTimepoints().begin(), instance.getTimepoints().end(), now));
            assert(curMemState - precision <= instance.getMemState(index) and instance.getMemState(index) <= curMemState + precision);
        }*/

        // Time for next step
        double next;
        next = cur->time;
        if(cur < eventEnd) {
            next = cur->time;
        } else {
            next = end;
        }

        // Add the "empty memory" events
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            if(getUsage(b) >= precision) {
                auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                if (t < next) {
                    next = t;
                }
            }
        }

        // Main loop
        while (cur < eventEnd or next < endTime) {

            // Jump to next (apply all now events)
            advanceToStop(next);

/*            // DEBUG
            if (std::binary_search(instance.getTimepoints().begin(), instance.getTimepoints().end(), next)) {
                double curMemState = 0.0;
                for(int b{0}; b < instance.getNumInstruments(); ++b) {
                    curMemState += current[b].usage;
                }
                auto index = std::distance(instance.getTimepoints().begin(), std::lower_bound(instance.getTimepoints().begin(), instance.getTimepoints().end(), next));
                assert(curMemState - precision <= instance.getMemState(index) and instance.getMemState(index) <= curMemState + precision);
            }*/

            // Apply all events up to now
            while (cur->time == now) {
                /*mLogger.Debug("Current event ", cur->time, ",", cur->buffer);*/
                // Check event
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;
                /*mLogger.Debug("Updating pointer for buffer ", e.buffer, " to ", current[e.buffer].pointer);*/
            }
            /*mLogger.Debug("Updating next to ", cur->time);*/
            next = cur->time;
            if(cur < eventEnd) {
                next = cur->time;
            } else {
                next = end;
            }

            computeRates(priority);

            // Add the "empty memory" events
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(getUsage(b) >= precision) {
                    auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                    if (t < next) {
                        next = t;
                    }
                }
            }
        }

        // Advance to the last timepoint
        advanceToStop(endTime);
    }


    /*
     * Run the simulation but do not exit when overflow is reached.
     * This allows to compute the violation for each buffer.
     */
    void Simulator::runWithEarliestStopViolation(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd,
                                                 double currentDumpRate, double startTime, double endTime, std::vector<std::vector<int>> &priority,
                                                 std::vector<double> &stop, std::vector<double> &target, std::vector<double> &htarget) {
        // Reset the stop times
        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            stop[b] = never;
        }
        // Update dumpRate
        dumpRate = currentDumpRate;
        // Target is the maximum peak we allow, for each buffer
        swap(capacity, target);
        // Handover constraints are the memory limit at the end of the instance
        swap(handover, htarget);
        // Now variable keeps track of the current timepoint
        now = startTime;
        // End time for the window
        end = endTime;
        // This is the POINTER to the current event in the event list
        auto cur{eventBegin};

        // Compute the initial rates
        computeRates(priority);

        while(cur->time == now and cur < eventEnd) {
            ++cur;
        }

        // Time for next step
        double next;
        next = cur->time;
        if(cur < eventEnd) {
            next = cur->time;
        } else {
            next = end;
        }

        // Stop all buffers that can be stopped from the beginning
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {

            auto t{getEarliestStopTime(b, next)};    // Check when we can stop the buffer
            if (t <= now) {
                stop[b] = now;
                deactivatedBuffers[b] = true;
                current[b].transferRate = 0;
            }
        }

        // Add the "empty memory" events
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            if(getUsage(b) >= precision) {
                auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                if (t < next) {
                    next = t;
                }
            }
        }
        // Add the early stop event
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {

            if (!deactivatedBuffers[b]) {
                auto t{getEarliestStopTime(b, next)};    // Check if when we can stop the buffer
                if (t < now) {
                    t = now;
                }
                if (t + precision < next) {
                    bufferIndex = b;
                    next = t;
                }
            }
        }

        while (cur < eventEnd or next < endTime) {

            // Jump to next (apply all now events)
            advanceToStop(next);

            // Stop the buffer that can be stopped the earliest
            if (bufferIndex >= 0) {
                stop[bufferIndex] = now;
                deactivatedBuffers[bufferIndex] = true;
                current[bufferIndex].transferRate = 0;
                bufferIndex = -1;
            }

            // Apply all events up to now
            while (cur->time == now) {
                // Check event
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;
            }
            next = cur->time;                   // Next timepoint
            if(cur < eventEnd) {
                next = cur->time;
            } else {
                next = end;
            }
            computeRates(priority);

            // Add the "empty memory" events
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(getUsage(b) >= precision) {
                    auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                    if (t < next) {
                        next = t;
                    }
                }
            }
            // Add the early stop event
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {

                if (!deactivatedBuffers[b]) {
                    auto t{getEarliestStopTime(b, next)};    // Check if when we can stop the buffer
                    if (t < now) {
                        t = now;
                    }
                    if (t + precision < next) {
                        bufferIndex = b;
                        next = t;
                    }
                }
            }
        }

        advanceToStop(endTime);

        swap(capacity, target);
        swap(handover, htarget);
    }


    // Main loop for the Simulator
    // Here, we interrupt the buffers as soon as possible
    int Simulator::runWithEarliestStop(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd,
                                       double currentDumpRate, double startTime, double endTime, std::vector<std::vector<int>> &priority,
                                       std::vector<double> &stop, std::vector<double> &target, std::vector<double> &hTarget, std::vector<double> &hConstraints) {

        // Reset the stop times
        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            stop[b] = never;
        }
        // Check if the initial memory state respect the rmax
        for(int b{0}; b < instance.getNumInstruments(); ++b) {
            if(getUsage(b) >= target[b] + precision) {
                return b;
            }
        }
        // Update dumpRate
        dumpRate = currentDumpRate;
        // Target is the maximum peak we allow, for each buffer
        swap(capacity, target);
        // Handover target are the memory target at the window's end
        swap(handover, hTarget);
        // Handover constraints are the memory state limit at the window's end
        swap(handoverCt, hConstraints);
        // Now variable keeps track of the current timepoint
        now = startTime;
        // End time for the window
        end = endTime;
        // This is the POINTER to the current event in the event list
        auto cur{eventBegin};

        // Compute the initial rates
        computeRates(priority);

        // Check all events until next time step
        while(cur->time == now and cur < eventEnd) {
            ++cur;
        }

        // Current time
        double next;
        if(cur < eventEnd) {
            next = cur->time;
        } else {
            next = end;
        }

        // Stop all buffers that can be stopped from the beginning
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            auto t{getEarliestStopTime(b, next)};    // Check when we can stop the buffer
            if (t <= now) {
                stop[b] = now;
                deactivatedBuffers[b] = true;
                current[b].transferRate = 0;
            }
        }

        // Add the "empty memory" events
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            if(getUsage(b) >= precision) {
                auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                if (t < next) {
                    next = t;
                }
            }
        }
        // Add the early stop event
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {

            if(not(deactivatedBuffers[b])) {
                auto t{getEarliestStopTime(b, next)};    // Check if when we can stop the buffer
                if(t < now) {
                    t = now;
                }
                if (t + precision < next) {
                    bufferIndex = b;
                    next = t;
                }
            }
        }

        while (cur < eventEnd or next < endTime) {

            // Jump to next (apply all now events)
            advanceToStop(next);

            // Return if we reach the first overflow
            if(timeOverflow < now) {
                swap(capacity, target);
                swap(handover, hTarget);
                swap(handoverCt, hConstraints);
                return bufferOverflow;
            }

            // Stop the buffer that can be stopped the earliest
            if (bufferIndex >= 0) {
                stop[bufferIndex] = now;
                deactivatedBuffers[bufferIndex] = true;
                current[bufferIndex].transferRate = 0;
                bufferIndex = -1;
            }

            // Apply all events up to now
            while (cur->time == now) {
                // Check event
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;
            }
            // Next timepoint
            if(cur < eventEnd) {
                next = cur->time;
            } else {
                next = end;
            }

            computeRates(priority);

            // Add the "empty memory" events
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(getUsage(b) >= precision) {
                    auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                    if (t < next) {
                        next = t;
                    }
                }
            }
            // Add the early stop event
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {

                if(not(deactivatedBuffers[b])) {
                    auto t{getEarliestStopTime(b, next)};    // Check if when we can stop the buffer
                    if(t < now) {
                        t = now;
                    }
                    if (t + precision < next) {
                        bufferIndex = b;
                        next = t;
                    }
                }
            }
        }

        advanceToStop(endTime);

        // Return if we reach the first overflow
        if(timeOverflow < now) {
            swap(capacity, target);
            swap(handover, hTarget);
            swap(handoverCt, hConstraints);
            return bufferOverflow;
        }

        // Return if we reach overflow
        if(timeOverflow + precision <= endTime) {
            swap(capacity, target);
            swap(handoverCt, hConstraints);
            swap(handover, hTarget);
            return bufferOverflow;
        }

        // Check handovers overflows
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {

            if (getUsage(b) >= handoverCt[b] + precision) {
                swap(capacity, target);
                swap(handover, hTarget);
                swap(handoverCt, hConstraints);
                return b;
            }
        }

        swap(capacity, target);
        swap(handover, hTarget);
        swap(handoverCt, hConstraints);

        return -1;
    }

    /*
     * Run with full bandwidth allocated to each buffer
     */
    void Simulator::runFullBandwidth(std::vector<Event>::const_iterator eventBegin, std::vector<Event>::const_iterator eventEnd,
                                     double currentDumpRate, double startTime, double endTime) {

        // Update dumpRate
        dumpRate = currentDumpRate;
        // Now variable keeps track of the current timepoint
        now = startTime;
        // End time for the window
        end = endTime;
        // This is the POINTER to the current event in the event list
        auto cur{eventBegin};

        while(cur->time == now and cur < eventEnd) {
            ++cur;
        }

        // Time for next step
        double next;
        next = cur->time;
        if(cur < eventEnd) {
            next = cur->time;
        } else {
            next = end;
        }

        // Compute the initial rates
        for(auto b{0}; b < instance.getNumInstruments(); ++b) {
            current[b].transferRate = dumpRate;
        }

        // Add the "empty memory" events
        for (auto b{0}; b < instance.getNumInstruments(); ++b) {
            if(getUsage(b) >= precision) {
                auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                if (t < next) {
                    next = t;
                }
            }
        }

        while (cur < eventEnd or next < endTime) {

            // Jump to next (apply all now events)
            advanceToStop(next);

            while (cur->time == now) {
                // Check event
                auto e{*cur};
                ++cur;
                ++current[e.buffer].pointer;
            }
            next = cur->time;                   // Next timepoint
            if(cur < eventEnd) {
                next = cur->time;
            } else {
                next = end;
            }

            // Add the "empty memory" events
            for (auto b{0}; b < instance.getNumInstruments(); ++b) {
                if(getUsage(b) >= precision) {
                    auto t{getEmptyTime(b)};        // Compute the timepoint for the buffer to empty
                    if (t < next) {
                        next = t;
                    }
                }
            }
        }
        advanceToStop(endTime);
    }
}

