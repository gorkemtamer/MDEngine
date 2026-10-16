// thermostat.cpp
//
// The part common to all thermostats.

#include <format>
#include <stdexcept>
#include "thermostat.h"

Thermostat::Thermostat(double T) : T_(T) {
    if (T <= 0) throw std::invalid_argument(std::format("Thermostat: temperature must be positive (got {:g})", T));
}