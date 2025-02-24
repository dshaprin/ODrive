#include <cmath>
#include "odrive_main.h"
#include "utils.hpp"

// Symbol                     Description
// Ta, Tv and Td              Duration of the stages of the AL profile
// Xi and Vi                  Adapted initial conditions for the AL profile
// Xf                         Position set-point
// s                          Direction (sign) of the trajectory
// Vmax, Amax, Dmax and jmax  Kinematic bounds
// Ar, Dr and Vr              Reached values of acceleration and velocity

bool SineAccelerationTrajectory::plan(float Xf, float Xi, float Vi,
                                            float Vmax, float Amax, float Dmax) {
                                                
    float dX = Xf - Xi;  // Distance to travel
    float stop_dist = (Vi * Vi) / (2.0f * Dmax); // Minimum stopping distance
    float dXstop = std::copysign(stop_dist, Vi); // Minimum stopping displacement
    float s = sign_hard(dX - dXstop); // Sign of coast velocity (if any)
    Ar_ = s * Amax;  // Maximum Acceleration (signed)
    Dr_ = -s * Dmax; // Maximum Deceleration (signed)
    Vr_ = s * Vmax;  // Maximum Velocity (signed)

    // If we start with a speed faster than cruising, then we need to decel instead of accel
    // aka "double deceleration move" in the paper
    if ((s * Vi) > (s * Vr_)) {
        Ar_ = -s * Amax;
    }

    // Time to accel/decel to/from Vr (cruise speed)
    Ta_ = (Vr_ - Vi) / Ar_;
    Td_ = -Vr_ / Dr_;

    // Integral of velocity ramps over the full accel and decel times to get
    // minimum displacement required to reach cuising speed
    float dXmin = 0.5f*Ta_*(Vr_ + Vi) + 0.5f*Td_*Vr_;

    // Are we displacing enough to reach cruising speed?
    if (s*dX < s*dXmin) {
        // Short move (triangle profile)
        Vr_ = s * std::sqrt(std::max((Dr_*SQ(Vi) + 2*Ar_*Dr_*dX) / (Dr_ - Ar_), 0.0f));
        Ta_ = std::max(0.0f, (Vr_ - Vi) / Ar_);
        Td_ = std::max(0.0f, -Vr_ / Dr_);
        Tv_ = 0.0f;
    } else {
        // Long move (trapezoidal profile)
        Tv_ = (dX - dXmin) / Vr_;
    }

    // Fill in the rest of the values used at evaluation-time
    Tf_ = Ta_ + Tv_ + Td_;
    Xi_ = Xi;
    Xf_ = Xf;
    Vi_ = Vi;
    yAccel_ = Xi + Vi*Ta_ + 0.5f*Ar_*SQ(Ta_); // pos at end of accel phase

    return true;
}


float smoothing_vel(float x) {
    return 0.5f * (1 - std::cos(x * PI));
}

float smoothing_s(float x) {
    return 0.5f * (x - std::sin(PI * x) / PI);
}

float smoothing_a(float x) {
    return std::sin(PI * x);
}

SineAccelerationTrajectory::Step_t SineAccelerationTrajectory::eval(float t) {
    Step_t trajStep;
    if (t < 0.0f) {  // Initial Condition
        trajStep.Y = Xi_;
        trajStep.Yd = Vi_;
        trajStep.Ydd = 0.0f;
    } else if (t < Ta_) {  // Accelerating
        float prop_t = t / Ta_;
        trajStep.Y = Xi_ + Vr_ * Ta_ * smoothing_s(prop_t) + Ta_ * Vi_ * 0.5f;
        trajStep.Yd = Vi_ + Ar_ * Ta_ * smoothing_vel(prop_t);
        trajStep.Ydd = 0.5f * PI * Ar_ * smoothing_a(prop_t);
    } else if (t < Ta_ + Tv_) {  // Coasting
        trajStep.Y = yAccel_ + Vr_ * (t - Ta_);
        trajStep.Yd = Vr_;
        trajStep.Ydd = 0.0f;
    } else if (t < Tf_) {  // Deceleration
        float td = t - Tf_;
        float delta_t = (Ta_ + Tv_ - Tf_);  // Duration of deceleration step
        float prop_t = td / delta_t;
        trajStep.Y = Xf_ + Dr_ * SQ(delta_t) * smoothing_s(prop_t);
        trajStep.Yd = Dr_ * delta_t * smoothing_vel(prop_t);
        trajStep.Ydd = 0.5f * PI * Dr_ * smoothing_a(prop_t);
    } else if (t >= Tf_) {  // Final Condition
        trajStep.Y = Xf_;
        trajStep.Yd = 0.0f;
        trajStep.Ydd = 0.0f;
    } else {
        // TODO: report error here
    }

    return trajStep;
}