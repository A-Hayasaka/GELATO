//
//  gravity.hpp
//  OpenTsiolkovsky
//
//  Created by Takahiro Inagawa on 2015/11/23.
//  Copyright © 2015 Takahiro Inagawa. All rights reserved.
//

/**
 * @file gravity.hpp
 * @brief Gravity model calculations for trajectory simulation
 * 
 * This header provides gravity acceleration calculations in the ECI frame.
 * Two models are available:
 * - Full model with J2 oblateness perturbation (gravityECI)
 * - Simplified two-body model (gravityECI_simple)
 * 
 * Based on WGS84 EGM96 gravitational parameters.
 */

#ifndef SRC_GRAVITY_HPP_
#define SRC_GRAVITY_HPP_

#include <Eigen/Core>
#include <cmath>

/**
 * @brief Calculate gravity acceleration in ECI frame with J2 perturbation
 * 
 * Computes the gravitational acceleration vector including the J2 oblateness
 * effect of the Earth. This provides higher accuracy than the simple two-body
 * model for trajectory calculations.
 * 
 * @param posECI_ Position vector in ECI frame [m]
 * @return Gravity acceleration vector in ECI frame [m/s²]
 */
Eigen::Vector3d gravityECI(Eigen::Vector3d posECI_);

/**
 * @brief Calculate gravity acceleration in ECI frame (simple two-body model)
 * 
 * Computes the gravitational acceleration vector using a simple inverse-square
 * law (two-body problem). Does not include Earth oblateness effects.
 * 
 * @param posECI_ Position vector in ECI frame [m]
 * @return Gravity acceleration vector in ECI frame [m/s²]
 */
Eigen::Vector3d gravityECI_simple(Eigen::Vector3d posECI_);

#endif  // SRC_GRAVITY_HPP_
