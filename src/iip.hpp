//
// The MIT License
//
// Copyright (c) 2024 Interstellar Technologies Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

/**
 * @file iip.hpp
 * @brief Instantaneous Impact Point (IIP) calculation
 * 
 * This header provides the IIP calculation function using the FAA
 * (Federal Aviation Administration) simplified method. The IIP is the
 * point on Earth's surface where a vehicle would impact if thrust were
 * terminated at the current instant, considering only ballistic motion.
 * reference: 14 CFR Appendix-B-to-Part-420(d)(3)(v) 
 * https://www.ecfr.gov/current/title-14/appendix-Appendix%20B%20to%20Part%20420#p-Appendix-B-to-Part-420(d)(3)(v)
 */

#ifndef SRC_IIP_HPP_
#define SRC_IIP_HPP_

#include <Eigen/Core>

/**
 * @brief Calculate Instantaneous Impact Point (IIP) using FAA method
 * 
 * Computes the instantaneous impact point, which is the point on the Earth's
 * surface where the vehicle would impact if thrust were terminated immediately.
 * Uses the FAA (Federal Aviation Administration) simplified method.
 * reference: 14 CFR Appendix-B-to-Part-420(d)(3)(v) 
 * https://www.ecfr.gov/current/title-14/appendix-Appendix%20B%20to%20Part%20420#p-Appendix-B-to-Part-420(d)(3)(v)
 * 
 * @param posECEF_ Current position vector in ECEF frame [m]
 * @param velECEF_ Current velocity vector in ECEF frame [m/s]
 * @return IIP position in geodetic coordinates (latitude [rad], longitude [rad], height [m])
 */
Eigen::Vector3d posLLH_IIP_FAA(Eigen::Vector3d posECEF_,
                               Eigen::Vector3d velECEF_);

#endif  // SRC_IIP_HPP_
