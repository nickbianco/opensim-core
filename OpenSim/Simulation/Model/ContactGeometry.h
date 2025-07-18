#ifndef OPENSIM_CONTACT_GEOMETRY_H_
#define OPENSIM_CONTACT_GEOMETRY_H_ 
/* -------------------------------------------------------------------------- *
 *                        OpenSim:  ContactGeometry.h                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2025 Stanford University and the Authors                *
 * Author(s): Peter Eastman                                                   *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <OpenSim/Simulation/osimSimulationDLL.h>
#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/PhysicalFrame.h"
#include "Appearance.h"

namespace OpenSim {

class ScaleSet;

/** 
 * \section ContactGeometry
 * A base class that represents the physical shape of an object for use in 
 * contact modeling and path wrapping.
 *
 * Concrete implementations of `ContactGeometry` define particular geometric
 * representations (e.g., spheres, cylinders, etc.). The geometry is attached to
 * a PhysicalFrame, which is specified using a Socket named "frame".
 *
 * @note ContactGeometry is not scaled with the Model.
 *
 * @author Peter Eastman
 */
class OSIMSIMULATION_API ContactGeometry : public ModelComponent {
OpenSim_DECLARE_ABSTRACT_OBJECT(ContactGeometry, ModelComponent);
public:
//=============================================================================
// SOCKETS
//=============================================================================
    OpenSim_DECLARE_SOCKET(frame, PhysicalFrame,
        "The frame to which this contact geometry is attached.");

//=============================================================================
// PROPERTIES
//=============================================================================
    OpenSim_DECLARE_PROPERTY(location, SimTK::Vec3,
        "The location of the contact geometry center in the PhysicalFrame.");
    OpenSim_DECLARE_PROPERTY(orientation, SimTK::Vec3,
        "The orientation of the contact geometry in the PhysicalFrame in "
        "body-fixed XYZ Euler angles.");
    OpenSim_DECLARE_UNNAMED_PROPERTY(Appearance,
        "The default appearance for this Geometry");

//=============================================================================
// METHODS
//=============================================================================

    // CONSTRUCTION
    /** Construct an empty ContactGeometry. */
    ContactGeometry();

    /** This constructor connects this ContactGeometry to the provided `frame`,
     * and uses the default location and orientation (both `SimTK::Vec3(0)`).
     *
     * @param frame        the PhysicalFrame this geometry is attached to.
     */
    explicit ContactGeometry(const PhysicalFrame& frame);

    /** This constructor connects this ContactGeometry to the provided `frame`,
     * and sets the location and orientation properties. 
     *
     * @param location     the location of the geometry expressed in `frame`.
     * @param orientation  the orientation of the geometry expressed in `frame`
     *                     as XYZ body-fixed Euler angles.
     * @param frame        the PhysicalFrame this geometry is attached to;
     *                     this constructor connects this ContactGeometry to
     *                     the provided `frame`.
     */
    ContactGeometry(const SimTK::Vec3& location,
                    const SimTK::Vec3& orientation,
                    const PhysicalFrame& frame);

    //** @name Accessors */
    // @{
    /** Get the PhysicalFrame this geometry is attached to. */
    const PhysicalFrame& getFrame() const;
    /** %Set the PhysicalFrame this geometry is attached to. */
    void setFrame(const PhysicalFrame& frame);

    /** Get a Transform representing the position and orientation of the
     * geometry relative to the PhysicalFrame `F` to which this geometry is
     * connected.
     *
     * If you want the transform of this geometry relative to the Frame (or
     * Ground) `B` in which this geometry is fixed, you can use the following
     * code:
     * @code{.cpp}
     * const auto& X_BF = geom.getFrame().findTransformInBaseFrame();
     * const auto X_FP = geom.getTransform();
     * const auto X_BP = X_BF * X_FP;
     * @endcode
     *
     * Prior to OpenSim 4.0, there wwas no intermediate PhysicalFrame `F`, so
     * this method essentially returned `X_BP`. */
    SimTK::Transform getTransform() const;

    /** Create a new SimTK::ContactGeometry based on this object. */
    SimTK::ContactGeometry createSimTKContactGeometry() const;

    /** Get a shared pointer to a SimTK::ContactGeometry based on this object. */
    std::shared_ptr<const SimTK::ContactGeometry> 
    getSimTKContactGeometryPtr() const;
    // @}

    /** @name Visualization */
    // @{
    void generateDecorations(bool fixed, const ModelDisplayHints& hints,
        const SimTK::State& s,
        SimTK::Array_<SimTK::DecorativeGeometry>& geometry) const override;
    // @}

    /** @name Deprecated */
    // @{
    /** <b>(Deprecated)</b> Use get_location() instead. */
    DEPRECATED_14("use get_location() instead")
    const SimTK::Vec3& getLocation() const;

    /** <b>(Deprecated)</b> Use set_location() instead. */
    DEPRECATED_14("use set_location() instead")
    void setLocation(const SimTK::Vec3& location);

    /** <b>(Deprecated)</b> Use get_orientation() instead. */
    DEPRECATED_14("use get_orientation() instead")
    const SimTK::Vec3& getOrientation() const;

    /** <b>(Deprecated)</b> Use set_orientation() instead. */
    DEPRECATED_14("use set_orientation() instead")
    void setOrientation(const SimTK::Vec3& orientation);

    /** <b>(Deprecated)</b> Use getFrame() instead.
     * Get the Body this geometry is attached to. */
    DEPRECATED_14("use getFrame() instead")
    const PhysicalFrame& getBody() const;

    /** <b>(Deprecated)</b> Use setFrame() instead.
     * %Set the Body this geometry is attached to. */
    DEPRECATED_14("use setFrame() instead")
    void setBody(const PhysicalFrame& body);
    // @}

protected:
    // Concrete implementations of ContactGeometry must implement this method.
    virtual SimTK::ContactGeometry createSimTKContactGeometryImpl() const = 0;

    // OBJECT INTERFACE
    void updateFromXMLNode(SimTK::Xml::Element& node, int versionNumber)
        override; 

private:
    mutable std::shared_ptr<const SimTK::ContactGeometry> _simTKContactGeometry;

    // INITIALIZATION
    void setNull();
    void constructProperties();

};

} // namespace OpenSim

#endif // OPENSIM_CONTACT_GEOMETRY_H_ 
