/* -------------------------------------------------------------------------- *
 *                   OpenSim:  Scholz2015GeometryPath.cpp                     *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2026 Stanford University and the Authors                *
 * Author(s): Nicholas Bianco                                                 *
 * Contributor(s): Pepijn van den Bos, Andreas Scholz                         *
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

#include "Scholz2015GeometryPath.h"

#include <OpenSim/Common/Constant.h>
#include <OpenSim/Simulation/Model/ForceConsumer.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimulationUtilities.h>

#include <algorithm>

using namespace OpenSim;

namespace {

// The location of a moving path point in the base frame B of its parent frame
// F, X_BF * (x(q), y(q), z(q)), from the point's scalar location functions.
class MovingPointFunction : public SimTK::Function_<SimTK::Vec3> {
public:
    MovingPointFunction(const SimTK::Transform& X_BF,
            const Function& x, const Function& y, const Function& z) :
            m_X_BF(X_BF),
            m_x(x.createSimTKFunction()),
            m_y(y.createSimTKFunction()),
            m_z(z.createSimTKFunction()) {}

    SimTK::Vec3 calcValue(const SimTK::Vector& q) const override {
        return m_X_BF * SimTK::Vec3(m_x->calcValue(q), m_y->calcValue(q),
                m_z->calcValue(q));
    }

    SimTK::Vec3 calcDerivative(const SimTK::Array_<int>& derivComponents,
            const SimTK::Vector& q) const override {
        return m_X_BF.R() * SimTK::Vec3(
                m_x->calcDerivative(derivComponents, q),
                m_y->calcDerivative(derivComponents, q),
                m_z->calcDerivative(derivComponents, q));
    }

    int getArgumentSize() const override { return m_x->getArgumentSize(); }

    int getMaxDerivativeOrder() const override {
        return std::min({m_x->getMaxDerivativeOrder(),
                m_y->getMaxDerivativeOrder(), m_z->getMaxDerivativeOrder()});
    }

private:
    SimTK::Transform m_X_BF;
    std::unique_ptr<SimTK::Function> m_x;
    std::unique_ptr<SimTK::Function> m_y;
    std::unique_ptr<SimTK::Function> m_z;
};

} // namespace

//=============================================================================
// SCHOLZ 2015 GEOMETRY PATH POINT
//=============================================================================
Scholz2015GeometryPathPoint::Scholz2015GeometryPathPoint() {
    constructProperty_PathPoint(PathPoint());
}

const PathPoint& Scholz2015GeometryPathPoint::getPathPoint() const {
    return get_PathPoint();
}

const PhysicalFrame& Scholz2015GeometryPathPoint::getParentFrame() const {
    return getPathPoint().getParentFrame();
}

SimTK::Vec3 Scholz2015GeometryPathPoint::calcLocation(
        const SimTK::State& s) const {
    return getPathPoint().getLocation(s);
}

//=============================================================================
// SCHOLZ 2015 GEOMETRY PATH MOVING POINT
//=============================================================================
Scholz2015GeometryPathMovingPoint::Scholz2015GeometryPathMovingPoint() {
    constructProperty_x_location(Constant(0.));
    constructProperty_y_location(Constant(0.));
    constructProperty_z_location(Constant(0.));
}

const PhysicalFrame& Scholz2015GeometryPathMovingPoint::getParentFrame() const {
    return getConnectee<PhysicalFrame>("parent_frame");
}

std::vector<const Coordinate*>
Scholz2015GeometryPathMovingPoint::getCoordinates() const {
    std::vector<const Coordinate*> coordinates;
    const auto& socket = getSocket<Coordinate>("coordinates");
    for (unsigned i = 0; i < socket.getNumConnectees(); ++i) {
        coordinates.push_back(&socket.getConnectee(i));
    }
    return coordinates;
}

SimTK::Vec3 Scholz2015GeometryPathMovingPoint::calcLocation(
        const SimTK::State& s) const {
    const std::vector<const Coordinate*> coordinates = getCoordinates();
    SimTK::Vector q((int)coordinates.size());
    for (int i = 0; i < (int)coordinates.size(); ++i) {
        q[i] = coordinates[i]->getValue(s);
    }
    return SimTK::Vec3(get_x_location().calcValue(q),
            get_y_location().calcValue(q), get_z_location().calcValue(q));
}

void Scholz2015GeometryPathMovingPoint::extendFinalizeFromProperties() {
    Super::extendFinalizeFromProperties();
    const int numCoordinates = (int)getSocket("coordinates").getNumConnectees();
    for (const Function* function : {&get_x_location(), &get_y_location(),
                 &get_z_location()}) {
        OPENSIM_THROW_IF_FRMOBJ(
                function->getArgumentSize() != numCoordinates, Exception,
                "Expected the location functions to take one argument per "
                "coordinate ({}), but function '{}' takes {}.",
                numCoordinates, function->getName(),
                function->getArgumentSize());
    }
}

//=============================================================================
// SCHOLZ 2015 GEOMETRY PATH OBSTACLE
//=============================================================================
Scholz2015GeometryPathObstacle::Scholz2015GeometryPathObstacle() {
    constructProperty_contact_hint(SimTK::Vec3(SimTK::NaN));
}

const ContactGeometry&
Scholz2015GeometryPathObstacle::getContactGeometry() const {
    return getSocket<ContactGeometry>("contact_geometry").getConnectee();
}

const SimTK::Vec3& Scholz2015GeometryPathObstacle::getContactHint() const {
    return get_contact_hint();
}

//=============================================================================
// CONSTRUCTION
//=============================================================================
Scholz2015GeometryPath::Scholz2015GeometryPath() : AbstractGeometryPath() {
    constructProperties();
}

//=============================================================================
// PATH CONFIGURATION
//=============================================================================
void Scholz2015GeometryPath::appendPathPoint(const PhysicalFrame& frame,
        const SimTK::Vec3& location) {

    // Construct a new path point.
    append_path_elements(Scholz2015GeometryPathPoint());
    auto& point = updElement<Scholz2015GeometryPathPoint>(
            getProperty_path_elements().size() - 1);
    point.setName("path_point_" + std::to_string(getNumPathPoints() - 1));

    // Set the path point's PathPoint.
    PathPoint& pathPoint = point.upd_PathPoint();
    pathPoint.set_location(location);
    pathPoint.setParentFrame(frame);

    // Call finalizeFromProperties() to ensure that the path point is recognized
    // as a subcomponent.
    finalizeFromProperties();
}

void Scholz2015GeometryPath::appendPathPoint(const PhysicalFrame& frame,
        const Function& xLocation, const Function& yLocation,
        const Function& zLocation,
        const std::vector<std::string>& coordinatePaths) {

    // Construct a new moving path point.
    append_path_elements(Scholz2015GeometryPathMovingPoint());
    auto& point = updElement<Scholz2015GeometryPathMovingPoint>(
            getProperty_path_elements().size() - 1);
    point.setName("path_point_" + std::to_string(getNumPathPoints() - 1));

    // Set the frame, location functions, and coordinates.
    point.connectSocket_parent_frame(frame);
    point.set_x_location(xLocation);
    point.set_y_location(yLocation);
    point.set_z_location(zLocation);
    for (const auto& path : coordinatePaths) {
        point.updSocket("coordinates").appendConnecteePath(path);
    }

    // Call finalizeFromProperties() to ensure that the path point is recognized
    // as a subcomponent.
    finalizeFromProperties();
}

const Scholz2015GeometryPathAbstractPoint&
Scholz2015GeometryPath::getOrigin() const {
    OPENSIM_THROW_IF_FRMOBJ(getNumPathPoints() == 0, Exception,
            "Tried retrieving the origin point, but the path contains no path "
            "points.");

    const auto* point = tryGetElement<Scholz2015GeometryPathAbstractPoint>(0);
    OPENSIM_THROW_IF_FRMOBJ(point == nullptr, Exception,
            "Tried retrieving the origin point, but the first element in the "
            "path is not a path point.");

    return *point;
}

const Scholz2015GeometryPathAbstractPoint&
Scholz2015GeometryPath::getInsertion() const {
    OPENSIM_THROW_IF_FRMOBJ(getNumPathPoints() < 2, Exception,
            "Tried retrieving the insertion point, but the path contains less "
            "than two path points.");

    const auto* point = tryGetElement<Scholz2015GeometryPathAbstractPoint>(
            getProperty_path_elements().size() - 1);
    OPENSIM_THROW_IF_FRMOBJ(point == nullptr, Exception,
            "Tried retrieving the insertion point, but the last element in the "
            "path is not a path point.");

    return *point;
}

const Scholz2015GeometryPathAbstractPoint&
Scholz2015GeometryPath::getPathPoint(int pathPointIndex) const {
    OPENSIM_THROW_IF_FRMOBJ(
            pathPointIndex < 0 || pathPointIndex >= getNumPathPoints(),
            Exception,
            "Index {} is out of range. There are {} path points in the path.",
            pathPointIndex, getNumPathPoints());

    int count = 0;
    for (int i = 0; i < getNumPathElements(); ++i) {
        if (const auto* point =
                    tryGetElement<Scholz2015GeometryPathAbstractPoint>(i)) {
            if (count == pathPointIndex) {
                return *point;
            }
            ++count;
        }
    }

    // We should never reach this point.
    OPENSIM_THROW_FRMOBJ(Exception,
            "Path point index {} could not be found.", pathPointIndex);
}

void Scholz2015GeometryPath::appendObstacle(
        const ContactGeometry& contactGeometry,
        const SimTK::Vec3& contactHint) {

    // Construct a new obstacle.
    append_path_elements(Scholz2015GeometryPathObstacle());
    auto& obstacle = updElement<Scholz2015GeometryPathObstacle>(
            getProperty_path_elements().size() - 1);
    obstacle.setName("obstacle_" + std::to_string(getNumObstacles() - 1));

    // Set the obstacle's ContactGeometry and contact hint.
    obstacle.connectSocket_contact_geometry(contactGeometry);
    obstacle.set_contact_hint(contactHint);
}

const ContactGeometry& Scholz2015GeometryPath::getContactGeometry(
        int obstacleIndex) const {
    const Scholz2015GeometryPathObstacle* obstacle = getObstacle(obstacleIndex);
    return obstacle->getContactGeometry();
}

const SimTK::Vec3& Scholz2015GeometryPath::getContactHint(
        int obstacleIndex) const {
    const Scholz2015GeometryPathObstacle* obstacle = getObstacle(obstacleIndex);
    return obstacle->getContactHint();
}

int Scholz2015GeometryPath::getNumPathPoints() const {
    int count = 0;
    for (int i = 0; i < getNumPathElements(); ++i) {
        if (tryGetElement<Scholz2015GeometryPathAbstractPoint>(i) != nullptr) {
            ++count;
        }
    }
    return count;
}

int Scholz2015GeometryPath::getNumObstacles() const {
    int count = 0;
    for (int i = 0; i < getNumPathElements(); ++i) {
        if (tryGetElement<Scholz2015GeometryPathObstacle>(i) != nullptr) {
            ++count;
        }
    }
    return count;
}

int Scholz2015GeometryPath::getNumPathElements() const {
    return getProperty_path_elements().size();
}

void Scholz2015GeometryPath::setUseWarmStart(bool useWarmStart) {
    if (_index.isValid()) {
        updCableSpan().setUseWarmStart(useWarmStart);
    }
    set_use_warm_start(useWarmStart);
}

bool Scholz2015GeometryPath::getUseWarmStart() const {
    return get_use_warm_start();
}

//=============================================================================
// ABSTRACT PATH INTERFACE
//=============================================================================
double Scholz2015GeometryPath::getLength(const SimTK::State& s) const {
    return getCableSpan().calcLength(s);
}

double Scholz2015GeometryPath::getLengtheningSpeed(
        const SimTK::State& s) const {
    return getCableSpan().calcLengthDot(s);
}

double Scholz2015GeometryPath::computeMomentArm(const SimTK::State& s,
        const Coordinate& coord) const {

    if (!_maSolver) {
        const_cast<Self*>(this)->_maSolver.reset(
            new MomentArmSolver(getModel()));
    }

    return _maSolver->solve(s, coord,  *this);
}

bool Scholz2015GeometryPath::isVisualPath() const {
    return true;
}

std::vector<ComponentPath>
Scholz2015GeometryPath::
findIndependentCoordinates(const SimTK::State& s) const {
    const PhysicalFrame& originFrame = getOrigin().getParentFrame();
    const PhysicalFrame& insertionFrame = getInsertion().getParentFrame();

    std::vector<SimTK::ReferencePtr<const Joint>>
    jointsBetweenFrames = findJointsBetweenPhysicalFrames(
            getModel(),
            originFrame.getAbsolutePathString(),
            insertionFrame.getAbsolutePathString());

    std::vector<ComponentPath> coordinates;
    for (const auto& joint : jointsBetweenFrames) {
        for (int i = 0; i < joint->numCoordinates(); ++i) {
            const Coordinate& coord = joint->get_coordinates(i);
            if (!coord.isConstrained(s)) {
                coordinates.push_back(
                    ComponentPath(coord.getAbsolutePathString()));
            }
        }
    }

    // Moving path points also depend on the coordinates of their location
    // functions.
    for (int ielt = 0; ielt < getNumPathElements(); ++ielt) {
        const auto* point =
                tryGetElement<Scholz2015GeometryPathMovingPoint>(ielt);
        if (!point) {
            continue;
        }
        for (const Coordinate* coord : point->getCoordinates()) {
            ComponentPath path(coord->getAbsolutePathString());
            if (!coord->isConstrained(s) &&
                    std::find(coordinates.begin(), coordinates.end(), path) ==
                            coordinates.end()) {
                coordinates.push_back(path);
            }
        }
    }
    return coordinates;
}

//=============================================================================
// FORCE PRODUCER INTERFACE
//=============================================================================
void Scholz2015GeometryPath::implForEachDecorativePathPoint(
    const SimTK::State& state,
    const ModelDisplayHints& hints,
    const std::function<void(const DecorativePathPoint&)>& callback) const
{
    const auto sink = [&callback](SimTK::Vec3 p) {
        callback(DecorativePathPoint{p});
    };

    if (hints.get_discretize_path()) {
        getCableSpan().calcResampledDecorativePathPoints(
            state,
            hints.get_num_samples_per_wrap_segment(),
            sink
        );
    } else {
        getCableSpan().calcDecorativePathPoints(state, sink);
    }
}

void Scholz2015GeometryPath::produceForces(const SimTK::State& state,
        double tension, ForceConsumer& forceConsumer) const {

    if (tension < 0.) {
        return;
    }

    const SimTK::CableSpan& cable = getCableSpan();
    const auto& origin = getOrigin();
    const auto& insertion = getInsertion();

    // Force applied at path origin point.
    // CableSpan::calcOriginTangentDirection() returns the tangent direction of
    // the cable at origin point, pointing away from the origin in the direction
    // of the tension force.
    forceConsumer.consumePointForce(state,
            origin.getParentFrame(),
            origin.calcLocation(state),
            tension * cable.calcOriginTangentDirection(state));


    // Forces applied to each obstacle body.
    for (const auto& [ix, eltIx] : _obstacleIndexes) {
        if (!cable.isInContactWithObstacle(state, ix)) {
            continue;
        }

        // Compute the initial and final Frenet frame for the curve segment.
        const SimTK::Transform X_GP =
                cable.calcCurveSegmentInitialFrenetFrame(state, ix);
        const SimTK::Transform X_GQ =
                cable.calcCurveSegmentFinalFrenetFrame(state, ix);

        // The tangent direction of the curve segment at the obstacle contact
        // point is the x-axis of the Frenet frame. The force applied to the
        // obstacle at the initial contact point is in the negative direction
        // of the tangent.
        const SimTK::UnitVec3& t_P = -X_GP.R().getAxisUnitVec(SimTK::XAxis);
        const SimTK::UnitVec3& t_Q = X_GQ.R().getAxisUnitVec(SimTK::XAxis);

        // Transform the Frenet frame positions from the ground frame G to the
        // obstacle frame O, to conform with the ForceConsumer interface.
        const Ground& ground = getModel().getGround();
        const auto& obstacle = getElement<Scholz2015GeometryPathObstacle>(eltIx);
        const auto& frame = obstacle.getContactGeometry().getFrame();
        const SimTK::Vec3& p_GP = X_GP.p();
        const SimTK::Vec3& p_GQ = X_GQ.p();
        const SimTK::Vec3 p_OP = ground.findStationLocationInAnotherFrame(
                state, p_GP, frame);
        const SimTK::Vec3 p_OQ = ground.findStationLocationInAnotherFrame(
                state, p_GQ, frame);

        // Consume the forces applied at the initial and final Frenet frames
        // associated with the obstacle.
        forceConsumer.consumePointForce(state, frame, p_OP, tension * t_P);
        forceConsumer.consumePointForce(state, frame, p_OQ, tension * t_Q);
    }

    // Forces applied to each via point.
    for (const auto& [ix, eltIx] : _viaPointIndexes) {
        // Calculate the via point's initial tangent direction (point from the
        // via point torwards the origin) and the final tangent direction
        // (point from the via point towards the insertion).
        const SimTK::UnitVec3 t_P =
                -cable.calcViaPointIncomingTangentDirection(state, ix);
        const SimTK::UnitVec3 t_Q =
                cable.calcViaPointOutgoingTangentDirection(state, ix);

        // Consume the forces applied at the via point.
        const auto& point =
                getElement<Scholz2015GeometryPathAbstractPoint>(eltIx);
        const auto& frame = point.getParentFrame();
        const SimTK::Vec3 location = point.calcLocation(state);
        forceConsumer.consumePointForce(state, frame, location, tension * t_P);
        forceConsumer.consumePointForce(state, frame, location, tension * t_Q);
    }

    // Force applied at path insertion point.
    // CableSpan::calcTerminationTangentDirection() returns the tangent
    // direction of the cable at insertion point, pointing towards the insertion
    // in the opposite direction of the tension force.
    forceConsumer.consumePointForce(state,
            insertion.getParentFrame(),
            insertion.calcLocation(state),
            -tension * cable.calcTerminationTangentDirection(state));

    // Moving path points move relative to their frames, so the tension also
    // does work through the coordinates they depend on.
    if (!_movingPointCoordinates.empty()) {
        const SimTK::SimbodyMatterSubsystem& matter =
                getModel().getMatterSubsystem();
        SimTK::Vector mobilityForces(state.getNU(), 0.);
        cable.applyMobilityForces(state, tension, mobilityForces);
        for (const Coordinate* coord : _movingPointCoordinates) {
            // Coordinates in OpenSim have qdot = u.
            const SimTK::UIndex ux(
                    matter.getMobilizedBody(coord->getBodyIndex())
                            .getFirstUIndex(state) +
                    coord->getMobilizerQIndex());
            forceConsumer.consumeGeneralizedForce(state, *coord,
                    mobilityForces[ux]);
        }
    }
}

//=============================================================================
// MODEL COMPONENT INTERFACE
//=============================================================================
void Scholz2015GeometryPath::extendConnectToModel(Model& model) {
    Super::extendConnectToModel(model);

    OPENSIM_THROW_IF_FRMOBJ(getNumPathPoints() < 2, Exception,
            "Expected at least two path points before finalizing the path, but "
            "{} path point(s) were found. Use addPathPoint() to add a path "
            "point to the path.", getNumPathPoints());
    OPENSIM_THROW_IF_FRMOBJ(
            tryGetElement<Scholz2015GeometryPathAbstractPoint>(0) == nullptr,
            Exception, "Expected the first element of the path to be a path "
            "point, but it is an obstacle.");
    OPENSIM_THROW_IF_FRMOBJ(
            tryGetElement<Scholz2015GeometryPathAbstractPoint>(
                    getNumPathElements() - 1) == nullptr,
            Exception, "Expected the last element of the path to be a path "
            "point, but it is an obstacle.");
}

void Scholz2015GeometryPath::extendAddToSystem(
        SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);

    // First, we need to clear the via point and obstacle indexes so that
    // repeated calls to addToSystem() do not add the same via points and
    // obstacles to the list of indexes multiple times.
    _viaPointIndexes.clear();
    _obstacleIndexes.clear();
    _movingPointCoordinates.clear();

    // Moving path points: the location function in the mobilized body frame
    // and the mobilized body and mobilizer q index of each coordinate.
    struct MovingStation {
        std::shared_ptr<const SimTK::Function_<SimTK::Vec3>> function;
        SimTK::Array_<SimTK::MobilizedBodyIndex> coordBodies;
        SimTK::Array_<SimTK::MobilizerQIndex> coordIndices;
    };
    auto createMovingStation =
            [&](const Scholz2015GeometryPathMovingPoint& point) {
        MovingStation station;
        station.function = std::make_shared<MovingPointFunction>(
                point.getParentFrame().findTransformInBaseFrame(),
                point.get_x_location(), point.get_y_location(),
                point.get_z_location());
        for (const Coordinate* coord : point.getCoordinates()) {
            station.coordBodies.push_back(coord->getBodyIndex());
            station.coordIndices.push_back(coord->getMobilizerQIndex());
            if (std::find(_movingPointCoordinates.begin(),
                        _movingPointCoordinates.end(), coord) ==
                    _movingPointCoordinates.end()) {
                _movingPointCoordinates.push_back(coord);
            }
        }
        return station;
    };

    // The station of a fixed origin or insertion; moving ones are set below.
    auto fixedStation = [&](int index) {
        const auto* point = tryGetElement<Scholz2015GeometryPathPoint>(index);
        return point ? point->getPathPoint().get_location() : SimTK::Vec3(0);
    };
    const int insertionIx = getNumPathElements() - 1;
    SimTK::CableSubsystem& cables = system.updCableSubsystem();
    SimTK::CableSpan cable(cables,
        getOrigin().getParentFrame().getMobilizedBodyIndex(),
        fixedStation(0),
        getInsertion().getParentFrame().getMobilizedBodyIndex(),
        fixedStation(insertionIx));

    if (const auto* point =
                tryGetElement<Scholz2015GeometryPathMovingPoint>(0)) {
        MovingStation station = createMovingStation(*point);
        cable.setOriginStation(station.function, station.coordBodies,
                station.coordIndices);
    }
    if (const auto* point =
                tryGetElement<Scholz2015GeometryPathMovingPoint>(insertionIx)) {
        MovingStation station = createMovingStation(*point);
        cable.setTerminationStation(station.function, station.coordBodies,
                station.coordIndices);
    }

    for (int ielt = 1; ielt < getProperty_path_elements().size() - 1; ++ielt) {
        if (const auto* point =
                    tryGetElement<Scholz2015GeometryPathPoint>(ielt)) {
            SimTK::CableSpanViaPointIndex ix = cable.addViaPoint(
                point->getPathPoint().getParentFrame().getMobilizedBodyIndex(),
                point->getPathPoint().get_location());
            _viaPointIndexes.emplace_back(std::make_pair(ix, ielt));
        } else if (const auto* point =
                    tryGetElement<Scholz2015GeometryPathMovingPoint>(ielt)) {
            MovingStation station = createMovingStation(*point);
            SimTK::CableSpanViaPointIndex ix = cable.addViaPoint(
                point->getParentFrame().getMobilizedBodyIndex(),
                station.function, station.coordBodies, station.coordIndices);
            _viaPointIndexes.emplace_back(std::make_pair(ix, ielt));
        } else if (const auto* obstacle =
                    tryGetElement<Scholz2015GeometryPathObstacle>(ielt)) {
            SimTK::CableSpanObstacleIndex ix = cable.addObstacle(
                obstacle->getContactGeometry().getFrame()
                                              .getMobilizedBodyIndex(),
                obstacle->getContactGeometry().getTransform(),
                std::make_shared<SimTK::ContactGeometry>(
                    obstacle->getContactGeometry().createSimTKContactGeometry()),
                obstacle->getContactHint());
            _obstacleIndexes.emplace_back(std::make_tuple(ix, ielt));
        }
    }

    cable.setSmoothnessTolerance(1e-5);
    cable.setCurveSegmentAccuracy(1e-10);
    cable.setSolverMaxIterations(50);
    cable.setAlgorithm(SimTK::CableSpanAlgorithm::Scholz2015);
    cable.setUseWarmStart(get_use_warm_start());
    _index = cable.getIndex();
}


void Scholz2015GeometryPath::extendPreScale(const SimTK::State& s,
        const ScaleSet& scaleSet) {
    Super::extendPreScale(s, scaleSet);
    setPreScaleLength(s, getLength(s));
}

void Scholz2015GeometryPath::extendPostScale(const SimTK::State& s,
        const ScaleSet& scaleSet) {
    Super::extendPostScale(s, scaleSet);
    getLength(s);
}

//=============================================================================
// CONVENIENCE METHODS
//=============================================================================
void Scholz2015GeometryPath::constructProperties() {
    constructProperty_path_elements();
    constructProperty_use_warm_start(false);
}

const Scholz2015GeometryPathObstacle* Scholz2015GeometryPath::getObstacle(
        int obstacleIndex) const {
    OPENSIM_THROW_IF_FRMOBJ(
            obstacleIndex < 0 || obstacleIndex >= getNumObstacles(), Exception,
            "Index {} is out of range. There are {} obstacles in the path.",
            obstacleIndex, getNumObstacles());

    int count = 0;
    for (int i = 0; i < getNumPathElements(); ++i) {
        if (const auto* obstacle =
                tryGetElement<Scholz2015GeometryPathObstacle>(i)) {
            if (count == obstacleIndex) {
                return obstacle;
            }
            ++count;
        }
    }

    // We should never reach this point.
    OPENSIM_THROW_FRMOBJ(Exception,
            "Obstacle index {} could not be found.", obstacleIndex);
}

const SimTK::CableSpan& Scholz2015GeometryPath::getCableSpan() const {
    return getModel().getMultibodySystem().getCableSubsystem().getCable(_index);
}

SimTK::CableSpan& Scholz2015GeometryPath::updCableSpan() {
    return getModel().updMultibodySystem().updCableSubsystem().updCable(_index);
}
