import * as THREE from 'three';
import proj4 from 'proj4';
import Coordinates from '../Core/Geographic/Coordinates';

const DEG2RAD = THREE.Math.DEG2RAD;

function quaternionIdentity(coordinates, target = new THREE.Quaternion()) {
    return coordinates ? target.set(0, 0, 0, 1) : quaternionIdentity;
}

/**
 * The OrientationUtils module provides methods to compute the quaternion that
 * model a rotation provided in various conventions, including between different
 * CRS.<br/>
 * CRS-related methods return an FunctionOrQuaternion if the origin is undefined,
 * so that this function can be create once and applied efficiently to many points of origin.
 * Otherwise, the target quaternion is returned at the provided origin coordinates.
 *
 * @example
 * // Compute the rotation around the point of origin from a frame aligned with Lambert93 axes (epsg:2154),
 * // to the geocentric frame (epsg:4978)
 * quat_crs2crs = OrientationUtils.quaternionFromCRSToCRS("EPSG:2154", "EPSG:4978")(origin);
 * // Compute the rotation of a sensor platform defined by its attitude
 * quat_attitude = OrientationUtils.quaternionFromAttitude(attitude);
 * // Compute the rotation from the sensor platform frame to the geocentric frame
 * quat = quat_crs2crs.multiply(quat_attitude);
 *
 * @module OrientationUtils
 */
export default {
    /**
     * @typedef {Object} Attitude
     * Properties are either defined as (omega, phi, kappa) or as (roll, pitch,
     * heading) or all `undefined`.
     *
     * @property {number} omega - angle in degrees
     * @property {number} phi - angle in degrees
     * @property {number} kappa - angle in degrees
     * @property {number} roll - angle in degrees
     * @property {number} pitch - angle in degrees
     * @property {number} heading - angle in degrees
     */

    /**
     * The transform from the platform frame to the local East, North, Up (ENU)
     * frame is `RotationZ(heading).RotationX(pitch).RotationY(roll)`
     *
     * @param {number} [roll=0] - angle in degrees
     * @param {number} [pitch=0] - angle in degrees
     * @param {number} [heading=0] - angle in degrees
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] - output Quaternion
     *
     * @return {THREE.Quaternion} target
     */
    quaternionFromRollPitchHeading(roll = 0, pitch = 0, heading = 0, target = new THREE.Quaternion()) {
        roll *= DEG2RAD;
        pitch *= DEG2RAD;
        heading *= DEG2RAD;
        // return this.setFromEuler(new THREE.Euler(pitch, roll, heading , 'ZXY')).conjugate();
        return target.setFromEuler(new THREE.Euler(-pitch, -roll, -heading, 'YXZ')); // optimized version of above
    },

    /**
     * From
     * [DocMicMac](https://github.com/micmacIGN/Documentation/raw/master/DocMicMac.pdf),
     * the transform from the platform frame to the local East, North, Up (ENU)
     * frame is:
     *
     * ```
     * RotationX(omega).RotationY(phi).RotationZ(kappa).RotationX(PI)
     * RotationX(PI) <=> Quaternion(1,0,0,0) : converts between the 2 conventions for the camera local frame:
     * X right, Y bottom, Z front : convention in photogrammetry and computer vision
     * X right, Y top,    Z back  : convention in webGL, threejs
     * ```
     *
     * @param {number} [omega=0] - angle in degrees
     * @param {number} [phi=0] - angle in degrees
     * @param {number} [kappa=0] - angle in degrees
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     *
     * @return {THREE.Quaternion} target
     */
    quaternionFromOmegaPhiKappa(omega = 0, phi = 0, kappa = 0, target = new THREE.Quaternion()) {
        omega *= DEG2RAD;
        phi *= DEG2RAD;
        kappa *= DEG2RAD;
        target.setFromEuler(new THREE.Euler(omega, phi, kappa, 'XYZ'));
        target.set(target.w, target.z, -target.y, -target.x); // <=> target.multiply(new THREE.Quaternion(1, 0, 0, 0));
        return target;
    },

    /**
     * Set the quaternion according to the rotation from the platform frame to
     * the local frame.
     *
     * @param {Attitude} attitude - Attitude
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     *
     * @return {THREE.Quaternion} target
     */
    quaternionFromAttitude(attitude, target = new THREE.Quaternion()) {
        if ((attitude.roll !== undefined) || (attitude.pitch !== undefined) || (attitude.heading !== undefined)) {
            return this.quaternionFromRollPitchHeading(attitude.roll, attitude.pitch, attitude.heading, target);
        }
        if ((attitude.omega !== undefined) || (attitude.phi !== undefined) || (attitude.kappa !== undefined)) {
            return this.quaternionFromOmegaPhiKappa(attitude.omega, attitude.phi, attitude.kappa, target);
        }
        return target.set(0, 0, 0, 1);
    },

    /**
     * FunctionOrQuaternion is either a THREE.Quaternion or a function that accepts
     * arguments `(coordinates, target)` and returns the quaternion that models a rotation
     * around the point of origin. If target is not provided, a new quaternion is
     * created and returned instead.
     *
     * @typedef {function|THREE.Quaternion} FunctionOrQuaternion
     *
     * @property {Coordinates} coordinates the origin of the local East North Up
     * (ENU) frame
     * @property {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion.
    */

    /**
     * A Projection object models a Coordinate Reference System (CRS).
     * Such an object is usually created with proj4 using `proj4.defs(crs);`
     *
     * @typedef {Object} Projection
     *
     * @property {string} projName
     */

    /**
     * Set the quaternion according to the rotation from the East North Up (ENU)
     * frame to the geocentric frame. The up direction of the ENU frame is
     * provided by the normalized geodetic normal of the provided coordinates
     * (geodeticNormal property).
     *
     * @param {Coordinates} [coordinates] the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionToGeocent(coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionToGeocent()(coordinates, target); }
        const matrix = new THREE.Matrix4();
        const north = new THREE.Vector3();
        const east = new THREE.Vector3();
        return (coordinates, target = new THREE.Quaternion()) => {
            const up = coordinates.geodesicNormal;
            if (up.x == 0 && up.y == 0) {
                return target.set(0, 0, 0, 1);
            }
            // this is an optimized version of matrix.lookAt(up, new THREE.Vector3(0, 0, 0), new THREE.Vector3(0, 0, 1));
            east.set(-up.y, up.x, 0).normalize();
            north.crossVectors(up, east);
            matrix.makeBasis(east, north, up);
            return target.setFromRotationMatrix(matrix);
        };
    },

    /**
     * Set the quaternion according to the rotation from the East North Up (ENU)
     * frame to the geocentric frame. The up direction of the ENU frame is
     * provided by the normalized geodetic normal of the provided coordinates
     * (geodeticNormal property).
     *
     * @param {Coordinates} [coordinates] the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionFromGeocent(coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionFromGeocent()(coordinates, target); }
        const toGeocent = this.quaternionToGeocent();
        return (coordinates, target = new THREE.Quaternion()) => toGeocent(coordinates, target).conjugate();
    },


    /**
     * Set the quaternion to correct for the
     * <a href="https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/alg0060.pdf">meridian convergence</a>
     * of the East North Up (ENU) frame to the Lambert Conformal Conic (LCC) frame.
     * This is a generally small rotation around Z.
     *
     * @param {Projection} proj the lcc projection (may be parsed using proj4)
     * @param {number} proj.lat0 - the latitude of origin
     * @param {number} proj.long0 - the longitude of the central meridian
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
    */
    quaternionFromLCC(proj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionFromLCC(proj)(coordinates, target); }
        const sinlat0 = Math.sin(proj.lat0);
        const axis = new THREE.Vector3().set(0, 0, 1);
        const coord = new Coordinates('EPSG:4326', 0, 0, 0);
        return (coordinates, target = new THREE.Quaternion()) => {
            const long = coordinates.as(coord.crs, coord).longitude * DEG2RAD;
            return target.setFromAxisAngle(axis, sinlat0 * (proj.long0 - long));
        };
    },

    /**
     * Set the quaternion to correct for the
     * <a href="https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/alg0060.pdf">meridian convergence</a>
     * of the East North Up (ENU) frame to the Lambert Conformal Conic (LCC) frame.
     * This is a generally small rotation around Z.
     *
     * @param {Projection} proj the lcc projection (may be parsed using proj4)
     * @param {number} proj.lat0 - the latitude of origin
     * @param {number} proj.long0 - the longitude of the central meridian
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
    */
    quaternionToLCC(proj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionToLCC(proj)(coordinates, target); }
        const fromLCC = this.quaternionFromLCC(proj);
        return (coordinates, target = new THREE.Quaternion()) => fromLCC(coordinates, target).conjugate();
    },

    /**
     * Set the quaternion to correct for the
     * <a href="https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/alg0061.pdf">meridian convergence</a>
     * of the East North Up (ENU) frame to the Transverse Mercator frame.
     * This is a generally small rotation around Z.
     *
     * @param {Projection} proj the tmerc projection (may be parsed using proj4)
     * @param {number} proj.e - the excentricity of the ellipsoid (supersedes
     * {proj.a} and {proj.b})
     * @param {number} proj.a - the semimajor radius of the ellipsoid axis
     * @param {number} proj.b - the semiminor radius of the ellipsoid axis
     * @param {number} proj.long0 - the longitude of the central meridian
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
    */
    quaternionFromTMerc(proj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionFromTMerc(proj)(coordinates, target); }
        const a2 = proj.a * proj.a;
        const b2 = proj.b * proj.b;
        const e2 = proj.e * proj.e;
        const eta0 = proj.e ? (e2 / (1 - e2)) : (a2 / b2 - 1);
        const axis = new THREE.Vector3().set(0, 0, 1);
        const coord = new Coordinates('EPSG:4326', 0, 0, 0);
        return (coordinates, target = new THREE.Quaternion()) => {
            coordinates.as(coord.crs, coord);
            const long = coord.longitude * DEG2RAD;
            const lat = coord.latitude * DEG2RAD;
            const dlong = proj.long0 - long;
            const coslat = Math.cos(lat);
            const sinlat = Math.sin(lat);
            const tanlat = sinlat / coslat;
            const coslat2 = coslat * coslat;
            const dl2 = dlong * dlong * coslat2;
            const eta2 = eta0 * coslat2;
            const gamma = dlong * sinlat * (1 + dl2 / 3 * (1 + 3 * eta2 + 2 * eta2 * eta2) + dl2 * dl2 * (2 - tanlat) / 15);
            return target.setFromAxisAngle(axis, gamma);
        };
    },

    /**
     * Set the quaternion to correct for the
     * <a href="https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/alg0061.pdf">meridian convergence</a>
     * of the East North Up (ENU) frame to the Transverse Mercator frame.
     * This is a generally small rotation around Z.
     *
     * @param {Projection} proj the tmerc projection (may be parsed using proj4)
     * @param {number} proj.e - the excentricity of the ellipsoid (supersedes
     * {proj.a} and {proj.b})
     * @param {number} proj.a - the semimajor radius of the ellipsoid axis
     * @param {number} proj.b - the semiminor radius of the ellipsoid axis
     * @param {number} proj.long0 - the longitude of the central meridian
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
    */
    quaternionToTMerc(proj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionToTMerc(proj)(coordinates, target); }
        const fromTMerc = this.quaternionFromTMerc(proj);
        return (coordinates, target = new THREE.Quaternion()) => fromTMerc(coordinates, target).conjugate();
    },

    /**
     * Returns the identity quaternion (0,0,0,1), as longlat and ENU frame are aligned (?).
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionFromLongLat(coordinates, target = new THREE.Quaternion()) {
        return quaternionIdentity(coordinates, target);
    },

    /**
     * Returns the identity quaternion (0,0,0,1), as longlat and ENU frame are aligned (?).
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionToLongLat(coordinates, target = new THREE.Quaternion()) {
        return quaternionIdentity(coordinates, target);
    },


    /**
     * Warns for an unimplemented projection, sets the quaternion to the
     * identity (0,0,0,1).
     *
     * @param {Projection} proj - the unimplemented projection (may be parsed
     * using proj4)
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionUnimplemented(proj, coordinates, target = new THREE.Quaternion()) {
        console.warn('quaternionToCRS is not implemented for projections of type', proj.projName);
        return quaternionIdentity(coordinates, target);
    },

    /**
     * Compute the quaternion that models the rotation from the local East North
     * Up (ENU) frame to the frame of the given crs.
     *
     * @param {string|Projection} crsOrProj - the CRS of the target frame or its
     * proj4-compatible object.
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionToCRS(crsOrProj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionToCRS(crsOrProj)(coordinates, target); }
        const proj = crsOrProj.projName ? crsOrProj : proj4.defs(crsOrProj);
        switch (proj.projName) {
            case 'geocent': return this.quaternionToGeocent();
            case 'lcc': return this.quaternionToLCC(proj);
            case 'tmerc': return this.quaternionToTMerc(proj);
            case 'longlat': return this.quaternionToLongLat();
            default: return this.quaternionUnimplemented(proj);
        }
    },

    /**
     * Compute the quaternion that models the rotation from the frame of the
     * given crs to the local East North Up (ENU) frame.
     *
     * @param {string|Projection} crsOrProj - the CRS of the source frame or its
     * proj4-compatible object.
     *
     * @param {Coordinates} [coordinates]  coordinates the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
     */
    quaternionFromCRS(crsOrProj, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionFromCRS(crsOrProj)(coordinates, target); }
        const proj = crsOrProj.projName ? crsOrProj : proj4.defs(crsOrProj);
        switch (proj.projName) {
            case 'geocent': return this.quaternionFromGeocent();
            case 'lcc': return this.quaternionFromLCC(proj);
            case 'tmerc': return this.quaternionFromTMerc(proj);
            case 'longlat': return this.quaternionFromLongLat();
            default: return this.quaternionUnimplemented(proj);
        }
    },

    /**
     * Return the function that computes the quaternion that represents a
     * rotation from coordinates in platform frame, defined using an `Attitude`
     * relative to a local `crsIn` frame, to coordinates expressed in the target
     * `crsOut` frame
     *
     * @param {string} crsIn - the CRS of the input frame.
     * @param {string} crsOut - the CRS of the output frame.
     * @param {Coordinates} [coordinates]  coordinates - the origin of the local East North Up
     * (ENU) frame
     * @param {THREE.Quaternion} [target=new THREE.Quaternion()] output Quaternion
     * @return {FunctionOrQuaternion} The target quaternion if coordinates is defined, otherwise, a function to compute it from coordinates.
    */
    quaternionFromCRSToCRS(crsIn, crsOut, coordinates, target = new THREE.Quaternion()) {
        if (coordinates) { return this.quaternionFromCRSToCRS(crsIn, crsOut)(coordinates, target); }
        if (crsIn == crsOut) {
            return (origin, target = new THREE.Quaternion()) => target.set(0, 0, 0, 1);
        }

        // get rotations from the local East/North/Up (ENU) frame to both CRS.
        const fromCrs = this.quaternionFromCRS(crsIn);
        const toCrs = this.quaternionToCRS(crsOut);
        const quat = new THREE.Quaternion();
        return (origin, target = new THREE.Quaternion()) =>
            toCrs(origin, target).multiply(fromCrs(origin, quat));
    },
};
