var geoConverter = geoConverter || {};

(function (GEO) {
    "use strict";

    var /* constants */
        pi = 3.14159265358979,
        /* Ellipsoid model constants (actual values here are for WGS84): */
        sm_a = 6378137.0,
        sm_b = 6356752.314,
        sm_EccSquared = 6.69437999013e-03,
        UTMScaleFactor = 0.9996,

        /* public functions */
        degToRad,
        radToDeg,
        arcLengthOfMeridian,
        footpointLatitude,
        utmCentralMeridian,
        mapLatLonToXY,
        mapXYToLatLon,
        latLonToUTMXYZone,
        utmXYZoneToLatLon;

    /*
    * degToRad
    *
    * Converts degrees to radians.
    *
    */
    degToRad = function (deg) {
        return (deg / 180.0 * pi);
    };

    /*
    * radToDeg
    *
    * Converts radians to degrees.
    *
    */
    radToDeg = function (rad) {
        return (rad / pi * 180.0);
    };

    /*
    * arcLengthOfMeridian
    *
    * Computes the ellipsoidal distance from the equator to a point at a
    * given latitude.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    * GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *     phi - Latitude of the point, in radians.
    *
    * Globals:
    *     sm_a - Ellipsoid model major axis.
    *     sm_b - Ellipsoid model minor axis.
    *
    * Returns:
    *     The ellipsoidal distance of the point from the equator, in meters.
    *
    */
    arcLengthOfMeridian = function (phi) {
        var alpha, beta, gamma, delta, epsilon, n;

        /* Precalculate n */
        n = (sm_a - sm_b) / (sm_a + sm_b);

        /* Precalculate alpha */
        alpha = ((sm_a + sm_b) / 2.0) * (1.0 + (Math.pow(n, 2.0) / 4.0) + (Math.pow(n, 4.0) / 64.0));

        /* Precalculate beta */
        beta = (-3.0 * n / 2.0) + (9.0 * Math.pow(n, 3.0) / 16.0) + (-3.0 * Math.pow(n, 5.0) / 32.0);

        /* Precalculate gamma */
        gamma = (15.0 * Math.pow(n, 2.0) / 16.0) + (-15.0 * Math.pow(n, 4.0) / 32.0);

        /* Precalculate delta */
        delta = (-35.0 * Math.pow(n, 3.0) / 48.0) + (105.0 * Math.pow(n, 5.0) / 256.0);

        /* Precalculate epsilon */
        epsilon = (315.0 * Math.pow(n, 4.0) / 512.0);

        /* Now calculate the sum of the series and return */
        return alpha * (phi + (beta * Math.sin(2.0 * phi)) + (gamma * Math.sin(4.0 * phi)) + (delta * Math.sin(6.0 * phi)) + (epsilon * Math.sin(8.0 * phi)));
    };

    /*
    * utmCentralMeridian
    *
    * Determines the central meridian for the given UTM zone.
    *
    * Inputs:
    *     zone - An integer value designating the UTM zone, range [1,60].
    *
    * Returns:
    *   The central meridian for the given UTM zone, in radians, or zero
    *   if the UTM zone parameter is outside the range [1,60].
    *   Range of the central meridian is the radian equivalent of [-177,+177].
    *
    */
    utmCentralMeridian = function (zone) {
        return degToRad(-183.0 + (zone * 6.0));
    };

    /*
    * footpointLatitude
    *
    * Computes the footpoint latitude for use in converting transverse
    * Mercator coordinates to ellipsoidal coordinates.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    *   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *   y - The UTM northing coordinate, in meters.
    *
    * Returns:
    *   The footpoint latitude, in radians.
    *
    */
    footpointLatitude = function (y) {
        var yy, alpha, beta, gamma, delta, epsilon, n;

        /* Precalculate n (Eq. 10.18) */
        n = (sm_a - sm_b) / (sm_a + sm_b);

        /* Precalculate alpha (Eq. 10.22) (Same as alpha in Eq. 10.17) */
        alpha = ((sm_a + sm_b) / 2.0) * (1 + (Math.pow(n, 2.0) / 4) + (Math.pow(n, 4.0) / 64));

        /* Precalculate y (Eq. 10.23) */
        yy = y / alpha;

        /* Precalculate beta (Eq. 10.22) */
        beta = (3.0 * n / 2.0) + (-27.0 * Math.pow(n, 3.0) / 32.0) + (269.0 * Math.pow(n, 5.0) / 512.0);

        /* Precalculate gamma (Eq. 10.22) */
        gamma = (21.0 * Math.pow(n, 2.0) / 16.0) + (-55.0 * Math.pow(n, 4.0) / 32.0);

        /* Precalculate delta (Eq. 10.22) */
        delta = (151.0 * Math.pow(n, 3.0) / 96.0) + (-417.0 * Math.pow(n, 5.0) / 128.0);

        /* Precalculate epsilon (Eq. 10.22) */
        epsilon = (1097.0 * Math.pow(n, 4.0) / 512.0);

        /* Now calculate the sum of the series (Eq. 10.21) */
        return yy + (beta    * Math.sin(2.0 * yy))
                  + (gamma   * Math.sin(4.0 * yy))
                  + (delta   * Math.sin(6.0 * yy))
                  + (epsilon * Math.sin(8.0 * yy));
    };

    /*
    * mapLatLonToXY
    *
    * Converts a latitude/longitude pair to x and y coordinates in the
    * Transverse Mercator projection.  Note that Transverse Mercator is not
    * the same as UTM; a scale factor is required to convert between them.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    * GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *    phi - Latitude of the point, in radians.
    *    lambda - Longitude of the point, in radians.
    *    lambda0 - Longitude of the central meridian to be used, in radians.
    *
    * Returns:
    *    An object containing the x and y coordinates of the computed point.
    *
    */
    mapLatLonToXY = function (phi, lambda, lambda0) {
        var N, nu2, ep2, t, t2, t4, t6, l,
            l3coef, l4coef, l5coef, l6coef, l7coef, l8coef,
            tmp, x, y;

        /* Precalculate ep2 */
        ep2 = (Math.pow(sm_a, 2.0) - Math.pow(sm_b, 2.0)) / Math.pow(sm_b, 2.0);

        /* Precalculate nu2 */
        nu2 = ep2 * Math.pow(Math.cos(phi), 2.0);

        /* Precalculate N */
        N = Math.pow(sm_a, 2.0) / (sm_b * Math.sqrt(1 + nu2));

        /* Precalculate t */
        t = Math.tan(phi);
        t2 = t * t;
        t4 = t2 * t2;
        t6 = t4 * t2;
        tmp = t6 - Math.pow(t, 6.0);

        /* Precalculate l */
        l = lambda - lambda0;

        /* Precalculate coefficients for l**n in the equations below
           so a normal human being can read the expressions for easting
           and northing
           -- l**1 and l**2 have coefficients of 1.0 */
        l3coef = 1.0 - t2 + nu2;
        l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
        l5coef = 5.0 - 18.0 * t2 + t4 + 14.0 * nu2 - 58.0 * t2 * nu2;
        l6coef = 61.0 - 58.0 * t2 + t4 + 270.0 * nu2 - 330.0 * t2 * nu2;
        l7coef = 61.0 - 479.0 * t2 + 179.0 * t4 - t6;
        l8coef = 1385.0 - 3111.0 * t2 + 543.0 * t4 - t6;

        /* Calculate easting (x) */
        x = N * Math.cos(phi) * l
            + (N / 6.0 * Math.pow(Math.cos(phi), 3.0) * l3coef * Math.pow(l, 3.0))
            + (N / 120.0 * Math.pow(Math.cos(phi), 5.0) * l5coef * Math.pow(l, 5.0))
            + (N / 5040.0 * Math.pow(Math.cos(phi), 7.0) * l7coef * Math.pow(l, 7.0));

        /* Calculate northing (y) */
        y = arcLengthOfMeridian(phi)
            + (t / 2.0 * N * Math.pow(Math.cos(phi), 2.0) * Math.pow(l, 2.0))
            + (t / 24.0 * N * Math.pow(Math.cos(phi), 4.0) * l4coef * Math.pow(l, 4.0))
            + (t / 720.0 * N * Math.pow(Math.cos(phi), 6.0) * l6coef * Math.pow(l, 6.0))
            + (t / 40320.0 * N * Math.pow(Math.cos(phi), 8.0) * l8coef * Math.pow(l, 8.0));

        return {'x': x, 'y': y};
    };

    /*
    * mapXYToLatLon
    *
    * Converts x and y coordinates in the Transverse Mercator projection to
    * a latitude/longitude pair.  Note that Transverse Mercator is not
    * the same as UTM; a scale factor is required to convert between them.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    *   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *   x - The easting of the point, in meters.
    *   y - The northing of the point, in meters.
    *   lambda0 - Longitude of the central meridian to be used, in radians.
    *
    * Returns:
    *   An object containing:
    *   lat - the latitude in radians
    *   lon - the longitude in radians
    *
    * Remarks:
    *   The local variables Nf, nuf2, tf, and tf2 serve the same purpose as
    *   N, nu2, t, and t2 in MapLatLonToXY, but they are computed with respect
    *   to the footpoint latitude phif.
    *
    *   x1frac, x2frac, x2poly, x3poly, etc. are to enhance readability and
    *   to optimize computations.
    *
    */
    mapXYToLatLon = function (x, y, lambda0) {
        var phif, Nf, Nfpow, nuf2, ep2, tf, tf2, tf4, cf,
            x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac,
            x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly,
            philambda = {};

        /* Get the value of phif, the footpoint latitude. */
        phif = footpointLatitude(y);

        /* Precalculate ep2 */
        ep2 = (Math.pow(sm_a, 2.0) - Math.pow(sm_b, 2.0)) / Math.pow(sm_b, 2.0);

        /* Precalculate cos (phif) */
        cf = Math.cos(phif);

        /* Precalculate nuf2 */
        nuf2 = ep2 * Math.pow(cf, 2.0);

        /* Precalculate Nf and initialize Nfpow */
        Nf = Math.pow(sm_a, 2.0) / (sm_b * Math.sqrt(1 + nuf2));
        Nfpow = Nf;

        /* Precalculate tf */
        tf = Math.tan(phif);
        tf2 = tf * tf;
        tf4 = tf2 * tf2;

        /* Precalculate fractional coefficients for x**n in the equations
           below to simplify the expressions for latitude and longitude. */
        x1frac = 1.0 / (Nfpow * cf);

        Nfpow *= Nf;   /* now equals Nf**2) */
        x2frac = tf / (2.0 * Nfpow);

        Nfpow *= Nf;   /* now equals Nf**3) */
        x3frac = 1.0 / (6.0 * Nfpow * cf);

        Nfpow *= Nf;   /* now equals Nf**4) */
        x4frac = tf / (24.0 * Nfpow);

        Nfpow *= Nf;   /* now equals Nf**5) */
        x5frac = 1.0 / (120.0 * Nfpow * cf);

        Nfpow *= Nf;   /* now equals Nf**6) */
        x6frac = tf / (720.0 * Nfpow);

        Nfpow *= Nf;   /* now equals Nf**7) */
        x7frac = 1.0 / (5040.0 * Nfpow * cf);

        Nfpow *= Nf;   /* now equals Nf**8) */
        x8frac = tf / (40320.0 * Nfpow);

        /* Precalculate polynomial coefficients for x**n.
           -- x**1 does not have a polynomial coefficient. */
        x2poly = -1.0 - nuf2;
        x3poly = -1.0 - 2 * tf2 - nuf2;
        x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2 - 3.0 * (nuf2 * nuf2) - 9.0 * tf2 * (nuf2 * nuf2);
        x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;
        x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2;
        x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);
        x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);

        /* Calculate latitude */
        philambda.lat = phif + x2frac * x2poly * (x * x)
            + x4frac * x4poly * Math.pow(x, 4.0)
            + x6frac * x6poly * Math.pow(x, 6.0)
            + x8frac * x8poly * Math.pow(x, 8.0);

        /* Calculate longitude */
        philambda.lon = lambda0 + x1frac * x
            + x3frac * x3poly * Math.pow(x, 3.0)
            + x5frac * x5poly * Math.pow(x, 5.0)
            + x7frac * x7poly * Math.pow(x, 7.0);

        return philambda;
    };

    /*
    * latLonToUTMXY
    *
    * Converts a latitude/longitude pair to x and y coordinates in the
    * Universal Transverse Mercator projection.
    *
    * Inputs:
    *   lat       - Latitude of the point, in radians.
    *   lon       - Longitude of the point, in radians.
    *   zoneParam - UTM zone to be used for calculating values for x and y.
    *               If zoneParam is less than 1 or greater than 60, or if zoneParam is not provided, the routine
    *               will determine the appropriate zone from the value of lon.
    *
    * Returns:
    *   An object containing three elements:
    *   x, y - the x and y coordinates
    *   zone - The UTM zone used for calculating the values of x and y.
    *
    */
    latLonToUTMXYZone = function (degLat, degLon, zoneParam) {
        var zone, utmXY, lat, lon;

		lat = degToRad(degLat);
		lon = degToRad(degLon);
		if (zoneParam && typeof zoneParam === 'number' && zoneParam >= 1 && zoneParam <= 60) {
			zone = zoneParam;
		} else {
	        // Compute the UTM zone.
            zone = Math.floor((degLon + 180.0) / 6) + 1;
        }

        utmXY = mapLatLonToXY(lat, lon, utmCentralMeridian(zone));

        /* Adjust easting and northing for UTM system. */
        utmXY.x = utmXY.x * UTMScaleFactor + 500000.0;
        utmXY.y = utmXY.y * UTMScaleFactor;
        if (utmXY.y < 0.0) {
            utmXY.y = utmXY.y + 10000000.0;
		}
        return {'x': utmXY.x, 'y': utmXY.y, 'zone': zone};
    };

    /*
    * utmXYZoneToLatLon
    *
    * Converts x and y coordinates in the Universal Transverse Mercator
    * projection to a latitude/longitude pair.
    *
    * Inputs:
    *    x - The easting of the point, in meters.
    *    y - The northing of the point, in meters.
    *    zone - The UTM zone in which the point lies.
    *    southhemi - True if the point is in the southern hemisphere;
    *               false otherwise.
    *
    * Outputs:
    *    latlon - A 2-element array containing the latitude and
    *            longitude of the point, in radians.
    *
    * Returns:
    *    The function does not return a value.
    *
    */
    utmXYZoneToLatLon = function (x, y, zone, southhemi) {
        var cmeridian, latlon;

        x -= 500000.0;
        x /= UTMScaleFactor;

        /* If in southern hemisphere, adjust y accordingly. */
        if (southhemi) {
            y -= 10000000.0;
        }
        y /= UTMScaleFactor;

        cmeridian = utmCentralMeridian(zone);
        latlon = mapXYToLatLon(x, y, cmeridian);
		latlon.lat = radToDeg(latlon.lat);
		latlon.lon = radToDeg(latlon.lon);
	
        return latlon;
    };

	GEO.latLonToUTMXYZone = latLonToUTMXYZone;
    GEO.utmXYZoneToLatLon = utmXYZoneToLatLon;
})(geoConverter);
