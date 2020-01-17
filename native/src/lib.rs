extern crate affine_transforms;

use affine_transforms::matrices::{ AffineMatrix };
use std::f64::consts::{ PI };

struct Body
{
    /// kilograms
    mass : f64,
    /// meters
    radius : f64,
    k : f64,
    G : f64,
}

impl Body
{
    /// Units m^3/s^2 (Mass * G)
    fn k(mass : f64, G : f64) -> f64 { mass * G }
}

impl Default for Body
{
    fn default() -> Self
    {
        let mass = 1.0;
        let radius = 1.0;
        let G = 6.67408e-11;

        Body {
            mass : radius,
            radius : radius,
            G : G,
            k : Body::k(0.0, G)
        }
    }
}

struct KeplerianElements
{    
    semimajor_axis : f64,
    eccentricity : f64,
    inclination : f64,
    ascending_node : f64,
    angle_of_periapsis : f64,
}

enum Orbit {
    Circular(KeplerianElements),
    Elliptical(KeplerianElements),
    Parabolic(KeplerianElements),
    Hyperbolic(KeplerianElements)
}


impl Orbit {

    /* Core Invariants - semi-major axis, eccentricity, inclination, longitude of ascending node, and angle of periapsis */

    fn eccentricity(&self) -> f64
    {
        match self
        {
            Orbit::Circular(_elements) => 0.0,
            Orbit::Elliptical(elements) => elements.eccentricity,
            Orbit::Parabolic(_elements) => 1.0,
            Orbit::Hyperbolic(elements) => elements.eccentricity
        }
    }

    fn semimajor_axis(&self) -> Option<f64>
    {
        match self
        {
            Orbit::Circular(elements) => Some(elements.semimajor_axis),
            Orbit::Elliptical(elements) => Some(elements.semimajor_axis),
            Orbit::Parabolic(_elements) => None,
            Orbit::Hyperbolic(elements) => Some(elements.semimajor_axis)
        }
    }

    fn inclination(&self) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => elements.inclination,
            Orbit::Elliptical(elements) => elements.inclination,
            Orbit::Parabolic(elements) => elements.inclination,
            Orbit::Hyperbolic(elements) => elements.inclination
        }
    }

    fn angle_of_ascending_node(&self) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => elements.ascending_node,
            Orbit::Elliptical(elements) => elements.ascending_node,
            Orbit::Parabolic(elements) => elements.ascending_node,
            Orbit::Hyperbolic(elements) => elements.ascending_node
        }
    }

    fn angle_of_periapsis(&self) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => elements.angle_of_periapsis,
            Orbit::Elliptical(elements) => elements.angle_of_periapsis,
            Orbit::Parabolic(elements) => elements.angle_of_periapsis,
            Orbit::Hyperbolic(elements) => elements.angle_of_periapsis
        }
    }

    /// q
    fn periapsis(&self) -> f64 // q
    {
        match self
        {
            Orbit::Circular(elements) => elements.semimajor_axis,
            Orbit::Elliptical(elements) => elements.semimajor_axis * (1.0 - elements.eccentricity),
            Orbit::Parabolic(elements) => elements.semimajor_axis,
            Orbit::Hyperbolic(elements) => elements.semimajor_axis * (1.0 - elements.eccentricity)
        }
    }

    /// p
    fn parameter(&self) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => elements.semimajor_axis,
            Orbit::Elliptical(elements) => elements.semimajor_axis * (1.0 - elements.eccentricity.powf(2.0)),
            Orbit::Parabolic(elements) => 2.0 * elements.semimajor_axis,
            Orbit::Hyperbolic(elements) => elements.semimajor_axis * (1.0 - elements.eccentricity.powf(2.0))
        }
    }

    /// E
    fn total_energy(&self, body : &Body) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => -body.k / 2.0 * elements.semimajor_axis,
            Orbit::Elliptical(elements) => -body.k / 2.0 * elements.semimajor_axis,
            Orbit::Parabolic(_elements) => 0.0,
            Orbit::Hyperbolic(elements) => -body.k / 2.0 * elements.semimajor_axis
        }
    }

    fn distance_from_parent(&self, body : &Body, anomaly : &Anomaly) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => elements.semimajor_axis,
            Orbit::Elliptical(elements) => self.parameter() / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos()),
            Orbit::Parabolic(_elements) => 2.0 * self.periapsis() / (1.0 + anomaly.true_anomaly.cos()),
            Orbit::Hyperbolic(elements) => self.parameter() / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos())
        }
    }

    fn velocity(&self, body : &Body, anomaly : &Anomaly) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => body.k * (1.0 / elements.semimajor_axis),
            Orbit::Elliptical(elements) => (body.k * (2.0 / self.distance_from_parent(body, anomaly) - 1.0 / elements.semimajor_axis)).sqrt(),
            Orbit::Parabolic(_elements) => (body.k * (2.0 / self.distance_from_parent(body, anomaly))).sqrt(),
            Orbit::Hyperbolic(elements) => (body.k * (2.0 / self.distance_from_parent(body, anomaly) - 1.0 / elements.semimajor_axis)).sqrt()
        }
    }

    /// v
    fn angle_of_velocity(&self, anomaly : &Anomaly) -> f64
    {
        match self
        {
            Orbit::Circular(_elements) => 0.0,
            Orbit::Elliptical(elements) => (elements.eccentricity * anomaly.true_anomaly.sin() / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos())).atan(),
            Orbit::Parabolic(_elements) => anomaly.true_anomaly / 2.0,
            Orbit::Hyperbolic(elements) => (elements.eccentricity * anomaly.true_anomaly.sin() / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos())).atan(),
        }
    }

    /// Vq
    fn velocity_at_periapsis(&self, body : &Body) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => (body.k / elements.semimajor_axis).sqrt(),
            Orbit::Elliptical(elements) => ((body.k / elements.semimajor_axis) * (1.0 + elements.eccentricity) / (1.0 - elements.eccentricity)).sqrt(),
            Orbit::Parabolic(_elements) => (body.k * self.periapsis() / 2.0).sqrt(),
            Orbit::Hyperbolic(elements) => ((body.k / elements.semimajor_axis) * (1.0 + elements.eccentricity) / (elements.eccentricity - 1.0)).sqrt(),
        }
    }

    /// A - rate of area swept by orbit
    fn areal_velocity(&self, body : &Body) -> f64
    {
        match self
        {
            Orbit::Circular(elements) => (body.k * elements.semimajor_axis).sqrt(),
            Orbit::Elliptical(elements) => ((body.k * elements.semimajor_axis) * (1.0 + elements.eccentricity) / (1.0 - elements.eccentricity)).sqrt(),
            Orbit::Parabolic(_elements) => (body.k * self.periapsis() / 2.0).sqrt(),
            Orbit::Hyperbolic(elements) => ((body.k * elements.semimajor_axis) * (1.0 + elements.eccentricity) / (elements.eccentricity - 1.0)).sqrt(),
        }
    }

    /// P
    fn orbital_period(&self, body : &Body) -> Option<f64>
    {
        match self
        {
            Orbit::Circular(elements) => Some(2.0 * PI * (elements.semimajor_axis.powf(3.0) / body.k).sqrt()),
            Orbit::Elliptical(elements) => Some(2.0 * PI * (elements.semimajor_axis.powf(3.0) / body.k).sqrt()),
            Orbit::Parabolic(_elements) => None,
            Orbit::Hyperbolic(_elements) => None,
        }
    }

    fn eccentric_anomaly(&self, anomaly : &Anomaly) -> f64
    {
        match self
        {
            Orbit::Circular(_elements) => anomaly.true_anomaly,
            Orbit::Elliptical(elements) => ( (elements.eccentricity + anomaly.true_anomaly.cos()) / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos()) ).acos(),
            Orbit::Parabolic(_elements) => (anomaly.true_anomaly / 2.0).tan(),
            Orbit::Hyperbolic(elements) => ( (elements.eccentricity + anomaly.true_anomaly.cos()) / (1.0 + elements.eccentricity * anomaly.true_anomaly.cos()) ).acosh(),
        }
    }

    fn mean_anomaly(&self, anomaly : &Anomaly) -> f64
    {
        match self
        {
            Orbit::Circular(_elements) => anomaly.true_anomaly,
            Orbit::Elliptical(elements) => {
                let E = self.eccentric_anomaly(anomaly);
                E - elements.eccentricity * E.sin()
            },
            Orbit::Parabolic(_elements) => {
                let D = self.eccentric_anomaly(anomaly);
                D + D.powf(3.0) / 3.0
            },
            Orbit::Hyperbolic(elements) => {
                let F = self.eccentric_anomaly(anomaly);
                elements.eccentricity * F.sinh() - F
            }
        }
    }

    fn next_anomaly(&self, dT : u64, body : &Body, anomaly : Anomaly) -> Anomaly
    {
        panic!("Not implemented")
    }
}

struct Anomaly {
    time_ms : u64,
    true_anomaly : f64,
    mean_anomaly : f64,
    eccentric_anomaly : f64
}

#[cfg(test)]
mod tests {

}
