extern crate affine_transforms;

use affine_transforms::matrices::{ AffineMatrix };
use std::f64::consts::{ PI };

struct Body {
    /// kilograms
    mass : f64, 
    /// meters
    radius : f64
}

impl Body {
    /// Units m^3/s^2 (Mass * G)
    fn k(&self, G : f64) -> f64 { 
        self.mass * G
    }
}

enum Orbit {
    Circular {
        angle_of_periapsis : f64,    
        ascending_node : f64,
        inclination : f64,
        semimajor_axis : f64,        
    },
    Elliptical {
        angle_of_periapsis : f64,    
        ascending_node : f64,
        inclination : f64,
        periapsis : f64,  
        eccentricity : f64,
    },
    Parabolic { 
        angle_of_periapsis : f64,    
        ascending_node : f64,
        inclination : f64,
        periapsis : f64,   
    },
    Hyperbolic { 
        angle_of_periapsis : f64,    
        ascending_node : f64,
        inclination : f64,
        periapsis : f64,
        eccentricity : f64,
    }
}


impl Orbit {
    
    fn eccentricity(&self) -> f64 {
        match *self {
            Orbit::Circular { .. }=> 0.0, 
            Orbit::Elliptical { eccentricity, .. } => eccentricity, 
            Orbit::Parabolic { .. } => 1.0,    
            Orbit::Hyperbolic { eccentricity, .. } => eccentricity  
        }
    }

    // fn get_angle_of_descending_node(&self) -> f64 { *self.ascending_node + PI }            
    // fn get_major_axis(&self) -> f64 { *self.periapsis + *self.apoapsis }
    // fn get_semimajor_axis(&self) -> f64 { self.get_major_axis() / 2 }

    // pub fn next_anomaly(&self, start : &Anomaly, time_delta_ms : f64) -> Anomaly {
    //     panic!("Not implememented");
    // }
    
    // pub fn reference_frame_at_anomaly(&self, anomaly : &Anomaly) -> AffineMatrix {
    //     panic!("Not implememented");
    // }

    // /// Distance moved 
    // pub fn velocity_at_anomaly(&self, anomaly : &Anomaly) -> f64 {
    //     panic!("Not implememented");
    // }
}

struct Anomaly {
    time_ms : u64,
    true_anomaly : f64,
    mean_anomaly : f64,    
    eccentric_anomaly : f64
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
